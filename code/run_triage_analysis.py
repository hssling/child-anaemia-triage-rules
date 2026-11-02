"""Confirmatory analysis for Manuscript 12.

Executed only after PROTOCOL.md was SHA-256 frozen (see
governance/PROTOCOL_FREEZE.txt). Produces every table reported in the paper.

Usage:
    python code/run_triage_analysis.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svydesign import SurveyDesign, auc_from_sorted, weighted_roc  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data" / "analysis_mother_child_pairs_v1.0.parquet"
TABLES = ROOT / "outputs" / "tables"

N_BOOT = 1000
BOOT_SEED = 20260805

# Pre-specified acceptability bar (PROTOCOL.md section 8).
BAR_SENS, BAR_NPV = 0.90, 0.80
BAR_SENS_LCL, BAR_NPV_LCL = 0.85, 0.75

RULE_LABELS = {
    "R1": "Mother anaemic",
    "R2": "Mother moderately or severely anaemic",
    "R3": "Child aged 6-23 months",
    "R4": "Mother anaemic OR child aged 6-23 months",
    "R5": "Mother anaemic OR poorest two wealth quintiles",
    "R6": "Mother anaemic OR child 6-23 months OR poorest two quintiles",
}


# --------------------------------------------------------------------------
def build_rules(d: pd.DataFrame) -> pd.DataFrame:
    r = pd.DataFrame(index=d.index)
    r["R1"] = d["mother_anaemic"]
    r["R2"] = d["mother_mod_sev"]
    r["R3"] = d["child_6_23"]
    r["R4"] = ((d["mother_anaemic"] == 1) | (d["child_6_23"] == 1)).astype(int)
    r["R5"] = ((d["mother_anaemic"] == 1) | (d["poorest2"] == 1)).astype(int)
    r["R6"] = (
        (d["mother_anaemic"] == 1)
        | (d["child_6_23"] == 1)
        | (d["poorest2"] == 1)
    ).astype(int)
    return r


def accuracy(design: SurveyDesign, test: np.ndarray, truth: np.ndarray) -> dict:
    """Design-based diagnostic accuracy for one binary rule."""
    test = np.asarray(test, dtype=float)
    truth = np.asarray(truth, dtype=float)

    tp = test * truth
    fp = test * (1 - truth)
    fn = (1 - test) * truth
    tn = (1 - test) * (1 - truth)
    ones = np.ones_like(truth)

    out = {}
    for name, num, den in [
        ("sensitivity", tp, truth),
        ("specificity", tn, 1 - truth),
        ("ppv", tp, test),
        ("npv", tn, 1 - test),
        ("tested", test, ones),
        ("prevalence", truth, ones),
    ]:
        p, se, lo, hi = design.proportion_ci(num, den)
        out[name] = p
        out[f"{name}_se"] = se
        out[f"{name}_lcl"] = lo
        out[f"{name}_ucl"] = hi

    # Population shares, used for net benefit and the per-100,000 translation.
    W = design.total(ones)
    out["tp_share"] = design.total(tp) / W
    out["fp_share"] = design.total(fp) / W
    out["fn_share"] = design.total(fn) / W
    out["tn_share"] = design.total(tn) / W
    return out


def likelihood_ratios(boot_weights, test, truth) -> dict:
    """LR+ and LR- with bootstrap percentile intervals.

    Ratios of ratios are not linear, so these use the replicate weights rather
    than Taylor linearisation.
    """
    test = np.asarray(test, dtype=float)
    truth = np.asarray(truth, dtype=float)
    tp, fp = test * truth, test * (1 - truth)
    fn, tn = (1 - test) * truth, (1 - test) * (1 - truth)

    def lrs(w):
        sens = (w * tp).sum() / (w * truth).sum()
        spec = (w * tn).sum() / (w * (1 - truth)).sum()
        lr_pos = sens / (1 - spec) if spec < 1 else np.inf
        lr_neg = (1 - sens) / spec if spec > 0 else np.inf
        return lr_pos, lr_neg

    reps = np.array([lrs(w) for w in boot_weights])
    return {
        "lr_pos_lcl": float(np.percentile(reps[:, 0], 2.5)),
        "lr_pos_ucl": float(np.percentile(reps[:, 0], 97.5)),
        "lr_neg_lcl": float(np.percentile(reps[:, 1], 2.5)),
        "lr_neg_ucl": float(np.percentile(reps[:, 1], 97.5)),
    }


def net_benefit(tp_share: float, fp_share: float, thresholds: np.ndarray) -> np.ndarray:
    odds = thresholds / (1 - thresholds)
    return tp_share - fp_share * odds


# --------------------------------------------------------------------------
def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    d = pd.read_parquet(DATA)
    design = SurveyDesign(d["w"], d["psu"], d["stratum"])

    manifest = {
        "manuscript": "12_child_anaemia_triage_rules",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": (ROOT / "governance" / "PROTOCOL_FREEZE.txt")
        .read_text(encoding="utf-8").split("SHA-256: ")[1].split("\n")[0],
        "n_children": int(len(d)),
        "n_mothers": int(d["caseid"].nunique()),
        "n_psu": int(design.n_psu),
        "n_strata": int(design.n_strata),
        "design_df": int(design.df),
        "singleton_strata_collapsed": int(design.n_singletons_collapsed),
        "bootstrap_replicates": N_BOOT,
        "bootstrap_seed": BOOT_SEED,
    }

    truth = d["child_anaemic"].to_numpy(float)
    rules = build_rules(d)

    # Replicate weights are generated once and reused for every quantity that
    # needs them, so all bootstrap intervals share the same resamples.
    print("generating bootstrap replicate weights ...", flush=True)
    boot = list(design.bootstrap_weights(N_BOOT, seed=BOOT_SEED))

    # ------------------------------------------------------- Table 1
    table1 = descriptives(d, design)
    table1.to_csv(TABLES / "table_1_sample_characteristics.csv", index=False)

    # ------------------------------------------------------- Table 2
    rows = []
    for rid, label in RULE_LABELS.items():
        a = accuracy(design, rules[rid].to_numpy(), truth)
        lr = likelihood_ratios(boot, rules[rid].to_numpy(), truth)
        sens, spec = a["sensitivity"], a["specificity"]
        rows.append({
            "rule": rid,
            "definition": label,
            "children_tested_pct": 100 * a["tested"],
            "children_tested_lcl": 100 * a["tested_lcl"],
            "children_tested_ucl": 100 * a["tested_ucl"],
            "sensitivity_pct": 100 * sens,
            "sensitivity_lcl": 100 * a["sensitivity_lcl"],
            "sensitivity_ucl": 100 * a["sensitivity_ucl"],
            "specificity_pct": 100 * spec,
            "specificity_lcl": 100 * a["specificity_lcl"],
            "specificity_ucl": 100 * a["specificity_ucl"],
            "ppv_pct": 100 * a["ppv"],
            "ppv_lcl": 100 * a["ppv_lcl"],
            "ppv_ucl": 100 * a["ppv_ucl"],
            "npv_pct": 100 * a["npv"],
            "npv_lcl": 100 * a["npv_lcl"],
            "npv_ucl": 100 * a["npv_ucl"],
            "lr_pos": sens / (1 - spec),
            "lr_pos_lcl": lr["lr_pos_lcl"],
            "lr_pos_ucl": lr["lr_pos_ucl"],
            "lr_neg": (1 - sens) / spec,
            "lr_neg_lcl": lr["lr_neg_lcl"],
            "lr_neg_ucl": lr["lr_neg_ucl"],
            "meets_bar": bool(
                sens >= BAR_SENS and a["npv"] >= BAR_NPV
                and a["sensitivity_lcl"] >= BAR_SENS_LCL
                and a["npv_lcl"] >= BAR_NPV_LCL
            ),
            "_tp_share": a["tp_share"],
            "_fp_share": a["fp_share"],
            "_fn_share": a["fn_share"],
            "_tn_share": a["tn_share"],
        })
    table2 = pd.DataFrame(rows)
    table2.drop(columns=[c for c in table2.columns if c.startswith("_")]).to_csv(
        TABLES / "table_2_rule_accuracy.csv", index=False
    )

    # ------------------------------------------------------- Table 3
    per = 100_000
    t3 = []
    for r in rows:
        t3.append({
            "rule": r["rule"],
            "definition": r["definition"],
            "tests_performed": round(per * (r["_tp_share"] + r["_fp_share"])),
            "tests_avoided": round(per * (r["_fn_share"] + r["_tn_share"])),
            "anaemic_detected": round(per * r["_tp_share"]),
            "anaemic_missed": round(per * r["_fn_share"]),
            "non_anaemic_tested": round(per * r["_fp_share"]),
            "anaemic_among_untested_pct": 100 * r["_fn_share"] / (r["_fn_share"] + r["_tn_share"]),
            "tests_avoided_per_anaemic_child_missed":
                (r["_fn_share"] + r["_tn_share"]) / r["_fn_share"] if r["_fn_share"] > 0 else np.nan,
        })
    prev = float((d["w"] * truth).sum() / d["w"].sum())
    t3.append({
        "rule": "TEST-ALL", "definition": "Universal testing of all children 6-59 months",
        "tests_performed": per, "tests_avoided": 0,
        "anaemic_detected": round(per * prev), "anaemic_missed": 0,
        "non_anaemic_tested": round(per * (1 - prev)),
        "anaemic_among_untested_pct": np.nan,
        "tests_avoided_per_anaemic_child_missed": np.nan,
    })
    t3.append({
        "rule": "TEST-NONE", "definition": "No child tested",
        "tests_performed": 0, "tests_avoided": per,
        "anaemic_detected": 0, "anaemic_missed": round(per * prev),
        "non_anaemic_tested": 0,
        "anaemic_among_untested_pct": 100 * prev,
        "tests_avoided_per_anaemic_child_missed": 1 / prev,
    })
    pd.DataFrame(t3).to_csv(TABLES / "table_3_per_100000_children.csv", index=False)

    # ------------------------------------------------------- Table 4 (states)
    t4 = []
    for state, g in d.groupby("state", observed=True):
        if len(g) < 100:
            continue
        sd = SurveyDesign(g["w"], g["psu"], g["stratum"])
        a = accuracy(sd, g["mother_anaemic"].to_numpy(float),
                     g["child_anaemic"].to_numpy(float))
        t4.append({
            "state": state,
            "n_pairs": int(len(g)),
            "child_anaemia_pct": 100 * a["prevalence"],
            "child_anaemia_lcl": 100 * a["prevalence_lcl"],
            "child_anaemia_ucl": 100 * a["prevalence_ucl"],
            "sensitivity_pct": 100 * a["sensitivity"],
            "sensitivity_lcl": 100 * a["sensitivity_lcl"],
            "sensitivity_ucl": 100 * a["sensitivity_ucl"],
            "npv_pct": 100 * a["npv"],
            "npv_lcl": 100 * a["npv_lcl"],
            "npv_ucl": 100 * a["npv_ucl"],
            "specificity_pct": 100 * a["specificity"],
            "children_tested_pct": 100 * a["tested"],
            "meets_bar": bool(
                a["sensitivity"] >= BAR_SENS and a["npv"] >= BAR_NPV
                and a["sensitivity_lcl"] >= BAR_SENS_LCL
                and a["npv_lcl"] >= BAR_NPV_LCL
            ),
        })
    table4 = pd.DataFrame(t4).sort_values("child_anaemia_pct")
    table4.to_csv(TABLES / "table_4_state_level_rule_R1.csv", index=False)
    manifest["states_meeting_bar"] = int(table4["meets_bar"].sum())
    manifest["states_evaluated"] = int(len(table4))

    # ------------------------------------------------------- R7: continuous
    roc = roc_analysis(d, design, boot)
    pd.DataFrame(roc["threshold_table"]).to_csv(
        TABLES / "table_5_maternal_hb_thresholds.csv", index=False
    )
    manifest["auc"] = roc["auc"]
    manifest["auc_lcl"] = roc["auc_lcl"]
    manifest["auc_ucl"] = roc["auc_ucl"]
    manifest["best_sens90_threshold"] = roc["best_sens90"]
    np.save(ROOT / "data" / "roc_curve.npy",
            np.vstack([roc["fpr"], roc["tpr"]]))

    # ------------------------------------------------------- decision curve
    thr = np.arange(0.01, 0.505, 0.005)
    dca = {"threshold": thr, "test_all": net_benefit(prev, 1 - prev, thr),
           "test_none": np.zeros_like(thr)}
    for r in rows:
        dca[r["rule"]] = net_benefit(r["_tp_share"], r["_fp_share"], thr)
    dca_df = pd.DataFrame(dca)
    dca_df.to_csv(TABLES / "table_6_decision_curve.csv", index=False)

    # Threshold at which any rule would overtake universal testing.
    beats = {}
    for r in rows:
        better = thr[dca[r["rule"]] > dca["test_all"] + 1e-12]
        beats[r["rule"]] = float(better.min()) if len(better) else None
    manifest["threshold_at_which_rule_beats_test_all"] = beats

    # ------------------------------------------------------- sensitivity
    sens_rows = sensitivity_analyses(d)
    pd.DataFrame(sens_rows).to_csv(
        TABLES / "table_7_sensitivity_analyses.csv", index=False
    )

    # ------------------------------------------------------- headline numbers
    r1 = next(r for r in rows if r["rule"] == "R1")
    manifest["headline"] = {
        "child_anaemia_prevalence_pct": round(100 * prev, 1),
        "maternal_anaemia_prevalence_pct": round(
            100 * float((d["w"] * d["mother_anaemic"]).sum() / d["w"].sum()), 1),
        "hb_correlation": round(float(weighted_corr(
            d["child_hb"], d["mother_hb"], d["w"])), 3),
        "R1_sensitivity_pct": round(r1["sensitivity_pct"], 1),
        "R1_npv_pct": round(r1["npv_pct"], 1),
        "R1_anaemic_among_untested_pct": round(
            100 * r1["_fn_share"] / (r1["_fn_share"] + r1["_tn_share"]), 1),
        "rules_meeting_bar": int(table2["meets_bar"].sum()),
    }

    (ROOT / "outputs" / "ANALYSIS_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, default=float), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, default=float))


# --------------------------------------------------------------------------
def descriptives(d: pd.DataFrame, design: SurveyDesign) -> pd.DataFrame:
    ones = np.ones(len(d))
    rows = []

    def add(group, level, ind):
        p, se, lo, hi = design.proportion_ci(np.asarray(ind, float), ones)
        sub = d[np.asarray(ind, bool)]
        sub_design = SurveyDesign(sub["w"], sub["psu"], sub["stratum"])
        ca, _, ca_lo, ca_hi = sub_design.proportion_ci(
            sub["child_anaemic"].to_numpy(float), np.ones(len(sub))
        )
        rows.append({
            "characteristic": group, "level": level,
            "n_unweighted": int(len(sub)),
            "weighted_pct_of_sample": 100 * p,
            "child_anaemia_pct": 100 * ca,
            "child_anaemia_lcl": 100 * ca_lo,
            "child_anaemia_ucl": 100 * ca_hi,
        })

    add("All children", "6-59 months", ones)
    for lo, hi, lab in [(6, 11, "6-11"), (12, 23, "12-23"), (24, 35, "24-35"),
                        (36, 47, "36-47"), (48, 59, "48-59")]:
        add("Child age (months)", lab, d["child_age_m"].between(lo, hi))
    add("Child sex", "Male", d["child_female"] == 0)
    add("Child sex", "Female", d["child_female"] == 1)
    add("Residence", "Urban", d["urban"] == 1)
    add("Residence", "Rural", d["urban"] == 0)
    for q, lab in enumerate(["Poorest", "Poorer", "Middle", "Richer", "Richest"], start=1):
        add("Wealth quintile", lab, d["wealth"] == q)
    for e, lab in enumerate(["No schooling", "Primary", "Secondary", "Higher"]):
        add("Mother's education", lab, d["mother_edu"] == e)
    add("Mother's anaemia status", "Anaemic", d["mother_anaemic"] == 1)
    add("Mother's anaemia status", "Not anaemic", d["mother_anaemic"] == 0)
    add("Mother pregnant at interview", "Yes", d["mother_pregnant"] == 1)
    add("Mother pregnant at interview", "No", d["mother_pregnant"] == 0)
    return pd.DataFrame(rows)


def roc_analysis(d: pd.DataFrame, design: SurveyDesign, boot) -> dict:
    """ROC for continuous maternal haemoglobin (threshold-relative scale)."""
    score = d["mother_hb_rel"].to_numpy(float)
    truth = d["child_anaemic"].to_numpy(float)
    w = d["w"].to_numpy(float)

    fpr, tpr, thr, auc = weighted_roc(score, truth, w, higher_is_positive=False)

    # Bootstrap the AUC by re-weighting a single pre-sorted marker.
    order = np.argsort(score, kind="mergesort")     # low Hb first = positive
    y_sorted = truth[order]
    aucs = np.array([auc_from_sorted(y_sorted, wb[order]) for wb in boot])

    # Performance at a grid of absolute maternal haemoglobin cut-offs, so a
    # reader can see that no cut-off rescues the rule.
    rows = []
    for cut in np.arange(8.0, 13.5, 0.5):
        test = (d["mother_hb"].to_numpy(float) < cut).astype(float)
        a = accuracy(design, test, truth)
        rows.append({
            "maternal_hb_cutoff_g_dl": cut,
            "children_tested_pct": 100 * a["tested"],
            "sensitivity_pct": 100 * a["sensitivity"],
            "specificity_pct": 100 * a["specificity"],
            "ppv_pct": 100 * a["ppv"],
            "npv_pct": 100 * a["npv"],
            "meets_bar": bool(a["sensitivity"] >= BAR_SENS and a["npv"] >= BAR_NPV),
        })
    rows_df = pd.DataFrame(rows)
    ok = rows_df[rows_df["sensitivity_pct"] >= 90]
    best = float(ok["maternal_hb_cutoff_g_dl"].min()) if len(ok) else None

    return {
        "auc": float(auc),
        "auc_lcl": float(np.percentile(aucs, 2.5)),
        "auc_ucl": float(np.percentile(aucs, 97.5)),
        "fpr": fpr, "tpr": tpr,
        "threshold_table": rows,
        "best_sens90": best,
        "npv_at_sens90": float(ok["npv_pct"].max()) if len(ok) else None,
    }


def sensitivity_analyses(d: pd.DataFrame) -> list[dict]:
    out = []

    def run(label: str, frame: pd.DataFrame, test_col, truth_col):
        sd = SurveyDesign(frame["w"], frame["psu"], frame["stratum"])
        a = accuracy(sd, np.asarray(test_col, float), np.asarray(truth_col, float))
        out.append({
            "analysis": label,
            "n_children": int(len(frame)),
            "prevalence_pct": 100 * a["prevalence"],
            "children_tested_pct": 100 * a["tested"],
            "sensitivity_pct": 100 * a["sensitivity"],
            "sensitivity_lcl": 100 * a["sensitivity_lcl"],
            "sensitivity_ucl": 100 * a["sensitivity_ucl"],
            "specificity_pct": 100 * a["specificity"],
            "npv_pct": 100 * a["npv"],
            "npv_lcl": 100 * a["npv_lcl"],
            "npv_ucl": 100 * a["npv_ucl"],
            "meets_bar": bool(a["sensitivity"] >= BAR_SENS and a["npv"] >= BAR_NPV),
        })

    run("Primary analysis (rule R1, all mothers)", d,
        d["mother_anaemic"], d["child_anaemic"])

    np_ = d[d["mother_pregnant"] == 0]
    run("S1: non-pregnant mothers only", np_,
        np_["mother_anaemic"], np_["child_anaemic"])

    run("S2: reference standard moderate/severe child anaemia (<10.0 g/dL)", d,
        d["mother_anaemic"], d["child_mod_sev"])

    older = d[d["child_age_m"] >= 12]
    run("S3: children aged 12-59 months only", older,
        older["mother_anaemic"], older["child_anaemic"])

    unadj = d[d["child_hb_unadj"].notna() & d["mother_hb_unadj"].notna()]
    run("S4: unadjusted haemoglobin for mother and child", unadj,
        unadj["mother_anaemic_unadj"], (unadj["child_hb_unadj"] < 11.0).astype(float))

    return out


def weighted_corr(x, y, w) -> float:
    x, y, w = np.asarray(x, float), np.asarray(y, float), np.asarray(w, float)
    mx = np.average(x, weights=w)
    my = np.average(y, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    return cov / np.sqrt(np.average((x - mx) ** 2, weights=w)
                         * np.average((y - my) ** 2, weights=w))


if __name__ == "__main__":
    main()
