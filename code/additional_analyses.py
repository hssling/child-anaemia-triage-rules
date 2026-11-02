"""Additional analyses for Manuscript 12, run after the confirmatory analysis.

None of these were pre-specified; each is reported as such in the manuscript.

1. WHO 2024 reference standard. The 2024 WHO guideline, adopted by Anemia
   Mukt Bharat Abhiyaan in 2026, lowers the anaemia cut-off to <10.5 g/dL at
   6-23 months and keeps <11.0 g/dL at 24-59 months. All six rules are
   re-evaluated against that reference standard with the same design-based
   estimators and acceptability bar.
2. Decision thresholds and testing yield. For each rule, the threshold
   probability above which it would have greater net benefit than universal
   testing (1 - NPV) and below which it beats testing no one (PPV), and the
   number of tests performed per anaemic child identified. These are exact
   algebraic consequences of the pre-specified estimates.
3. Correction of the maternal-child haemoglobin correlation for measurement
   error (Spearman), across assumed reliabilities.

Usage:
    python code/additional_analyses.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svydesign import SurveyDesign  # noqa: E402
from run_triage_analysis import (BAR_NPV, BAR_NPV_LCL, BAR_SENS,  # noqa: E402
                                 BAR_SENS_LCL, BOOT_SEED, N_BOOT, RULE_LABELS,
                                 accuracy, build_rules, likelihood_ratios,
                                 weighted_corr)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data" / "analysis_mother_child_pairs_v1.0.parquet"
TABLES = ROOT / "outputs" / "tables"


def who2024_reference(d: pd.DataFrame) -> pd.Series:
    cut = np.where(d["child_age_m"] < 24, 10.5, 11.0)
    return (d["child_hb"] < cut).astype(float)


def main() -> None:
    d = pd.read_parquet(DATA)
    design = SurveyDesign(d["w"], d["psu"], d["stratum"])
    boot = list(design.bootstrap_weights(N_BOOT, seed=BOOT_SEED))
    rules = build_rules(d)

    # ---------------------------------------------- 1. WHO 2024 cut-offs
    truth = who2024_reference(d).to_numpy()
    rows = []
    for rule in RULE_LABELS:
        test = rules[rule].to_numpy(float)
        a = accuracy(design, test, truth)
        lr = likelihood_ratios(boot, test, truth)
        sens, spec = a["sensitivity"], a["specificity"]
        meets = (sens >= BAR_SENS and a["npv"] >= BAR_NPV
                 and a["sensitivity_lcl"] > BAR_SENS_LCL
                 and a["npv_lcl"] > BAR_NPV_LCL)
        rows.append({
            "rule": rule,
            "definition": RULE_LABELS[rule],
            "prevalence_pct": 100 * a["prevalence"],
            "children_tested_pct": 100 * a["tested"],
            "sensitivity_pct": 100 * sens,
            "sensitivity_lcl": 100 * a["sensitivity_lcl"],
            "sensitivity_ucl": 100 * a["sensitivity_ucl"],
            "npv_pct": 100 * a["npv"],
            "npv_lcl": 100 * a["npv_lcl"],
            "npv_ucl": 100 * a["npv_ucl"],
            "lr_neg": (1 - sens) / spec,
            "lr_neg_lcl": lr["lr_neg_lcl"],
            "lr_neg_ucl": lr["lr_neg_ucl"],
            "meets_bar": bool(meets),
        })
    t9 = pd.DataFrame(rows)
    t9.to_csv(TABLES / "table_9_who2024_reference_standard.csv", index=False)

    prev = float(t9["prevalence_pct"].iloc[0]) / 100
    odds = prev / (1 - prev)
    lr_required = (1 - BAR_NPV) / BAR_NPV / odds

    # ---------------------------------------------- 2. thresholds and yield
    t2 = pd.read_csv(TABLES / "table_2_rule_accuracy.csv")
    t3 = pd.read_csv(TABLES / "table_3_per_100000_children.csv")
    all_row = t3[t3["rule"] == "TEST-ALL"].iloc[0]
    yrows = []
    for _, r in t3.iterrows():
        tested, found = r["tests_performed"], r["anaemic_detected"]
        row = {"strategy": r["rule"],
               "tests_per_anaemic_child_detected":
                   tested / found if found else np.nan}
        if r["rule"] not in ("TEST-ALL", "TEST-NONE"):
            acc = t2[t2["rule"] == r["rule"]].iloc[0]
            extra_tests = all_row["tests_performed"] - tested
            extra_found = all_row["anaemic_detected"] - found
            row.update({
                "additional_tests_per_additional_child_if_all_tested":
                    extra_tests / extra_found,
                "threshold_rule_beats_test_all_pct": 100 - acc["npv_pct"],
                "threshold_rule_beats_test_all_lcl": 100 - acc["npv_ucl"],
                "threshold_rule_beats_test_all_ucl": 100 - acc["npv_lcl"],
                "threshold_rule_beats_test_none_pct": acc["ppv_pct"],
                "harm_benefit_ratio_needed_to_prefer_rule":
                    (100 - acc["npv_pct"]) / acc["npv_pct"],
            })
        yrows.append(row)
    t10 = pd.DataFrame(yrows)
    t10.to_csv(TABLES / "table_10_thresholds_and_testing_yield.csv", index=False)

    # Consistency check against the pre-specified decision-curve grid.
    dca = pd.read_csv(TABLES / "table_6_decision_curve.csv")
    for rule in RULE_LABELS:
        thr = float(t10.loc[t10["strategy"] == rule,
                            "threshold_rule_beats_test_all_pct"].iloc[0]) / 100
        assert (dca[rule] <= dca["test_all"] + 1e-12).all(), rule
        assert thr > dca["threshold"].max(), rule

    # ---------------------------------------------- 3. disattenuation
    r_obs = float(weighted_corr(d["child_hb"], d["mother_hb"], d["w"]))
    dis = [{"assumed_reliability_each": rel,
            "disattenuated_correlation": r_obs / rel}
           for rel in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)]
    pd.DataFrame(dis).to_csv(TABLES / "table_11_disattenuated_correlation.csv",
                             index=False)

    summary = {
        "who2024_prevalence_pct": 100 * prev,
        "who2024_lr_neg_required_for_npv80": lr_required,
        "who2024_rules_meeting_bar": int(t9["meets_bar"].sum()),
        "observed_correlation": r_obs,
    }
    (ROOT / "outputs" / "ADDITIONAL_ANALYSES.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    pd.set_option("display.width", 200)
    print(t9.round(3).to_string())
    print(t10.round(3).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
