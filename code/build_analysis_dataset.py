"""Build the mother-child analysis dataset for Manuscript 12.

Runs from the original DHS ZIP archive with no manual steps. Writes a derived
parquet file, a variable dictionary and a construction manifest recording every
exclusion.

Usage:
    python code/build_analysis_dataset.py
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                      # manuscripts/12_child_anaemia_triage_rules
REPO = ROOT.parent.parent               # d:/Anemia research
SOURCE_ZIP = REPO / "NFHS 5 microdata" / "IAKR7EDT.zip"
MEMBER = "IAKR7EFL.DTA"
OUT_DATA = ROOT / "data"

KEEP = [
    "caseid", "midx",
    "v001", "v002", "v005", "v021", "v022", "v023", "v024", "v025",
    "v012", "v106", "v190", "v213",
    "v453", "v456", "v457",
    "hw1", "hw53", "hw55", "hw56", "hw57",
    "b4", "b5", "b19", "bord",
]

# DHS special codes for haemoglobin: 994-999 are non-numeric result codes.
HB_VALID_MAX = 900          # deci-g/dL, i.e. 90.0 g/dL — well above any real value


def main() -> None:
    OUT_DATA.mkdir(parents=True, exist_ok=True)
    log: dict = {
        "manuscript": "12_child_anaemia_triage_rules",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "source_zip": str(SOURCE_ZIP),
        "source_zip_sha256": sha256_file(SOURCE_ZIP),
        "source_member": MEMBER,
    }

    tmp = Path(tempfile.mkdtemp(prefix="m12_"))
    try:
        with zipfile.ZipFile(SOURCE_ZIP) as z:
            z.extract(MEMBER, tmp)
        dta = tmp / MEMBER
        log["source_member_sha256"] = sha256_file(dta)

        df, meta = pyreadstat.read_dta(str(dta), usecols=KEEP)
        state_labels = meta.variable_value_labels.get("v024", {})
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    log["rows_in_children_recode"] = int(len(df))

    # ---------------------------------------------------------------- filters
    steps = []

    def drop(mask: pd.Series, label: str, frame: pd.DataFrame) -> pd.DataFrame:
        kept = frame[mask]
        steps.append({"step": label, "excluded": int(len(frame) - len(kept)),
                      "remaining": int(len(kept))})
        return kept

    d = df.copy()
    d = drop(d["b5"] == 1, "child alive", d)
    d = drop(d["b19"].between(6, 59), "child aged 6-59 completed months", d)
    d = drop(d["hw56"].notna() & (d["hw56"] < HB_VALID_MAX),
             "valid altitude-adjusted child haemoglobin", d)
    d = drop(d["v456"].notna() & (d["v456"] < HB_VALID_MAX),
             "valid adjusted maternal haemoglobin", d)
    d = drop(d["v005"].notna() & (d["v005"] > 0), "positive sampling weight", d)

    log["exclusion_steps"] = steps
    log["analytic_children"] = int(len(d))
    log["analytic_mothers"] = int(d["caseid"].nunique())
    log["analytic_psus"] = int(d["v021"].nunique())

    # ------------------------------------------------------------- variables
    out = pd.DataFrame(index=d.index)
    out["caseid"] = d["caseid"].astype(str)
    out["childid"] = d["caseid"].astype(str) + "_" + d["midx"].astype(int).astype(str)

    out["w"] = d["v005"] / 1e6
    out["psu"] = d["v021"].astype(int)
    out["stratum"] = d["v022"].astype(int)

    out["state_code"] = d["v024"].astype(int)
    out["state"] = out["state_code"].map(
        {int(k): str(v) for k, v in state_labels.items()}
    ).fillna("Unknown")
    out["urban"] = (d["v025"] == 1).astype(int)
    out["wealth"] = d["v190"].astype("Int64")
    out["poorest2"] = d["v190"].isin([1, 2]).astype(int)
    out["mother_edu"] = d["v106"].astype("Int64")
    out["mother_age"] = d["v012"].astype("Int64")

    out["child_age_m"] = d["b19"].astype(int)
    out["child_female"] = (d["b4"] == 2).astype(int)
    out["birth_order"] = d["bord"].astype("Int64")

    # Haemoglobin, converted from deci-g/dL to g/dL.
    out["child_hb"] = d["hw56"] / 10.0
    out["mother_hb"] = d["v456"] / 10.0
    out["child_hb_unadj"] = np.where(
        d["hw53"].notna() & (d["hw53"] < HB_VALID_MAX), d["hw53"] / 10.0, np.nan
    )
    out["mother_hb_unadj"] = np.where(
        d["v453"].notna() & (d["v453"] < HB_VALID_MAX), d["v453"] / 10.0, np.nan
    )

    out["mother_pregnant"] = (d["v213"] == 1).astype(int)

    # ------------------------------------------------- target condition
    out["child_anaemic"] = (out["child_hb"] < 11.0).astype(int)
    out["child_mod_sev"] = (out["child_hb"] < 10.0).astype(int)
    out["child_severity"] = pd.cut(
        out["child_hb"],
        bins=[-np.inf, 7.0, 10.0, 11.0, np.inf],
        labels=["severe", "moderate", "mild", "none"],
        right=False,
    )

    # ------------------------------------------------- maternal anaemia
    # Pregnancy-appropriate thresholds, per WHO and the AMB programme.
    preg = out["mother_pregnant"] == 1
    out["mother_anaemic"] = np.where(
        preg, out["mother_hb"] < 11.0, out["mother_hb"] < 12.0
    ).astype(int)
    out["mother_mod_sev"] = np.where(
        preg, out["mother_hb"] < 10.0, out["mother_hb"] < 9.0
    ).astype(int)
    out["mother_anaemic_unadj"] = np.where(
        preg, out["mother_hb_unadj"] < 11.0, out["mother_hb_unadj"] < 12.0
    ).astype(float)

    # Threshold-relative maternal haemoglobin: distance below the applicable
    # cut-off, so pregnant and non-pregnant mothers sit on one continuous scale
    # for the ROC analysis.
    out["mother_hb_rel"] = out["mother_hb"] - np.where(preg, 11.0, 12.0)

    out["child_6_23"] = (out["child_age_m"] < 24).astype(int)

    out = out.reset_index(drop=True)

    path = OUT_DATA / "analysis_mother_child_pairs_v1.0.parquet"
    out.to_parquet(path, index=False)
    log["output"] = str(path.relative_to(ROOT))
    log["output_sha256"] = sha256_file(path)
    log["output_rows"] = int(len(out))

    # Weighted prevalence check against the published NFHS-5 national figure.
    W = out["w"].sum()
    log["weighted_child_anaemia_pct"] = round(
        100 * float((out["w"] * out["child_anaemic"]).sum() / W), 2
    )
    log["weighted_maternal_anaemia_pct"] = round(
        100 * float((out["w"] * out["mother_anaemic"]).sum() / W), 2
    )
    log["published_nfhs5_child_anaemia_pct"] = 67.1

    (ROOT / "governance" / "DATASET_BUILD_MANIFEST.json").write_text(
        json.dumps(log, indent=2), encoding="utf-8"
    )
    print(json.dumps(log, indent=2))


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


if __name__ == "__main__":
    main()
