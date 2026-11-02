# Household triage rules for child haemoglobin testing in India

Analysis code for:

> **Can a mother's anaemia decide which child gets tested? Diagnostic accuracy and net benefit of household triage rules for child haemoglobin testing in India**
> Siddalingaiah H S. Manuscript submitted for publication, 2026.

The study asks whether a mother's anaemia status, alone or with routine household
characteristics, can decide which children aged 6–59 months receive a haemoglobin test.
It evaluates six pre-specified triage rules as diagnostic tests in 182,209 mother–child
pairs from India's National Family Health Survey 2019–21 (NFHS-5). It uses design-based
diagnostic accuracy and decision-curve analysis, nationally and in 36 states and union
territories. No rule met the pre-specified bar (sensitivity ≥90% and negative predictive
value ≥80%), and universal testing had the greatest net benefit at every threshold
examined.

## What this repository reproduces

Every table and figure in the article, from the original DHS Program file, with no manual
steps. A clean re-run reproduces the files in `outputs/tables/` byte for byte
(see [Verifying a run](#verifying-a-run)).

| Path | Contents |
|---|---|
| `PROTOCOL.md` | Analysis plan v1.0, written before the confirmatory analysis |
| `governance/PROTOCOL_FREEZE.txt` | SHA-256 of the analysis plan, recorded before the analysis was run |
| `code/svydesign.py` | Design-based estimators: Taylor linearisation, rescaled bootstrap, weighted ROC |
| `code/build_analysis_dataset.py` | Builds the mother–child analysis file from the NFHS-5 Children's Recode |
| `code/run_triage_analysis.py` | Confirmatory analysis: accuracy, likelihood ratios, state-level, decision curves, sensitivity analyses |
| `code/required_discrimination.py` | Prevalence–NPV arithmetic (negative likelihood ratio required for NPV 80%) |
| `code/additional_analyses.py` | Analyses not pre-specified: WHO 2024 cut-offs, decision thresholds and testing yield, correction for attenuation |
| `code/make_figures.py` | Figures 1–3 |
| `outputs/tables/` | All aggregate result tables, as reported |
| `figures/` | Figures 1–3 (PNG, 300 dpi) |

| Result table | Article |
|---|---|
| `table_1` … `table_4` | Tables 1–4 |
| `table_5` | Supplementary Table S1 |
| `table_6` | Supplementary Table S2 and Figure 3a |
| `table_7`, `table_9` | Supplementary Table S3 |
| `table_8`, `table_8b` | Supplementary Table S4 |
| `table_10` | Table 3 (tests per anaemic child) and decision thresholds |
| `table_11` | Correction for attenuation (Discussion) |

## Data

The NFHS-5 microdata are **not** included and must not be committed here. They are
available free of charge from [The DHS Program](https://dhsprogram.com/data/) after
registration and approval of a data request. Download the India Children's Recode in Stata
format (`IAKR7EDT.zip`) and place it two levels above the repository folder:

```
<parent>/
├── NFHS 5 microdata/
│   └── IAKR7EDT.zip
└── <any folder>/
    └── child-anaemia-triage-rules/   ← this repository
```

## Running

Python 3.14 was used. Install the pinned dependencies and run the five scripts in order
from the repository root:

```bash
pip install -r requirements.txt
python code/build_analysis_dataset.py      # mother-child analysis file
python code/run_triage_analysis.py         # confirmatory analysis (about 2 minutes)
python code/required_discrimination.py     # prevalence-NPV arithmetic
python code/additional_analyses.py         # analyses not pre-specified (about 1 minute)
python code/make_figures.py                # Figures 1-3
```

The scripts write to `data/`, `governance/`, `outputs/` and `figures/`. The `data/` folder
holds individual-level records derived from the microdata and is excluded by `.gitignore`.

## Verifying a run

```bash
python code/verify_outputs.py
```

This compares every regenerated result table with `EXPECTED_SHA256.txt`, and checks
`PROTOCOL.md` against the hash recorded in `governance/PROTOCOL_FREEZE.txt`
(`4d012f463a94eab3542cbc241c83535fd7e4f0008a46d2e0c133dd56d4ec6703`). Table hashes are
taken after normalising line endings, because pandas writes CRLF on Windows and LF
elsewhere. All lines should report `OK`.

## Licence and citation

Code is released under the MIT License. Use of the NFHS-5 data is governed by the DHS
Program terms of use. If you use this code, please cite the article and this repository
(see `CITATION.cff`).
