"""Why no triage rule can work: the prevalence-NPV arithmetic.

Pre-specified in PROTOCOL.md section 11 ("the paper reports the arithmetic
reason for the result"). This is a deterministic identity, not an estimate, so
it carries no confidence interval.

For a binary triage rule, the odds of the target condition among children the
rule declines to test are

    post-test odds = pre-test odds x LR-

so an NPV of at least `bar` requires

    LR- <= [(1 - bar) / bar] / pre-test odds

Usage:
    python code/required_discrimination.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "outputs" / "tables"
BAR_NPV = 0.80


def required_lr_neg(prevalence: np.ndarray, bar: float = BAR_NPV) -> np.ndarray:
    pre_odds = prevalence / (1 - prevalence)
    target_odds = (1 - bar) / bar
    return target_odds / pre_odds


def max_npv(prevalence: np.ndarray, lr_neg: float) -> np.ndarray:
    post_odds = (prevalence / (1 - prevalence)) * lr_neg
    return 1 / (1 + post_odds)


def main() -> None:
    obs = pd.read_csv(TABLES / "table_2_rule_accuracy.csv")
    best_lr_neg = float(obs["lr_neg"].min())
    best_rule = obs.loc[obs["lr_neg"].idxmin(), "rule"]

    states = pd.read_csv(TABLES / "table_4_state_level_rule_R1.csv")

    prev_grid = np.arange(0.20, 0.96, 0.05)
    grid = pd.DataFrame({
        "child_anaemia_prevalence": prev_grid,
        "lr_neg_required_for_npv80": required_lr_neg(prev_grid),
        "npv_attainable_at_best_observed_lr_neg": max_npv(prev_grid, best_lr_neg),
        "sensitivity_required_at_observed_specificity_49pct":
            1 - required_lr_neg(prev_grid) * 0.4928,
    })
    grid["sensitivity_required_at_observed_specificity_49pct"] = grid[
        "sensitivity_required_at_observed_specificity_49pct"
    ].clip(lower=np.nan, upper=1.0)
    grid.to_csv(TABLES / "table_8_required_discrimination.csv", index=False)

    national_prev = float(
        pd.read_csv(TABLES / "table_3_per_100000_children.csv")
        .query("rule == 'TEST-ALL'")["anaemic_detected"].iloc[0]
    ) / 100_000

    summary = {
        "national_child_anaemia_prevalence": round(national_prev, 4),
        "pre_test_odds": round(national_prev / (1 - national_prev), 3),
        "lr_neg_required_for_npv80_national": round(
            float(required_lr_neg(np.array([national_prev]))[0]), 4),
        "best_observed_lr_neg": round(best_lr_neg, 3),
        "best_observed_lr_neg_rule": str(best_rule),
        "shortfall_factor": round(
            best_lr_neg / float(required_lr_neg(np.array([national_prev]))[0]), 1),
        "lowest_state_prevalence": round(float(states["child_anaemia_pct"].min()) / 100, 4),
        "lowest_state": str(states.loc[states["child_anaemia_pct"].idxmin(), "state"]),
        "lr_neg_required_in_lowest_state": round(float(required_lr_neg(
            np.array([states["child_anaemia_pct"].min() / 100]))[0]), 4),
        "max_prevalence_at_which_best_rule_reaches_npv80": round(float(
            max(prev for prev in np.arange(0.01, 0.99, 0.001)
                if max_npv(np.array([prev]), best_lr_neg)[0] >= BAR_NPV)), 3),
    }
    pd.DataFrame([summary]).to_csv(
        TABLES / "table_8b_arithmetic_summary.csv", index=False
    )
    for k, v in summary.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
