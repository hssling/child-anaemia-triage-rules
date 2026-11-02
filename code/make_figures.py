"""Print-ready figures for Manuscript 12.

Three figures, 300 dpi, written as both PNG and TIFF. Colour is validated
(dataviz six-checks, all-pairs mode) but every series also carries a distinct
line style and marker, so the figures survive greyscale journal reproduction.

Usage:
    python code/make_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "outputs" / "tables"
FIGS = ROOT / "figures"

# Validated categorical slots (dataviz all-pairs PASS) + neutral reference ink.
BLUE, ORANGE, VIOLET = "#2a78d6", "#eb6834", "#4a3aa7"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d8d7d2"

MM = 1 / 25.4
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "axes.titleweight": "bold",
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def save(fig, name: str) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "tiff"):
        fig.savefig(FIGS / f"{name}.{ext}", format=ext, dpi=300,
                    facecolor="white", pil_kwargs={"compression": "tiff_lzw"}
                    if ext == "tiff" else None)
    plt.close(fig)
    print(f"wrote {name}.png / .tiff")


def panel(ax, letter: str) -> None:
    ax.text(-0.16, 1.04, f"({letter})", transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", color=INK, ha="left", va="bottom")


def recede(ax, y_grid=True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if y_grid:
        ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)


# --------------------------------------------------------------------- Fig 1
def figure1_flow() -> None:
    """STARD participant flow."""
    import json
    log = json.loads((ROOT / "governance" / "DATASET_BUILD_MANIFEST.json")
                     .read_text(encoding="utf-8"))
    steps = log["exclusion_steps"]

    fig, ax = plt.subplots(figsize=(120 * MM, 150 * MM))
    ax.set_xlim(0, 10.6)
    ax.set_ylim(0, 10)
    ax.axis("off")

    boxes = [
        (f"Children in NFHS-5 Children's Recode\nn = {log['rows_in_children_recode']:,}", None),
        (f"Living children\nn = {steps[0]['remaining']:,}",
         f"Deceased children excluded\nn = {steps[0]['excluded']:,}"),
        (f"Aged 6–59 completed months\nn = {steps[1]['remaining']:,}",
         f"Outside age window\nn = {steps[1]['excluded']:,}"),
        (f"Valid child haemoglobin\n(reference standard)\nn = {steps[2]['remaining']:,}",
         f"Not measured, absent, refused\nor invalid result\nn = {steps[2]['excluded']:,}"),
        (f"Valid maternal haemoglobin\n(index test)\nn = {steps[3]['remaining']:,}",
         f"Mother not measured or\ninvalid result\nn = {steps[3]['excluded']:,}"),
        (f"ANALYTIC SAMPLE\n{log['analytic_children']:,} children of "
         f"{log['analytic_mothers']:,} mothers\nin {log['analytic_psus']:,} sampling clusters", None),
    ]

    y = 9.3
    dy = 1.62
    for i, (main, excl) in enumerate(boxes):
        final = i == len(boxes) - 1
        w, x0 = 5.5, 0.2
        h = 0.86 if not final else 1.0
        ax.add_patch(Rectangle((x0, y - h), w, h, facecolor="#f4f4f2" if not final else "#e8eef8",
                               edgecolor=INK if final else MUTED,
                               linewidth=1.0 if final else 0.7, zorder=2))
        ax.text(x0 + w / 2, y - h / 2, main, ha="center", va="center",
                fontsize=7.2, color=INK, zorder=3,
                fontweight="bold" if final else "normal")

        if excl is not None:
            # Branch off the vertical connector midway between the two boxes,
            # so the arrow never starts inside the box above.
            y_branch = y + (dy - h) / 2
            ax.add_patch(Rectangle((6.35, y_branch - 0.41), 4.15, 0.82, facecolor="white",
                                   edgecolor=MUTED, linewidth=0.6,
                                   linestyle=(0, (2, 1.5)), zorder=2))
            ax.text(6.35 + 4.15 / 2, y_branch, excl, ha="center", va="center",
                    fontsize=6.6, color=MUTED, zorder=3)
            ax.add_patch(FancyArrowPatch((x0 + w / 2, y_branch), (6.35, y_branch),
                                         arrowstyle="-|>", mutation_scale=7,
                                         color=MUTED, linewidth=0.7, zorder=3))

        if i < len(boxes) - 1:
            ax.add_patch(FancyArrowPatch((x0 + w / 2, y - h), (x0 + w / 2, y - dy),
                                         arrowstyle="-|>", mutation_scale=8,
                                         color=INK, linewidth=0.8, zorder=3))
        y -= dy

    save(fig, "figure_1_participant_flow")


# --------------------------------------------------------------------- Fig 2
def figure2_roc_and_frontier() -> None:
    roc = np.load(ROOT / "data" / "roc_curve.npy")
    fpr, tpr = roc[0], roc[1]
    import json
    man = json.loads((ROOT / "outputs" / "ANALYSIS_MANIFEST.json").read_text(encoding="utf-8"))
    t3 = pd.read_csv(TABLES / "table_3_per_100000_children.csv")

    fig, axes = plt.subplots(1, 2, figsize=(170 * MM, 78 * MM))

    # ---- (a) ROC
    ax = axes[0]
    ax.plot([0, 1], [0, 1], color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)), zorder=1)
    ax.text(0.55, 0.47, "no discrimination", fontsize=6.8, color=MUTED,
            rotation=39, ha="center", va="center")
    ax.plot(fpr, tpr, color=BLUE, linewidth=2.0, zorder=3)
    ax.text(0.40, 0.20,
            f"AUC {man['auc']:.3f}\n(95% CI {man['auc_lcl']:.3f}–{man['auc_ucl']:.3f})",
            fontsize=7.6, color=INK, ha="left", va="center")
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    recede(ax, y_grid=False)
    ax.grid(color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    panel(ax, "a")

    # ---- (b) tests-avoided vs anaemic children missed
    ax = axes[1]
    rules = t3[~t3["rule"].isin(["TEST-ALL", "TEST-NONE"])]
    x = rules["tests_avoided"] / 1000
    y = rules["anaemic_missed"] / 1000

    all_row = t3[t3["rule"] == "TEST-ALL"].iloc[0]
    none_row = t3[t3["rule"] == "TEST-NONE"].iloc[0]

    # The line every rule would have to beat: proportional loss.
    ax.plot([0, none_row["tests_avoided"] / 1000],
            [0, none_row["anaemic_missed"] / 1000],
            color=MUTED, linewidth=0.9, linestyle=(0, (3, 2)), zorder=1)
    ax.text(84, 9, "children missed in exact\nproportion to tests avoided",
            fontsize=6.8, color=MUTED, ha="center", va="center")

    ax.scatter(x, y, s=46, color=ORANGE, edgecolor="white", linewidth=1.2, zorder=4)
    # Offsets are nudged per rule: R4 and R5 sit almost on top of each other.
    offsets = {"R1": (7, 2), "R2": (2, -13), "R3": (7, 2),
               "R4": (8, 0), "R5": (-19, 0), "R6": (8, -3)}
    for _, r in rules.iterrows():
        ax.annotate(r["rule"], (r["tests_avoided"] / 1000, r["anaemic_missed"] / 1000),
                    textcoords="offset points", xytext=offsets[r["rule"]],
                    fontsize=7.2, color=INK, fontweight="bold")

    ax.scatter([0], [0], s=52, marker="s", color=VIOLET,
               edgecolor="white", linewidth=1.2, zorder=4)
    ax.annotate("Test all", (0, 0), textcoords="offset points", xytext=(7, 1),
                fontsize=7.2, color=VIOLET, fontweight="bold")
    ax.scatter([none_row["tests_avoided"] / 1000], [none_row["anaemic_missed"] / 1000],
               s=52, marker="D", color=INK, edgecolor="white", linewidth=1.2, zorder=4)
    ax.annotate("Test none", (none_row["tests_avoided"] / 1000,
                              none_row["anaemic_missed"] / 1000),
                textcoords="offset points", xytext=(-8, 6), ha="right",
                fontsize=7.2, color=INK, fontweight="bold")

    ax.set_xlabel("Haemoglobin tests avoided\n(thousands per 100 000 children)")
    ax.set_ylabel("Anaemic children missed\n(thousands per 100 000)")
    ax.set_xlim(-6, 112); ax.set_ylim(-5, 78)
    recede(ax)
    panel(ax, "b")

    fig.subplots_adjust(wspace=0.30)
    save(fig, "figure_2_discrimination_and_frontier")


# --------------------------------------------------------------------- Fig 3
def figure3_decision_and_states() -> None:
    dca = pd.read_csv(TABLES / "table_6_decision_curve.csv")
    states = pd.read_csv(TABLES / "table_4_state_level_rule_R1.csv")
    arith = pd.read_csv(TABLES / "table_8_required_discrimination.csv")
    summ = pd.read_csv(TABLES / "table_8b_arithmetic_summary.csv").iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(170 * MM, 78 * MM))

    # ---- (a) decision curve
    ax = axes[0]
    ax.plot(dca["threshold"], dca["test_all"], color=VIOLET, linewidth=2.0,
            label="Test all children", zorder=4)
    ax.plot(dca["threshold"], dca["test_none"], color=INK, linewidth=1.2,
            linestyle=(0, (1, 1.6)), label="Test no children", zorder=3)
    ax.plot(dca["threshold"], dca["R1"], color=BLUE, linewidth=1.6,
            linestyle=(0, (4, 1.6)), label="R1  mother anaemic", zorder=3)
    ax.plot(dca["threshold"], dca["R6"], color=ORANGE, linewidth=1.6,
            linestyle=(0, (5, 1.4, 1, 1.4)), label="R6  broadest rule", zorder=3)

    ax.set_xlabel("Threshold probability of child anaemia\nat which testing is worthwhile")
    ax.set_ylabel("Net benefit")
    ax.set_xlim(0.01, 0.50)
    ax.legend(loc="lower left", bbox_to_anchor=(0.015, 0.02),
              handlelength=2.4, frameon=True, facecolor="white",
              framealpha=0.96, edgecolor="none", borderpad=0.45,
              labelspacing=0.3)
    recede(ax)
    panel(ax, "a")

    # ---- (b) state NPV against prevalence, with the attainable ceiling
    ax = axes[1]
    prev = arith["child_anaemia_prevalence"] * 100
    ceiling = arith["npv_attainable_at_best_observed_lr_neg"] * 100
    ax.plot(prev, ceiling, color=MUTED, linewidth=1.2, linestyle=(0, (3, 2)), zorder=2)
    ax.text(85, 71, "Best attainable NPV\nat strongest observed rule\n"
            f"(LR− {summ['best_observed_lr_neg']:.2f})",
            fontsize=6.8, color=MUTED, ha="center", va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.92,
                      pad=2.0), zorder=6)

    ax.axhline(80, color=INK, linewidth=1.0, zorder=3)
    ax.text(97, 81.5, "required NPV ≥ 80%", fontsize=7.2, color=INK,
            ha="right", va="bottom", fontweight="bold")

    ax.scatter(states["child_anaemia_pct"], states["npv_pct"], s=26,
               color=BLUE, edgecolor="white", linewidth=0.9, zorder=5)
    label_offsets = {
        "kerala": (8, 7),
        "ladakh": (8, 7),
        "gujarat": (8, -3),
    }
    for label in ["kerala", "ladakh", "gujarat"]:
        row = states[states["state"] == label]
        if len(row):
            r = row.iloc[0]
            ax.annotate(r["state"].title(), (r["child_anaemia_pct"], r["npv_pct"]),
                        textcoords="offset points", xytext=label_offsets[label],
                        fontsize=6.9, color=INK, zorder=7,
                        bbox=dict(facecolor="white", edgecolor="none",
                                  alpha=0.86, pad=1.1))

    ax.set_xlabel("Child anaemia prevalence in the state (%)")
    ax.set_ylabel("Negative predictive value of\n“mother anaemic” (%)")
    ax.set_xlim(18, 100); ax.set_ylim(0, 100)
    recede(ax)
    panel(ax, "b")

    fig.subplots_adjust(wspace=0.32)
    save(fig, "figure_3_net_benefit_and_states")


if __name__ == "__main__":
    figure1_flow()
    figure2_roc_and_frontier()
    figure3_decision_and_states()
