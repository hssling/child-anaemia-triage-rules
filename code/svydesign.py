"""Design-based estimation for stratified multistage survey data.

Implements the ultimate-cluster variance estimator with Taylor linearisation
for ratio estimators, matching the behaviour of Stata's `svy: ratio` and
`svy: proportion` (without finite population correction).

Written for Manuscript 12 of the Anaemia in India Research Programme. Kept
free of any study-specific logic so it can be reused and audited independently.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

Z_975 = 1.959963984540054


# --------------------------------------------------------------------------
# Design container
# --------------------------------------------------------------------------
class SurveyDesign:
    """Stratified, single-stage-of-clustering (ultimate cluster) design.

    Parameters
    ----------
    weights : array of individual sampling weights
    psu     : array of primary sampling unit identifiers
    strata  : array of stratum identifiers

    Singleton strata (one PSU) have no within-stratum degrees of freedom. They
    are collapsed into the nearest non-singleton stratum, and the number
    collapsed is retained in ``n_singletons_collapsed`` for reporting.
    """

    def __init__(self, weights, psu, strata):
        w = np.asarray(weights, dtype=float)
        psu = pd.Series(np.asarray(psu)).astype("category").cat.codes.to_numpy()
        strata = pd.Series(np.asarray(strata)).astype("category").cat.codes.to_numpy()

        if not (len(w) == len(psu) == len(strata)):
            raise ValueError("weights, psu and strata must be the same length")
        if np.any(w < 0):
            raise ValueError("negative sampling weights")

        strata, self.n_singletons_collapsed = _collapse_singletons(strata, psu)

        self.w = w
        self.psu = psu
        self.strata = strata

        # Unique (stratum, PSU) cell index, used to aggregate to PSU totals.
        cell = pd.MultiIndex.from_arrays([strata, psu])
        self._cell_codes, self._cells = pd.factorize(cell)
        cell_frame = pd.DataFrame({"h": strata, "cell": self._cell_codes}).drop_duplicates()
        self._cell_stratum = (
            cell_frame.sort_values("cell")["h"].to_numpy()
        )
        self._n_cells = self._cell_codes.max() + 1

        # PSUs per stratum, and degrees of freedom.
        h_counts = np.bincount(self._cell_stratum)
        self._psu_per_stratum = h_counts
        self.n_strata = len(h_counts)
        self.n_psu = self._n_cells
        self.df = self.n_psu - self.n_strata

    # ---------------------------------------------------------------- totals
    def total(self, x) -> float:
        return float(np.sum(self.w * np.asarray(x, dtype=float)))

    def _psu_totals(self, z) -> np.ndarray:
        """Aggregate a linearised variate to (stratum, PSU) totals."""
        return np.bincount(self._cell_codes, weights=z, minlength=self._n_cells)

    def var_total(self, z) -> float:
        """Ultimate-cluster variance of a weighted total of variate ``z``.

        ``z`` must already include the sampling weight.
        """
        t = self._psu_totals(np.asarray(z, dtype=float))
        h = self._cell_stratum
        n_h = self._psu_per_stratum

        sum_h = np.bincount(h, weights=t, minlength=self.n_strata)
        mean_h = np.divide(
            sum_h, n_h, out=np.zeros_like(sum_h, dtype=float), where=n_h > 0
        )
        dev2 = (t - mean_h[h]) ** 2
        ss_h = np.bincount(h, weights=dev2, minlength=self.n_strata)

        scale = np.divide(
            n_h, (n_h - 1), out=np.zeros(self.n_strata), where=n_h > 1
        )
        return float(np.sum(scale * ss_h))

    # ---------------------------------------------------------------- ratios
    def ratio(self, num, den):
        """Design-based ratio of weighted totals, with linearised SE.

        Returns ``(estimate, standard_error)``. Returns ``(nan, nan)`` when the
        denominator has no weighted mass.
        """
        x = np.asarray(num, dtype=float)
        y = np.asarray(den, dtype=float)
        X = self.total(x)
        Y = self.total(y)
        if Y <= 0:
            return float("nan"), float("nan")
        R = X / Y
        z = self.w * (x - R * y)          # linearised variate
        se = np.sqrt(max(self.var_total(z), 0.0)) / Y
        return float(R), float(se)

    def proportion_ci(self, num, den, alpha=0.05):
        """Ratio estimate with a logit-transformed confidence interval.

        The logit transform keeps the interval inside (0, 1), which matters for
        sensitivity and predictive values estimated close to their bounds.
        """
        p, se = self.ratio(num, den)
        if not np.isfinite(p) or not np.isfinite(se):
            return p, se, float("nan"), float("nan")
        if p <= 0:
            return p, se, 0.0, 0.0
        if p >= 1:
            return p, se, 1.0, 1.0

        z_crit = _t_crit(self.df, alpha)
        se_logit = se / (p * (1 - p))
        lo_l = np.log(p / (1 - p)) - z_crit * se_logit
        hi_l = np.log(p / (1 - p)) + z_crit * se_logit
        return p, se, float(_expit(lo_l)), float(_expit(hi_l))

    # ------------------------------------------------------------- bootstrap
    def bootstrap_weights(self, n_replicates: int, seed: int = 20260805):
        """Rescaled bootstrap replicate weights (Rao–Wu–Yue), PSUs within strata.

        Yields one weight vector per replicate. PSUs are resampled with
        replacement within each stratum; a stratum with n_h PSUs draws n_h - 1
        so the replicate variance is unbiased for the ultimate-cluster
        estimator.
        """
        rng = np.random.default_rng(seed)
        h = self._cell_stratum
        n_h = self._psu_per_stratum

        # Order cells by stratum once, so each stratum occupies a contiguous
        # block. The whole resample is then a single vectorised bincount
        # instead of a Python loop over ~1,000 strata per replicate.
        order = np.argsort(h, kind="mergesort")
        h_sorted = h[order]
        offsets = np.r_[0, np.cumsum(n_h)][:-1]          # block start per stratum
        eligible = n_h >= 2
        m_h = np.where(eligible, n_h - 1, 0)             # draws per stratum

        draw_stratum = np.repeat(np.arange(self.n_strata), m_h)
        draw_n = n_h[draw_stratum].astype(float)
        draw_offset = offsets[draw_stratum]
        total_draws = int(m_h.sum())

        # Rao-Wu-Yue rescaling constants, one per stratum.
        with np.errstate(divide="ignore", invalid="ignore"):
            lam = np.sqrt(np.where(eligible, m_h / (n_h - 1), 0.0))
            scale = np.where(eligible, lam * n_h / np.maximum(m_h, 1), 0.0)
        base = np.where(eligible, 1.0 - lam, 1.0)

        base_by_cell = base[h_sorted]
        scale_by_cell = scale[h_sorted]

        inv_order = np.empty_like(order)
        inv_order[order] = np.arange(len(order))

        for _ in range(n_replicates):
            picks = draw_offset + np.floor(rng.random(total_draws) * draw_n).astype(np.int64)
            counts = np.bincount(picks, minlength=self._n_cells).astype(float)
            adj_sorted = base_by_cell + scale_by_cell * counts
            adj = adj_sorted[inv_order]
            yield self.w * adj[self._cell_codes]


# --------------------------------------------------------------------------
# Weighted ROC / AUC
# --------------------------------------------------------------------------
def weighted_roc(score, positive, weights, higher_is_positive=False):
    """Design-weighted ROC curve.

    ``score`` is the continuous marker (maternal haemoglobin). By default a
    *lower* score indicates the target condition, which is the case here, so
    the marker is negated internally.

    Returns ``(fpr, tpr, thresholds, auc)`` with the curve traced from
    (0, 0) to (1, 1).
    """
    s = np.asarray(score, dtype=float)
    if not higher_is_positive:
        s = -s
    y = np.asarray(positive, dtype=float)
    w = np.asarray(weights, dtype=float)

    order = np.argsort(-s, kind="mergesort")
    s, y, w = s[order], y[order], w[order]

    tp = np.cumsum(w * y)
    fp = np.cumsum(w * (1 - y))
    P, N = tp[-1], fp[-1]

    # Keep only the last index of each tied score block.
    keep = np.r_[np.diff(s) != 0, True]
    tpr = np.r_[0.0, tp[keep] / P]
    fpr = np.r_[0.0, fp[keep] / N]
    thr = np.r_[np.inf, s[keep]]

    auc = float(np.trapezoid(tpr, fpr))
    if not higher_is_positive:
        thr = -thr
    return fpr, tpr, thr, auc


def auc_from_sorted(y_sorted, w, P_index=None):
    """Fast AUC for a pre-sorted marker, used inside the bootstrap loop.

    ``y_sorted`` must already be ordered by the marker (descending in the
    positive-indicating direction). Only the weights change between bootstrap
    replicates, so the sort is done once by the caller.
    """
    tp = np.cumsum(w * y_sorted)
    fp = np.cumsum(w * (1 - y_sorted))
    P, N = tp[-1], fp[-1]
    if P <= 0 or N <= 0:
        return float("nan")
    tpr = np.r_[0.0, tp / P]
    fpr = np.r_[0.0, fp / N]
    return float(np.trapezoid(tpr, fpr))


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _collapse_singletons(strata, psu):
    """Merge strata containing a single PSU into the neighbouring stratum."""
    frame = pd.DataFrame({"h": strata, "psu": psu}).drop_duplicates()
    counts = frame.groupby("h")["psu"].nunique()
    singles = counts[counts < 2].index.to_numpy()
    if len(singles) == 0:
        return strata, 0

    strata = strata.copy()
    valid = np.sort(counts[counts >= 2].index.to_numpy())
    for s in singles:
        # Nearest non-singleton stratum by code, which preserves the survey's
        # own state-then-residence ordering of stratum identifiers.
        target = valid[np.argmin(np.abs(valid - s))]
        strata[strata == s] = target
    return strata, int(len(singles))


def _expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def _t_crit(df, alpha=0.05):
    """Student-t critical value; falls back to the normal for large df."""
    if df is None or df <= 0 or df > 2000:
        return Z_975
    from scipy import stats

    return float(stats.t.ppf(1 - alpha / 2, df))
