"""Read-depth bias correction and the depth statistics around it.

Last update: 2026-08-08

Functions:
- correct_readcount_lowess: HMMcopy-style LOWESS correction on GC, MAP, REPLI
- correct_readcount_quadreg: median quadratic-regression correction, the default
- correct_readcount_by_target_sites: run either corrector on- and off-target separately
- normalize_library_by_target: one library-size factor per (dataset, capture group)
- compute_gc_rd_stats: GC-vs-depth correlation and spread, before and after
- compute_depth_statistics: per-dataset depth summary written to depth_statistics.tsv
"""

import logging

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.interpolate import interp1d
from scipy.stats import pearsonr, spearmanr
from statsmodels.nonparametric.smoothers_lowess import lowess

from utils import sort_chroms


def correct_readcount_lowess(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
    lowess_frac_tight=0.03,
    lowess_frac_smooth=0.3,
    grid_size=1001,
    seed=42,
):
    """HMMcopy-style correctReadcount: GC + optional mappability + optional
    replication timing stages (ACEseq-style sequential LOWESS correction).

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per fixed bin (1-D).
    gc : np.ndarray
        Per-bin GC fraction in [0, 1].
    mappability : np.ndarray or None
        Per-bin mappability values in [0, 1]. If None, stage 2 is skipped.
    repliseq : np.ndarray or None
        Per-bin consensus replication timing score. If None, stage 3 is
        skipped.  Higher values = earlier replication = higher expected
        coverage in cycling cells.
    samplesize : int
        Max number of ideal bins for loess fitting.
    routlier : float
        Upper quantile for read-count outlier removal.
    doutlier : float
        Quantile for GC/mappability domain outlier removal (top/bottom).
    min_mappability : float
        Minimum mappability for ideal bins.
    lowess_frac_tight : float
        Bandwidth fraction for the first (tight) LOWESS pass.
    lowess_frac_smooth : float
        Bandwidth fraction for the second (smoothing) LOWESS pass.
    grid_size : int
        Number of grid points for LOWESS interpolation.
    seed : int
        Base random seed for subsampling ideal bins (incremented per stage).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*). Zero depth stays 0.0; NaN where
        correction is not possible.
    float
        RMSE from the GC LOWESS fit (root mean squared error between raw
        depth and LOWESS prediction at valid bins).
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    def _dedup_sorted(xy):
        """Remove duplicate x-values from sorted LOWESS output to avoid
        division-by-zero in interp1d."""
        _, idx = np.unique(xy[:, 0], return_index=True)
        return xy[idx]

    def _dedup_input(y, x):
        """Average y-values for duplicate x-values to avoid degenerate
        local regressions (division by zero) inside LOWESS."""
        ux, inv = np.unique(x, return_inverse=True)
        if len(ux) == len(x):
            return y, x
        uy = np.zeros(len(ux), dtype=np.float64)
        np.add.at(uy, inv, y)
        counts = np.bincount(inv).astype(np.float64)
        uy /= counts
        return uy, ux

    def _fit_lowess_interp(y, x, grid):
        """Tight LOWESS -> grid smooth -> final interpolator, no extrapolation.

        Neither stage predicts outside the fitted covariate range: R's ``predict.loess``
        returns NA there, and the grid spans the whole covariate domain, so extrapolating
        would smooth a fabricated tail back into the data edge.
        """
        y, x = _dedup_input(y, x)
        if len(x) < 2:
            return None
        s1 = lowess(y, x, frac=lowess_frac_tight, return_sorted=True)
        s1 = _dedup_sorted(s1)
        if len(s1) < 2:
            return None
        s1_fn = interp1d(
            s1[:, 0],
            s1[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        on_grid = s1_fn(grid)
        in_support = np.isfinite(on_grid)
        if in_support.sum() < 2:
            return None
        s2 = lowess(
            on_grid[in_support],
            grid[in_support],
            frac=lowess_frac_smooth,
            return_sorted=True,
        )
        s2 = _dedup_sorted(s2)
        if len(s2) < 2:
            return None
        return interp1d(
            s2[:, 0],
            s2[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

    def _apply_stage(
        prev,
        covariate,
        cov_lo,
        cov_hi,
        grid,
        seed,
        stage_name="",
        extra_ideal_mask=None,
    ):
        """Apply one LOWESS correction stage and return (corrected, rmse)."""
        valid = (prev > 0) & np.isfinite(covariate)
        val_hi = np.nanquantile(prev[valid], 1.0 - routlier)
        ideal = valid & (prev <= val_hi) & (covariate >= cov_lo) & (covariate <= cov_hi)
        if extra_ideal_mask is not None:
            ideal &= extra_ideal_mask

        ideal_idx = np.where(ideal)[0]
        if ideal_idx.size < 10:
            logging.warning(
                f"    {stage_name}: only {ideal_idx.size} ideal bins; skipping correction"
            )
            return prev.copy(), 0.0
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(seed)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        interp_fn = _fit_lowess_interp(prev[ideal_idx], covariate[ideal_idx], grid)
        if interp_fn is None:
            logging.warning(
                f"    {stage_name}: LOWESS fit failed (too few distinct values); skipping correction"
            )
            return prev.copy(), 0.0
        predicted = interp_fn(covariate)

        # RMSE between raw values and LOWESS prediction (valid bins only)
        valid_pred = (predicted > 0) & (prev > 0)
        n_valid_pred = int(valid_pred.sum())
        rmse = (
            float(np.sqrt(np.mean((prev[valid_pred] - predicted[valid_pred]) ** 2)))
            if n_valid_pred > 0
            else 0.0
        )

        with np.errstate(invalid="ignore", divide="ignore"):
            corrected = np.where(predicted > 0, prev / predicted, np.nan)
        positive = (predicted > 0) & (prev > 0)
        if positive.any():
            scale = np.median(prev[positive]) / np.median(corrected[positive])
            corrected[np.isfinite(corrected)] *= scale

        n_valid = int(valid.sum())
        n_ideal = ideal_idx.size
        n_nan = int(np.isnan(corrected).sum())
        ideal_pct = n_ideal / max(n_valid, 1) * 100
        nan_pct = n_nan / max(n, 1) * 100
        logging.info(
            f"    {stage_name:<5s}  {n_ideal:>8d}/{n_valid} ({ideal_pct:5.1f}%) fit,  "
            f"{n_nan:>8d}/{n} ({nan_pct:5.1f}%) NaN"
        )
        return corrected, rmse

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)
    grid_01 = np.linspace(0, 1, grid_size)
    extra = (mappability >= min_mappability) if mappability is not None else None
    cor, gc_rmse = _apply_stage(
        reads,
        gc,
        gc_lo,
        gc_hi,
        grid_01,
        seed=seed,
        stage_name="GC",
        extra_ideal_mask=extra,
    )

    if mappability is not None:
        map_lo = np.nanquantile(mappability, doutlier)
        map_hi = np.nanquantile(mappability, 1.0 - doutlier)
        cor, _ = _apply_stage(
            cor,
            mappability,
            map_lo,
            map_hi,
            grid_01,
            seed=seed + 1,
            stage_name="MAP",
        )

    if repliseq is not None:
        repli_finite = np.isfinite(repliseq)
        repli_lo = np.nanquantile(repliseq[repli_finite], doutlier)
        repli_hi = np.nanquantile(repliseq[repli_finite], 1.0 - doutlier)
        grid_repli = np.linspace(repli_lo, repli_hi, grid_size)
        cor, _ = _apply_stage(
            cor,
            repliseq,
            repli_lo,
            repli_hi,
            grid_repli,
            seed=seed + 2,
            stage_name="REPLI",
        )

    return cor.astype(np.float32), gc_rmse


def correct_readcount_quadreg(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    doutlier=0.001,
    min_mappability=0.9,
    eps_quantile=0.01,
):
    """Quadratic median quantile regression bias correction.

    Fits RD ~ GC + GC**2 via median quantile regression on valid bins,
    then divides raw depth by predicted and rescales to preserve median.
    When ``repliseq`` is provided, RT + RT**2 terms are added to the model.
    Mappability is used as a filter (not a covariate).

    Returns
    -------
    np.ndarray
        Corrected read counts, 0.0 where depth is 0. NaN where the correction is
        undefined; the mappability floor is applied by the caller.
    float
        RMSE from the GC fit.
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)

    valid = (reads > 0) & np.isfinite(gc) & (gc >= gc_lo) & (gc <= gc_hi)
    if mappability is not None:
        valid &= mappability >= min_mappability

    if repliseq is not None:
        valid &= np.isfinite(repliseq)

    if repliseq is not None:
        fit_df = pd.DataFrame(
            {"RD": reads[valid], "GC": gc[valid], "RT": repliseq[valid]}
        )
        pred_df = pd.DataFrame({"GC": gc, "RT": np.nan_to_num(repliseq, nan=0.0)})
        formula = "RD ~ GC + I(GC**2) + RT + I(RT**2)"
    else:
        fit_df = pd.DataFrame({"RD": reads[valid], "GC": gc[valid]})
        pred_df = pd.DataFrame({"GC": gc})
        formula = "RD ~ GC + I(GC**2)"

    n_fit = len(fit_df)
    logging.info(f"    MEDIAN  {n_fit:>8d}/{n} fitting bins, formula: {formula}")

    if n_fit < 10:
        logging.warning("    MEDIAN: too few fitting bins; skipping correction")
        return reads.astype(np.float32), 0.0

    res = smf.quantreg(formula, data=fit_df).fit(q=0.5)
    predicted = res.predict(pred_df).to_numpy()

    all_valid = (reads > 0) & np.isfinite(gc)
    n_all_valid = int(all_valid.sum())
    rmse = (
        float(np.sqrt(np.mean((reads[all_valid] - predicted[all_valid]) ** 2)))
        if n_all_valid > 0
        else 0.0
    )

    eps = (
        np.nanquantile(predicted[predicted > 0], eps_quantile)
        if (predicted > 0).any()
        else 1.0
    )
    den = np.clip(np.nan_to_num(predicted, nan=eps), eps, None)

    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = reads / den

    positive = np.isfinite(corrected) & (reads > 0)
    if mappability is not None:
        positive &= mappability >= min_mappability
    if positive.any():
        scale = np.median(reads[positive]) / np.median(corrected[positive])
        corrected[np.isfinite(corrected)] *= scale

    n_nan = int(np.isnan(corrected).sum())
    logging.info(f"    MEDIAN  {n_nan:>8d}/{n} ({n_nan / max(n, 1) * 100:5.1f}%) NaN")

    return corrected.astype(np.float32), rmse


def correct_readcount_by_target_sites(correct_fn, reads, gc, target_sites, **kwargs):
    """Fit and apply *correct_fn* on- and off-target separately.

    Capture makes depth bimodal, and exons are GC-rich, so a pooled fit reads the capture
    split as a GC effect and divides it out. It does so asymmetrically between tumor and
    normal, which have different sets of non-zero windows to fit on. Fitting the two apart
    removes the confounding; each keeps its own depth scale, since every corrector
    rescales to the median of the bins it fit.

    Args:
        correct_fn: ``correct_readcount_lowess`` or ``correct_readcount_quadreg``.
        reads: Raw per-window depth, 1-D.
        gc: Per-window GC fraction, aligned to *reads*.
        target_sites: Bool per window, True where the window overlaps a capture target.
        **kwargs: Passed through; any array-valued entry the length of *reads* is sliced
            to the windows being fit, everything else is passed whole.

    Returns:
        ``(corrected, rmse)``: the corrected depth over every window, and the
        window-count-weighted mean of the two GC RMSEs.
    """
    out = np.full(len(reads), np.nan, dtype=np.float32)
    rmses, weights = [], []
    for label, keep in (("off-target", ~target_sites), ("on-target", target_sites)):
        sub = {
            k: (v[keep] if isinstance(v, np.ndarray) and v.shape == reads.shape else v)
            for k, v in kwargs.items()
        }
        logging.info(f"    {label}: {int(keep.sum())} windows")
        out[keep], rmse = correct_fn(reads[keep], gc[keep], **sub)
        rmses.append(rmse)
        weights.append(int(keep.sum()))
    return out, float(np.average(rmses, weights=weights))


def normalize_library_by_target(mat, target_sites, bin_lengths):
    """Rescale every dataset to a common depth level inside each capture group.

    Hybrid capture gives a library two depth scales, on- and off-target, and their ratio
    is that library's own capture efficiency. A single library-size factor per dataset
    cannot centre both groups at once, so whichever one it centres, the other sits off by
    the efficiency difference; a bb's tumor/normal ratio then slides with its own
    on-target composition. Writing the expected depth of window ``w`` in dataset ``s`` as
    ``alpha_s * kappa_s(g) * c(w) * C_s(w)`` - the dataset's depth scale, its capture
    efficiency in group ``g``, the window-intrinsic efficiency, the copy number -
    rescaling each group to the across-dataset mean of that group divides
    ``alpha_s * kappa_s(g)`` out of every pairwise ratio, leaving ``c(w)`` to cancel as it
    already does. Total sequenced bases does not stand in for ``alpha_s * kappa_s(g)``: it
    says how much was sequenced, not how it split between captured and uncaptured
    sequence, which is the whole quantity at issue.

    What the rescale cannot separate is copy number from capture density: the estimated
    level of group ``g`` carries the genome-average copy number seen through that group's
    windows. The result is a copy-number ratio up to one global constant only while that
    average is the same for both groups, i.e. while copy number is uncorrelated with
    capture-target density genome-wide. Estimating the on-minus-off level gap once
    genome-wide and once over regions believed diploid bounds the departure: the two agree
    when the gap is efficiency rather than copy number.

    Rescales every column, not only the ``bulkWES`` ones: the efficiency term cancels in a
    ratio only when numerator and denominator are both rescaled, so a WES tumor whose RDR
    base is a WGS normal needs the normal rescaled too. A dataset with no capture has the
    same expected level in both groups, so its factors come out near 1. The caller decides
    whether any capture is present at all.

    Args:
        mat: ``(n_bins, n_datasets)`` corrected depth, NaN where undefined.
        target_sites: Bool per bin, True where the bin overlaps a capture target.
        bin_lengths: Per-bin length, aligned to the rows of *mat*.

    Returns:
        ``(scaled, factors)``: the rescaled depth, and ``{group label: (n_datasets,)
        factor}``. A group with no finite bin in a column leaves that column untouched and
        its factor NaN.
    """
    scaled = mat.copy()
    factors = {}
    for label, keep in (("off-target", ~target_sites), ("on-target", target_sites)):
        levels = np.full(mat.shape[1], np.nan, dtype=np.float64)
        for s in range(mat.shape[1]):
            finite = keep & np.isfinite(mat[:, s])
            span = bin_lengths[finite].sum()
            if span > 0:
                levels[s] = float(mat[finite, s] @ bin_lengths[finite] / span)
        ref = np.nanmean(levels) if np.isfinite(levels).any() else np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            fac = np.where(np.isfinite(levels) & (levels > 0), ref / levels, np.nan)
        for s in np.flatnonzero(np.isfinite(fac)):
            scaled[keep, s] = mat[keep, s] * fac[s]
        factors[label] = fac
    return scaled, factors


def compute_gc_rd_stats(mat, gc_vals, labels, n_gc_bins=100):
    """Compute per-label Pearson/Spearman corr(RD, GC) and std of binned median RD.

    Parameters
    ----------
    mat : np.ndarray
        (n_bins, n_samples) depth matrix.
    gc_vals : np.ndarray
        Per-bin GC fraction (same length as mat rows).
    labels : list[str]
        Column labels (sample/dataset_id IDs).
    n_gc_bins : int
        Number of equal-width GC bins in [0, 1].

    Returns
    -------
    gc_corr : dict[str, tuple[float, float]]
        {label: (pearson_r, spearman_r)}
    gc_bin_median_std : dict[str, float]
        {label: std of per-GC-bin median RD (A_GC)}
    """
    gc_bins = np.linspace(0, 1, n_gc_bins + 1)
    gc_corr = {}
    gc_bin_median_std = {}

    for i, label in enumerate(labels):
        v = mat[:, i] if mat.ndim == 2 else mat
        valid = np.isfinite(v) & np.isfinite(gc_vals)

        if valid.sum() > 2:
            pr_val, _ = pearsonr(gc_vals[valid], v[valid])
            sr_val, _ = spearmanr(gc_vals[valid], v[valid])
            gc_corr[label] = (pr_val, sr_val)
        else:
            gc_corr[label] = (np.nan, np.nan)

        gc_bin_idx = np.digitize(gc_vals, gc_bins) - 1
        gc_bin_idx = np.clip(gc_bin_idx, 0, n_gc_bins - 1)
        medians = []
        for b in range(n_gc_bins):
            mask = (gc_bin_idx == b) & np.isfinite(v)
            if mask.any():
                medians.append(np.median(v[mask]))
        gc_bin_median_std[label] = float(np.std(medians)) if medians else np.nan

    return gc_corr, gc_bin_median_std


def compute_depth_statistics(dp_raw, win_df, sample_ids):
    """Compute per-chromosome and whole-genome mean/median depth per sample.

    Returns a DataFrame with columns: SAMPLE, #CHR, mean_depth, median_depth.
    """
    chroms = win_df["#CHR"].to_numpy()
    sorted_chroms = sort_chroms(win_df["#CHR"].unique().tolist())
    rows = []
    for chrom in sorted_chroms:
        mask = chroms == chrom
        for s in range(len(sample_ids)):
            vals = dp_raw[mask, s]
            rows.append(
                [sample_ids[s], chrom, float(np.mean(vals)), float(np.median(vals))]
            )
    for s in range(len(sample_ids)):
        vals = dp_raw[:, s]
        rows.append(
            [sample_ids[s], "TOTAL", float(np.mean(vals)), float(np.median(vals))]
        )
    return pd.DataFrame(rows, columns=["SAMPLE", "#CHR", "mean_depth", "median_depth"])
