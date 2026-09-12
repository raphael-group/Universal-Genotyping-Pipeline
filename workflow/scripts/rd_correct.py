"""Bulk: GC/mappability/replication-timing bias correction of per-fixed-bin depth.

Last update: 2026-08-12

Inputs:
- pileup_dir/{assay}/out_mosdepth/{dataset_id}.regions.bed.gz: per-dataset per-bin depth
- aux_dir/windows.bed.gz: the shared fixed bins with GC, MAP, REPLI, region_id
- aux_dir/window.target.npz: optional per-window capture-target fraction; when present,
  a bulkWES dataset is corrected on- and off-target separately and every dataset is then
  rescaled to a common depth level within each group
- genome_size, region_bed, blacklist_bed: QC plot axis and shading
Outputs:
- pileup_dir/bulk/window.raw.dp.npz: raw mosdepth depth, windows x all bulk datasets,
  unmasked and uncorrected
- pileup_dir/bulk/window.dp.npz: corrected depth, windows x all bulk datasets,
  0.0 where depth is 0 and NaN below min_mappability or where the fit is undefined
- pileup_dir/bulk/depth_statistics.tsv: per-dataset depth summary
- qc_dir/rd_correction.bulk.pdf: depth scatter before/after plus covariate KDE

Both matrices are float32, row-aligned to the window BED filtered to `chroms`, and
column-ordered as `params.dataset_ids`, which is the `SAMPLE` order of
depth_statistics.tsv. Neither npz stores column labels.
"""

import logging

snakemake_handle = snakemake

from utils import (
    set_omp_threads,
    setup_logging,
    maybe_path,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np

from io_utils import read_mosdepth_bed, read_window_bed
from rd_correct_utils import (
    compute_depth_statistics,
    compute_gc_rd_stats,
    correct_readcount_by_target_sites,
    correct_readcount_lowess,
    correct_readcount_quadreg,
    normalize_library_by_target,
)
from plot_count_reads import plot_rd_1d_scatter, plot_rd_2d_kde

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages


# inputs
mosdepth_files = list(snakemake_handle.input["mosdepth_files"])
window_bed = snakemake_handle.input["window_bed"]
window_target = maybe_path(snakemake_handle.input["window_target"])
genome_size = snakemake_handle.input["genome_size"]
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])

# parameters
sample_id = snakemake_handle.params["sample_id"]
dataset_ids = list(snakemake_handle.params["dataset_ids"])
dataset_assays = list(snakemake_handle.params["dataset_assays"])
sample_types = list(snakemake_handle.params["sample_types"])
chroms = list(snakemake_handle.params["chroms"])
samplesize = int(snakemake_handle.params["samplesize"])
routlier = float(snakemake_handle.params["routlier"])
doutlier = float(snakemake_handle.params["doutlier"])
min_mappability = float(snakemake_handle.params["min_mappability"])
gc_correct = bool(snakemake_handle.params["gc_correct"])
rd_correct_method = snakemake_handle.params["rd_correct_method"]
rt_correct = bool(snakemake_handle.params["rt_correct"])

# outputs
out_depth_stats = snakemake_handle.output["depth_stats"]
out_dp_raw = snakemake_handle.output["dp_raw"]
out_dp_corrected = snakemake_handle.output["dp_corrected"]
out_qc_pdf = snakemake_handle.output["qc_pdf"]

n_samples = len(dataset_ids)
sample_ids = [f"{sample_id}_{dataset_id}" for dataset_id in dataset_ids]
join_keys = ["#CHR", "START", "END"]

logging.info(
    f"rd_correct: {n_samples} bulk datasets across assays={sorted(set(dataset_assays))}, "
    f"{len(chroms)} chroms"
)

logging.info("load the shared fixed bins and the per-dataset mosdepth depth")
bin_df = read_window_bed(window_bed, chroms=chroms, keep_covariates=True)
assert "GC" in bin_df.columns, (
    f"window_bed, missing `GC` column: {bin_df.columns.tolist()}"
)
n_bins = len(bin_df)
logging.info(f"{n_bins} fixed bins across {len(chroms)} chromosomes")

dp_raw = np.zeros((n_bins, n_samples), dtype=np.float32)
for i, (dataset_id, mos_file) in enumerate(zip(dataset_ids, mosdepth_files)):
    mos_df = read_mosdepth_bed(mos_file)
    depth = bin_df[join_keys].merge(mos_df, on=join_keys, how="left", sort=False)
    n_missing = int(depth["DEPTH"].isna().sum())
    assert n_missing == 0, (
        f"{dataset_id}: {n_missing}/{n_bins} fixed bins absent from {mos_file}"
    )
    dp_raw[:, i] = depth["DEPTH"].to_numpy(dtype=np.float32)

np.savez_compressed(out_dp_raw, mat=dp_raw)
logging.info(f"wrote raw depth to {out_dp_raw}")

depth_stats = compute_depth_statistics(dp_raw, bin_df, sample_ids)
depth_stats.to_csv(out_depth_stats, sep="\t", index=False)
logging.info(f"wrote depth statistics to {out_depth_stats}")
for _, row in depth_stats[depth_stats["#CHR"] == "TOTAL"].iterrows():
    logging.info(
        f"  {row['SAMPLE']}: mean={row['mean_depth']:.2f}, median={row['median_depth']:.2f}"
    )

gc_vals = bin_df["GC"].to_numpy()

rd_raw_ylim = max(np.nanquantile(dp_raw, 0.99), 1.0) * 1.1
gc_corr_before, gc_std_before = compute_gc_rd_stats(dp_raw, gc_vals, dataset_ids)

logging.info(f"{n_bins} fixed bins for bias correction")

map_vals = bin_df["MAP"].to_numpy() if "MAP" in bin_df.columns else None
repli_vals = (
    bin_df["REPLI"].to_numpy(dtype=np.float64)
    if rt_correct and "REPLI" in bin_df.columns
    else None
)
if repli_vals is not None:
    n_repli_finite = int(np.isfinite(repli_vals).sum())
    logging.info(
        f"REPLI column: {n_repli_finite}/{n_bins} "
        f"({n_repli_finite / max(n_bins, 1) * 100:.1f}%) finite"
    )
else:
    logging.info("no REPLI column; skipping replication timing correction")

target_sites = None
if window_target is not None:
    target_frac = np.load(window_target)["mat"]
    assert len(target_frac) == n_bins, (
        f"{window_target}: {len(target_frac)} rows for {n_bins} windows in {window_bed}"
    )
    target_sites = target_frac > 0
    logging.info(
        f"capture targets: {int(target_sites.sum())}/{n_bins} on-target windows; "
        "bulkWES datasets are corrected on- and off-target separately"
    )

gc_rmse_list = None
if gc_correct:
    dp_corrected = np.zeros_like(dp_raw, dtype=np.float32)
    gc_rmse_list = []
    if rd_correct_method == "median":
        correct_readcount = correct_readcount_quadreg
        extra_kwargs = {}
    else:
        correct_readcount = correct_readcount_lowess
        extra_kwargs = {"samplesize": samplesize, "routlier": routlier}

    logging.info(f"applying {correct_readcount.__name__} per sample")
    for i, dataset_id in enumerate(dataset_ids):
        by_target = target_sites is not None and dataset_assays[i] == "bulkWES"
        logging.info(
            f"correcting {dataset_id}"
            + (" on- and off-target separately" if by_target else "")
        )
        kwargs = dict(
            mappability=map_vals,
            repliseq=repli_vals,
            doutlier=doutlier,
            min_mappability=min_mappability,
            **extra_kwargs,
        )
        if by_target:
            dp_corrected[:, i], gc_rmse = correct_readcount_by_target_sites(
                correct_readcount, dp_raw[:, i], gc_vals, target_sites, **kwargs
            )
        else:
            dp_corrected[:, i], gc_rmse = correct_readcount(
                dp_raw[:, i], gc_vals, **kwargs
            )
        gc_rmse_list.append(gc_rmse)
else:
    logging.info("gc_correct=False; skipping bias correction")
    dp_corrected = dp_raw.copy()

rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0) * 1.1

with PdfPages(out_qc_pdf) as pdf:
    plot_rd_1d_scatter(
        bin_df,
        dp_raw,
        dp_corrected,
        dataset_ids,
        dataset_assays,
        sample_types,
        genome_size,
        pdf,
        ylim_before=rd_raw_ylim,
        ylim_after=rd_ylim,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
    )
    plot_rd_2d_kde(
        gc_vals,
        dp_raw,
        dp_corrected,
        dataset_ids,
        dataset_assays,
        sample_types,
        pdf,
        gc_rmse=gc_rmse_list,
        mappability=map_vals,
        repliseq=repli_vals,
    )

if map_vals is not None:
    low_map = map_vals < min_mappability
    dp_corrected[low_map, :] = np.nan
    logging.info(
        f"mappability filter: {int(low_map.sum())}/{n_bins} fixed bins below "
        f"{min_mappability}, NaN for every dataset"
    )

if target_sites is not None and "bulkWES" in dataset_assays:
    logging.info(
        "library-size normalization per capture group (on- and off-target apart), "
        "every bulk dataset in the run"
    )
    bin_lengths = (bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.float64)
    dp_corrected, target_factors = normalize_library_by_target(
        dp_corrected, target_sites, bin_lengths
    )
    for label, fac in target_factors.items():
        pretty = ", ".join(f"{d}={f:.4f}" for d, f in zip(dataset_ids, fac))
        logging.info(f"    {label}: {pretty}")

for i, dataset_id in enumerate(dataset_ids):
    n_nan = int(np.isnan(dp_corrected[:, i]).sum())
    logging.info(
        f"  {dataset_id}: {n_nan}/{n_bins} fixed bins NaN after correction "
        f"({n_nan / max(n_bins, 1) * 100:.1f}%)"
    )

np.savez_compressed(out_dp_corrected, mat=dp_corrected)
logging.info(f"wrote corrected depth to {out_dp_corrected}")

logging.info("finished rd_correct.")
