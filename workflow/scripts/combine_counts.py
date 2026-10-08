"""Bulk: SNP-informed adaptive binning over all bulk assays, plus depth, RDR and reads.

Last update: 2026-09-25

Inputs and outputs: see docs/reference.md#unit-level and docs/reference.md#final-bins.

Read starts are additive, so every level sums them from the windows; depth, a per-base
mean, is length-weighted instead. Both cover every dataset, in sample_ids.tsv order.

A bb closes when its SNP reads meet min_snp_reads in every tumor AND its read starts meet
min_total_reads in every dataset. Inside a clonal-LOH region the SNP criterion is dropped
(no germline het survives there) and the bb's allele counts are stamped: B=0 in a tumor,
T=the reads the hets it would have carried at the neutral rate, so its weight in a
downstream beta-binomial matches a measured bb of the same span.

Under detect_loh_tumor_cell_line those regions are called here from the collapse in
het-SNP density per tile, over the run's shared SNP set and every bulk dataset's
assayable span; see script_utils/loh_utils.py. The genotyping step never sees them,
since one site's counts cannot tell a gHET inside LOH from a gHOM.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, log_hist, log_ratios, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from segmentation_utils import build_adaptive_bins, sum_features_to_bbs
from range_utils import (
    assign_pos_to_range,
    assign_range_to_range,
    merge_ranges_to_clusters,
)
from feature_utils import explode_feature_ids, merge_feature_ids
from loh_utils import build_loh_tiles, fit_snp_rate, loh_intervals, viterbi_loh
from io_utils import (
    read_mosdepth_bed,
    read_snp_mats,
    read_window_bed,
    write_bb_file,
)
from combine_counts_utils import (
    summarize_read_depth_bb,
    summarize_rdr_bb,
    tumor_observation_indices,
)
from phase_em import correct_bin_phases
from phasing_utils import (
    apply_phase_to_mat,
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_between_bbs,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_combine_counts import (
    plot_loh_density,
    plot_rdr_baf,
    plot_rdr_baf_2d,
    plot_segmentation_qc,
)

##################################################

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
dp_corrected_file = snakemake_handle.input["dp_corrected"]
rdcount_files = list(snakemake_handle.input["rdcount_files"])
window_bed = snakemake_handle.input["window_bed"]
sample_file = snakemake_handle.input["sample_file"]
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])
genome_size = snakemake_handle.input["genome_size"]

# parameters
sample_id = snakemake_handle.params["sample_id"]
assay_types = list(snakemake_handle.params["assay_types"])
dp_dataset_ids = list(snakemake_handle.params["dataset_ids"])
dp_dataset_assays = list(snakemake_handle.params["dataset_assays"])
chroms = list(snakemake_handle.params["chroms"])

phase_em = bool(snakemake_handle.params["phase_em"])
phase_em_tau = float(snakemake_handle.params["phase_em_tau"])
phase_em_min_llr = float(snakemake_handle.params["phase_em_min_llr"])
phase_em_n_grid = int(snakemake_handle.params["phase_em_n_grid"])
phase_em_grid_eps = float(snakemake_handle.params["phase_em_grid_eps"])
phase_em_n_restarts = int(snakemake_handle.params["phase_em_n_restarts"])
phase_em_max_iter = int(snakemake_handle.params["phase_em_max_iter"])
phase_em_tol = float(snakemake_handle.params["phase_em_tol"])
phase_em_min_snps = int(snakemake_handle.params["phase_em_min_snps"])

gene_aware_binning = bool(snakemake_handle.params["gene_aware_binning"])
phaseset_aware_binning = bool(snakemake_handle.params["phaseset_aware_binning"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_total_reads = int(snakemake_handle.params["min_total_reads"])

detect_loh_tumor_cell_line = bool(snakemake_handle.params["detect_loh_tumor_cell_line"])
loh_tile_size = int(snakemake_handle.params["loh_tile_size"])
loh_rate_ratio = float(snakemake_handle.params["loh_rate_ratio"])
loh_breakpoint_rate = float(snakemake_handle.params["loh_breakpoint_rate"])
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
nsnp_multi = int(snakemake_handle.params["nsnp_multi"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])

# outputs
out_bb_file = list(snakemake_handle.output["bb_file"])
out_tot_mtx_bb = list(snakemake_handle.output["tot_mtx_bb"])
out_a_mtx_bb = list(snakemake_handle.output["a_mtx_bb"])
out_b_mtx_bb = list(snakemake_handle.output["b_mtx_bb"])
out_dp_mtx_bb = list(snakemake_handle.output["dp_mtx_bb"])
out_rdr_mtx_bb = list(snakemake_handle.output["rdr_mtx_bb"])
out_rdcount_mtx_bb = list(snakemake_handle.output["rdcount_mtx_bb"])
out_sample_file = list(snakemake_handle.output["sample_file"])
out_unit_snp_file = snakemake_handle.output["unit_snp_file"]
out_unit_tot_mtx = snakemake_handle.output["unit_tot_mtx"]
out_unit_a_mtx = snakemake_handle.output["unit_a_mtx"]
out_unit_b_mtx = snakemake_handle.output["unit_b_mtx"]
out_unit_window_file = snakemake_handle.output["unit_window_file"]
out_unit_dp_mtx = snakemake_handle.output["unit_dp_mtx"]
out_unit_rdcount_mtx = snakemake_handle.output["unit_rdcount_mtx"]
out_unit_sample_file = snakemake_handle.output["unit_sample_file"]
out_multi_bb_file = snakemake_handle.output["multi_bb_file"]
out_multi_tot_mtx = snakemake_handle.output["multi_tot_mtx"]
out_multi_a_mtx = snakemake_handle.output["multi_a_mtx"]
out_multi_b_mtx = snakemake_handle.output["multi_b_mtx"]
out_multi_dp_mtx = snakemake_handle.output["multi_dp_mtx"]
out_multi_rdr_mtx = snakemake_handle.output["multi_rdr_mtx"]
out_multi_rdcount_mtx = snakemake_handle.output["multi_rdcount_mtx"]
out_multi_sample_file = snakemake_handle.output["multi_sample_file"]
out_qc_stats_pdf = list(snakemake_handle.output["qc_stats_pdf"])
out_qc_1d2d_pdf = list(snakemake_handle.output["qc_1d2d_pdf"])
out_flip_tsv_bb = list(snakemake_handle.output["flip_tsv_bb"])
out_loh_pdf = snakemake_handle.output["loh_pdf"] if detect_loh_tumor_cell_line else None

##################################################
# load inputs
sample_df = pd.read_table(sample_file)
bin_df = read_window_bed(window_bed, chroms=chroms)
snps, tot_mtx, a_mtx, b_mtx = read_snp_mats(
    snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp, mat_dtype=np.int32
)
dp_corrected = np.load(dp_corrected_file)["mat"]
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None

##################################################
# observation order, labels and RDR bases
num_datasets = len(sample_df)
dataset_assays = sample_df["assay_type"].to_numpy()
dataset_ids = sample_df["dataset_id"].tolist()
tumor_dataset_indices = tumor_observation_indices(sample_df)
tumor_dataset_ids = [dataset_ids[c] for c in tumor_dataset_indices]
tumor_assays = [dataset_assays[c] for c in tumor_dataset_indices]

assert dp_dataset_ids == dataset_ids and dp_dataset_assays == list(dataset_assays), (
    f"depth columns {list(zip(dp_dataset_ids, dp_dataset_assays))} do not match "
    f"{sample_file}: {list(zip(dataset_ids, dataset_assays))}"
)
assert dp_corrected.shape == (len(bin_df), num_datasets), (
    f"{dp_corrected_file}: shape {dp_corrected.shape}, expected "
    f"({len(bin_df)}, {num_datasets}) from {window_bed}"
)
assert len(rdcount_files) == num_datasets, (
    f"{len(rdcount_files)} read-start files for {num_datasets} datasets"
)

# read starts, joined onto the window grid by coordinate as rd_correct does for depth
join_keys = ["#CHR", "START", "END"]
rdcount = np.zeros((len(bin_df), num_datasets), dtype=np.int32)
for i, (dataset_id, rdc_file) in enumerate(zip(dataset_ids, rdcount_files)):
    counts = bin_df[join_keys].merge(
        read_mosdepth_bed(rdc_file, value_col="COUNT"),
        on=join_keys,
        how="left",
        sort=False,
    )
    n_missing = int(counts["COUNT"].isna().sum())
    assert n_missing == 0, (
        f"{dataset_id}: {n_missing}/{len(bin_df)} windows absent from {rdc_file}"
    )
    rdcount[:, i] = counts["COUNT"].to_numpy(dtype=np.int32)
logging.info(
    f"read starts: {rdcount.sum(axis=0).tolist()} over {len(bin_df)} windows, "
    f"columns {dataset_ids}"
)

get_rdr_base_dataset_id = {}
if "rdr_base_dataset_id" in sample_df.columns:
    obs_of = {rid: i for i, rid in enumerate(dataset_ids)}
    for i in tumor_dataset_indices:
        base = sample_df["rdr_base_dataset_id"].iloc[i]
        if pd.notna(base):
            get_rdr_base_dataset_id[i] = obs_of[base]
is_tumor_col = np.zeros(num_datasets, dtype=bool)
is_tumor_col[tumor_dataset_indices] = True
if detect_loh_tumor_cell_line and len(tumor_dataset_indices) > 1:
    logging.warning(
        f"{len(tumor_dataset_indices)} tumor columns {tumor_dataset_ids} share one "
        f"clonal-LOH segmentation, called from the run's shared het-SNP set"
    )

tumor_base_dataset_ids = [
    dataset_ids[get_rdr_base_dataset_id[c]] if c in get_rdr_base_dataset_id else None
    for c in tumor_dataset_indices
]

logging.info(
    f"combine_counts\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#SNPs={len(snps)}\n"
    f"#windows={len(bin_df)}\n"
    f"#datasets={num_datasets}\n"
    f"#tumor_datasets={len(tumor_dataset_indices)}"
)

##################################################
# assign SNPs to windows, then drop the misses from the matrices too
tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_dataset_indices], dtype=np.float64)
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
snp_spans = (snps["END"] - snps["START"]).to_numpy()
log_ratios(
    "SNPs outside every window", len(off_idx), len(snps), snp_spans[off_idx], "snp"
)
if len(off_idx):
    log_hist(tot_tumor[off_idx].sum(axis=1), "depth of SNPs outside every bin")
log_hist(
    snps_binned.groupby("bin_id").size().reindex(range(len(bin_df)), fill_value=0),
    "SNPs per window",
)
keep_snps = np.ones(len(snps), dtype=bool)
keep_snps[off_idx] = False
tot_mtx, a_mtx, b_mtx = tot_mtx[keep_snps], a_mtx[keep_snps], b_mtx[keep_snps]
tot_tumor = np.ascontiguousarray(tot_tumor[keep_snps])

##################################################
# clonal LOH: the het-SNP density collapses there, so those windows bin on read starts
# alone and their allele counts are stamped rather than measured
if detect_loh_tumor_cell_line:
    # a window with no assayable span emits nothing, so a coverage hole cannot read as
    # depletion; depth enters here only, never as a copy-number signal
    assayable = np.isfinite(dp_corrected).all(axis=1) & (rdcount.sum(axis=1) > 0)
    tiles = build_loh_tiles(
        bin_df,
        snps_binned["bin_id"].to_numpy(),
        assayable,
        loh_tile_size,
        cluster_cols=("seg_id",) if "seg_id" in bin_df.columns else (),
    )
    exposure_mb = tiles["exposure_bp"].to_numpy() / 1e6
    rate_ref, nb_size = fit_snp_rate(tiles["n_het"].to_numpy(), exposure_mb)
    rate_loh = loh_rate_ratio * rate_ref
    tiles["is_loh"] = viterbi_loh(
        tiles["n_het"].to_numpy(),
        exposure_mb,
        tiles["mid"].to_numpy(),
        pd.factorize(tiles["region_id"])[0],
        rate_ref,
        rate_loh,
        nb_size,
        loh_breakpoint_rate,
    )
    logging.info(
        f"clonal-LOH decode: {rate_ref:.1f} het/Mb neutral (NB size {nb_size:.2f}) vs "
        f"{rate_loh:.1f} LOH (ratio {loh_rate_ratio}), "
        f"{int(tiles['is_loh'].sum())}/{len(tiles)} tiles LOH"
    )
    loh_df = loh_intervals(tiles, tiles["is_loh"].to_numpy())
    if len(loh_df):
        flagged, _ = assign_range_to_range(
            bin_df[["#CHR", "START", "END"]],
            loh_df.assign(loh_region=np.arange(len(loh_df))),
            "loh_region",
            rule="midpoint",
        )
        bin_df["is_loh"] = flagged["loh_region"].notna().to_numpy()
    else:
        bin_df["is_loh"] = False

    win_is_loh = bin_df["is_loh"].to_numpy()
    snp_in_loh = win_is_loh[snps_binned["bin_id"].to_numpy()]
    log_ratios(
        f"windows in a clonal-LOH region ({len(loh_df)} intervals)",
        int(win_is_loh.sum()),
        len(bin_df),
        (bin_df["END"] - bin_df["START"]).to_numpy()[win_is_loh],
        "window",
    )
    # the residual het calls inside LOH: sequencing error the cutoff let through, kept
    # in the unit-level table but overwritten by the bb stamp
    log_ratios(
        "SNPs inside an LOH window (their counts are overwritten by the stamp)",
        int(snp_in_loh.sum()),
        len(snps_binned),
        np.zeros(int(snp_in_loh.sum())),
        "snp",
    )
    with PdfPages(out_loh_pdf) as loh_pdf:
        plot_loh_density(
            tiles,
            rate_ref,
            rate_loh,
            loh_df,
            genome_size,
            sample_id,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=loh_pdf,
        )
    logging.info(f"saved clonal-LOH QC to {out_loh_pdf}")
else:
    bin_df["is_loh"] = False

# one key per maximal run, so a bb never spans an LOH boundary in either direction
bin_df["loh_id"] = (bin_df["is_loh"] != bin_df["is_loh"].shift()).cumsum()
win_is_loh = bin_df["is_loh"].to_numpy()
# the LOH key is defined on windows; a SNP inherits its window's
snps_binned["loh_id"] = bin_df["loh_id"].to_numpy()[snps_binned["bin_id"].to_numpy()]
snps_binned["is_loh"] = win_is_loh[snps_binned["bin_id"].to_numpy()]

##################################################
# unit level: the two grids the binning consumes, before any merge
snps_binned.drop(
    columns=["bin_id", "loh_id"] + ([] if detect_loh_tumor_cell_line else ["is_loh"])
).to_csv(out_unit_snp_file, sep="\t", index=False)
np.savez_compressed(out_unit_tot_mtx, mat=tot_mtx)
np.savez_compressed(out_unit_a_mtx, mat=a_mtx)
np.savez_compressed(out_unit_b_mtx, mat=b_mtx)
bin_df.drop(
    columns=["bin_id", "loh_id"] + ([] if detect_loh_tumor_cell_line else ["is_loh"])
).to_csv(out_unit_window_file, sep="\t", index=False)
np.savez_compressed(out_unit_dp_mtx, mat=dp_corrected)
np.savez_compressed(out_unit_rdcount_mtx, mat=rdcount)
sample_df.to_csv(out_unit_sample_file, sep="\t", index=False)
logging.info(
    f"unit level: {len(snps_binned)} SNPs to {out_unit_snp_file}, "
    f"{len(bin_df)} windows to {out_unit_window_file}"
)

##################################################
# adaptive segmentation bounderies
cluster_cols = ["region_id", "seg_id", "loh_id"]

if phaseset_aware_binning and "PS" in snps.columns:
    assert snps["PS"].notna().all(), "SNP file, `PS` column has NaNs"
    cluster_cols.append("PS")
    modal = snps_binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
    bin_df["PS"] = bin_df["bin_id"].map(modal).ffill().bfill().fillna(1)

if gene_aware_binning:
    gene_spans = (
        explode_feature_ids(snps_binned, cols=["bin_id"])
        .groupby("feature_id")["bin_id"]
        .agg(["min", "max"])
    )
    bin_df["gene_cluster"] = merge_ranges_to_clusters(
        len(bin_df), zip(gene_spans["min"].to_numpy(), gene_spans["max"].to_numpy() + 1)
    )

logging.info(
    f"binning cluster keys={cluster_cols}, gene_aware_binning={gene_aware_binning}, "
    f"phaseset_aware_binning={phaseset_aware_binning}"
)

##################################################
# multi-SNP groups: nsnp_multi SNPs each, independent of the binning sweep
multi_cols = ["#CHR", "START", "END", "region_id"] + (
    ["seg_id"] if "seg_id" in snps_binned.columns else []
)
multi_bins = snps_binned[multi_cols].reset_index(drop=True)
multi_bins["bin_id"] = np.arange(len(multi_bins))
multi_snps_in, off_multi = assign_pos_to_range(
    snps_binned, multi_bins, ref_id="bin_id", dropna=True
)
keep_multi = np.ones(len(multi_bins), dtype=bool)
keep_multi[off_multi] = False
multi_bin_spans = (multi_bins["END"] - multi_bins["START"]).to_numpy()
log_ratios(
    "SNPs with no per-SNP range, dropped from the multi-SNP grouping",
    len(off_multi),
    len(multi_bins),
    multi_bin_spans[off_multi],
    "snp",
)
if len(off_multi):
    log_hist(tot_tumor[off_multi].sum(axis=1), "depth of SNPs with no per-SNP range")

multi_bbs, snps_multi = build_adaptive_bins(
    multi_bins,
    multi_snps_in,
    np.ascontiguousarray(tot_tumor[keep_multi]),
    0,
    nsnp_multi,
    cluster_cols=[c for c in ("region_id", "seg_id") if c in multi_bins.columns],
    gene_aware=False,
)
num_multi = len(multi_bbs)
multi_ids = snps_multi["bb_id"].to_numpy()
multi_tot = sum_features_to_bbs(tot_mtx[keep_multi], multi_ids, num_multi)
multi_a = sum_features_to_bbs(a_mtx[keep_multi], multi_ids, num_multi)
multi_b = sum_features_to_bbs(b_mtx[keep_multi], multi_ids, num_multi)

# the groups tile each region, so every window falls in one; depth follows by midpoint
multi_spans = (
    multi_bbs.groupby("bb_id", sort=True)
    .agg(
        **{"#CHR": ("#CHR", "first"), "START": ("START", "min"), "END": ("END", "max")}
    )
    .reset_index()
)
mapped, off_win = assign_range_to_range(
    bin_df[["#CHR", "START", "END"]], multi_spans, "bb_id", rule="midpoint"
)
win_multi = (
    pd.to_numeric(mapped["bb_id"], errors="coerce").fillna(-1).to_numpy(np.int64)
)
# an LOH window carries no measured allele signal, so it joins no multi-SNP group
inside = (win_multi >= 0) & ~win_is_loh
win_spans = (bin_df["END"] - bin_df["START"]).to_numpy()
logging.info(f"multi-SNP groups: {num_multi} over {int(keep_multi.sum())} SNPs")
log_ratios(
    "bins outside every multi-SNP group",
    int((~inside).sum()),
    len(bin_df),
    win_spans[~inside],
    "bin",
)
multi_bin_df = bin_df.loc[inside].assign(bb_id=win_multi[inside])
multi_dp, multi_bases = summarize_read_depth_bb(
    multi_bin_df, dp_corrected[inside], num_multi, dataset_ids=dataset_ids
)
multi_rdr = summarize_rdr_bb(
    multi_bin_df,
    dp_corrected[inside],
    multi_dp,
    tumor_dataset_indices,
    get_rdr_base_dataset_id,
    dataset_ids,
)
multi_rdcount = sum_features_to_bbs(rdcount[inside], win_multi[inside], num_multi)

multi_baf = np.divide(
    multi_b,
    multi_tot,
    where=multi_tot > 0,
    out=np.full_like(multi_b, np.nan, dtype=np.float32),
)
multi_nan_baf = np.isnan(multi_baf).any(axis=1)
multi_nan_dp = np.isnan(multi_dp).any(axis=1)
multi_nan_rdr = np.isnan(multi_rdr).any(axis=1)
multi_nan = multi_nan_baf | multi_nan_dp | multi_nan_rdr
multi_spans_bp = multi_bbs["BLOCKSIZE"].to_numpy()
log_ratios(
    "multi-SNP NaN row filter, dropped",
    int(multi_nan.sum()),
    num_multi,
    multi_spans_bp[multi_nan],
    "bb",
)
for _cause, _m in (
    ("BAF", multi_nan_baf),
    ("depth", multi_nan_dp),
    ("RDR", multi_nan_rdr),
):
    log_ratios(
        f"NaN by cause (overlapping) {_cause}",
        int(_m.sum()),
        num_multi,
        multi_spans_bp[_m],
        "bb",
        prefix="  ",
    )
valid_multi = ~multi_nan
multi_bbs = multi_bbs.loc[valid_multi].reset_index(drop=True)
multi_tot, multi_a, multi_b = (
    multi_tot[valid_multi],
    multi_a[valid_multi],
    multi_b[valid_multi],
)
multi_dp, multi_rdr = multi_dp[valid_multi], multi_rdr[valid_multi]
multi_rdcount = multi_rdcount[valid_multi]
old_to_new_multi = {old: new for new, old in enumerate(np.where(valid_multi)[0])}
snps_multi_valid = snps_multi[snps_multi["bb_id"].isin(old_to_new_multi)].copy()
snps_multi_valid["bb_id"] = snps_multi_valid["bb_id"].map(old_to_new_multi)
multi_bbs["bb_id"] = np.arange(len(multi_bbs))

if genetic_map is not None:
    multi_bbs["switchprobs"] = estimate_switchprobs_cM(
        interp_cM_between_bbs(
            multi_bbs, snps_multi_valid, genetic_map, bb_id_col="bb_id"
        ),
        nu=nu,
        min_switchprob=min_switchprob,
    )
else:
    multi_bbs["switchprobs"] = estimate_switchprobs_PS(multi_bbs, switchprob_ps)
multi_bbs["feature_id"] = (
    multi_bbs["bb_id"]
    .map(snps_multi_valid.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
    .fillna("intergenic")
)

write_bb_file(multi_bbs, out_multi_bb_file)
np.savez_compressed(out_multi_tot_mtx, mat=multi_tot)
np.savez_compressed(out_multi_a_mtx, mat=multi_a)
np.savez_compressed(out_multi_b_mtx, mat=multi_b)
np.savez_compressed(out_multi_dp_mtx, mat=multi_dp)
np.savez_compressed(out_multi_rdr_mtx, mat=multi_rdr)
np.savez_compressed(out_multi_rdcount_mtx, mat=multi_rdcount)
sample_df.to_csv(out_multi_sample_file, sep="\t", index=False)
logging.info(f"wrote {len(multi_bbs)} multi-SNP groups to {out_multi_bb_file}")

##################################################
# one adaptive binning per min_snp_reads, on the shared fixed bins
for (
    msr,
    out_flip_tsv,
    out_bb,
    out_tot,
    out_a,
    out_b,
    out_dp,
    out_rdr,
    out_rdcount,
    out_samp,
    out_stats_pdf,
    out_1d2d_pdf,
) in zip(
    msr_list,
    out_flip_tsv_bb,
    out_bb_file,
    out_tot_mtx_bb,
    out_a_mtx_bb,
    out_b_mtx_bb,
    out_dp_mtx_bb,
    out_rdr_mtx_bb,
    out_rdcount_mtx_bb,
    out_sample_file,
    out_qc_stats_pdf,
    out_qc_1d2d_pdf,
):
    min_snp_reads_vec = np.full(len(tumor_dataset_indices), msr, dtype=np.float64)
    bbs, snps_bb = build_adaptive_bins(
        bin_df,
        snps_binned,
        tot_tumor,
        min_snp_reads_vec,
        min_snp_per_bin,
        cluster_cols=cluster_cols,
        gene_aware=gene_aware_binning,
        bin_total_reads=rdcount,
        min_total_reads=min_total_reads,
        is_loh_col="is_loh",
    )
    num_bbs = len(bbs)
    bb_spans = bbs["BLOCKSIZE"].to_numpy()
    loh_bb = np.zeros(num_bbs, dtype=bool)
    loh_bb[bin_df["bb_id"].to_numpy()[win_is_loh]] = True
    # carried on the frame so the NaN filter below keeps it aligned to the rows
    bbs["is_loh"] = loh_bb
    empty_bb = ((bbs["#SNPS"] == 0).to_numpy()) & ~loh_bb
    n_empty = int(empty_bb.sum())
    logging.info(f"MSR={msr}: {num_bbs} bbs from {len(bin_df)} bins")
    log_ratios(
        "SNP-free bbs (all-zero allele rows, dropped below)",
        n_empty,
        num_bbs,
        bb_spans[empty_bb],
        "bb",
    )
    if detect_loh_tumor_cell_line:
        log_ratios(
            "clonal-LOH bbs (binned on read starts, allele counts stamped)",
            int(loh_bb.sum()),
            num_bbs,
            bb_spans[loh_bb],
            "bb",
        )

    bb_ids = snps_bb["bb_id"].to_numpy()
    # a stamped bb has no measured orientation, so its residual hets are not fitted
    em_bb_ids = np.where(snps_bb["is_loh"].to_numpy(), -1, bb_ids)
    flips = (
        correct_bin_phases(
            snps_bb.assign(bb_id=em_bb_ids),
            a_mtx,
            b_mtx,
            tumor_dataset_indices,
            phase_em_tau,
            phase_em_min_llr,
            n_grid=phase_em_n_grid,
            grid_eps=phase_em_grid_eps,
            n_restarts=phase_em_n_restarts,
            max_iter=phase_em_max_iter,
            tol=phase_em_tol,
            min_snps=phase_em_min_snps,
        )
        if phase_em and (em_bb_ids >= 0).any()
        else np.zeros(len(snps_bb), dtype=np.int8)
    )
    a_msr, b_msr = apply_phase_to_mat(tot_mtx, b_mtx, a_mtx, 1 - flips.astype(np.int64))
    a_mtx_bb = sum_features_to_bbs(a_msr, bb_ids, num_bbs)
    b_mtx_bb = sum_features_to_bbs(b_msr, bb_ids, num_bbs)
    tot_mtx_bb = sum_features_to_bbs(tot_mtx, bb_ids, num_bbs)
    # the sites a downstream re-phasing of the SNP VCF has to swap
    snps_binned.loc[flips == 1, ["#CHR", "POS"]].to_csv(
        out_flip_tsv, sep="\t", index=False, header=False
    )
    if phase_em:
        bb_flipped = np.zeros(num_bbs, dtype=bool)
        bb_flipped[bb_ids[flips == 1]] = True
        log_ratios(
            f"MSR={msr}: bbs re-oriented by phase_em",
            int(bb_flipped.sum()),
            num_bbs,
            bb_spans[bb_flipped],
            "bb",
        )
        log_ratios(
            "SNPs re-oriented inside their bb",
            int(flips.sum()),
            len(flips),
            np.zeros(int(flips.sum())),
            "snp",
        )

    logging.info("aggregating corrected fixed-bin depth into bbs")
    bb_dp, bb_bases = summarize_read_depth_bb(
        bin_df, dp_corrected, num_bbs, dataset_ids=dataset_ids
    )

    # LOH bbs hold no callable het, so their allele counts are stamped: the tumor's minor
    # haplotype is absent (B = 0), a diploid normal is still het (B = T/2). T is the
    # counterfactual total - the reads the hets a neutral bb of this span would have
    # carried - not the bare depth, so a downstream beta-binomial weights an LOH bb on
    # the same scale as a measured one
    if loh_bb.any():
        stamp_dp = (
            pd.DataFrame(dp_corrected)
            .groupby(bin_df["bb_id"].to_numpy(), sort=True)
            .median()
            .reindex(range(num_bbs))
            .to_numpy()
        )
        bb_span_mb = (
            np.bincount(
                bin_df["bb_id"].to_numpy(),
                weights=(bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.float64),
                minlength=num_bbs,
            )
            / 1e6
        )
        # rate_ref is per Mb of assayable span, so pair it with the summed window span
        expected_het = rate_ref * bb_span_mb[loh_bb]
        stamp_T = np.maximum(
            np.rint(np.nan_to_num(stamp_dp[loh_bb]) * expected_het[:, None]), 1.0
        )
        stamp_B = np.where(is_tumor_col[None, :], 0.0, stamp_T // 2)
        tot_mtx_bb[loh_bb] = stamp_T
        b_mtx_bb[loh_bb] = stamp_B
        a_mtx_bb[loh_bb] = stamp_T - stamp_B
        logging.info(
            f"stamped {int(loh_bb.sum())} clonal-LOH bbs: T={rate_ref:.0f} het/Mb x span "
            f"x median depth (median {np.median(stamp_T):.0f}, median "
            f"{np.median(expected_het):.0f} het/bb), "
            f"B=0 for {[dataset_ids[c] for c in tumor_dataset_indices]}, T//2 otherwise"
        )

    baf_mtx_bb = np.divide(
        b_mtx_bb,
        tot_mtx_bb,
        where=tot_mtx_bb > 0,
        out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
    )
    bb_rdcount = sum_features_to_bbs(rdcount, bin_df["bb_id"].to_numpy(), num_bbs)

    logging.info(
        f"compute bb RDR, {len(get_rdr_base_dataset_id)}/{len(tumor_dataset_indices)} tumors with RDR base"
    )
    bb_rdr = summarize_rdr_bb(
        bin_df,
        dp_corrected,
        bb_dp,
        tumor_dataset_indices,
        get_rdr_base_dataset_id,
        dataset_ids,
    )

    # drop bbs carrying a NaN, before the QC plots, so the PDF shows what is written
    nan_baf = np.isnan(baf_mtx_bb).any(axis=1)
    nan_dp = np.isnan(bb_dp).any(axis=1)
    nan_rdr = np.isnan(bb_rdr).any(axis=1)
    nan_mask = nan_baf | nan_dp | nan_rdr
    n_nan_rows = int(nan_mask.sum())
    n_valid = num_bbs - n_nan_rows
    log_ratios("NaN row filter, dropped", n_nan_rows, num_bbs, bb_spans[nan_mask], "bb")
    log_ratios("NaN row filter, kept", n_valid, num_bbs, bb_spans[~nan_mask], "bb")
    for _cause, _m in (("BAF", nan_baf), ("depth", nan_dp), ("RDR", nan_rdr)):
        log_ratios(
            f"NaN by cause (overlapping) {_cause}",
            int(_m.sum()),
            num_bbs,
            bb_spans[_m],
            "bb",
            prefix="  ",
        )
    _zero_read = nan_baf & ~empty_bb
    log_ratios("BAF no SNP", n_empty, num_bbs, bb_spans[empty_bb], "bb", prefix="    ")
    log_ratios(
        "BAF zero-read column",
        int(_zero_read.sum()),
        num_bbs,
        bb_spans[_zero_read],
        "bb",
        prefix="    ",
    )

    if n_nan_rows > 0:
        valid = ~nan_mask
        bbs = bbs.loc[valid].reset_index(drop=True)
        tot_mtx_bb = tot_mtx_bb[valid]
        a_mtx_bb = a_mtx_bb[valid]
        b_mtx_bb = b_mtx_bb[valid]
        baf_mtx_bb = baf_mtx_bb[valid]
        bb_dp = bb_dp[valid]
        bb_bases = bb_bases[valid]
        bb_rdr = bb_rdr[valid]
        bb_rdcount = bb_rdcount[valid]
    num_bbs = len(bbs)

    kept_bb_ids = np.where(~nan_mask)[0]
    old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
    snps_valid = snps_bb[snps_bb["bb_id"].isin(old_to_new)].copy()
    snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
    bbs["bb_id"] = np.arange(num_bbs)

    if genetic_map is not None:
        dist_cms = interp_cM_between_bbs(
            bbs, snps_valid, genetic_map, bb_id_col="bb_id"
        )
        bbs["switchprobs"] = estimate_switchprobs_cM(
            dist_cms,
            nu=nu,
            min_switchprob=min_switchprob,
        )
    else:
        bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_valid.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

    depth_tumor = bb_dp[:, tumor_dataset_indices]
    depth_normal = np.full_like(depth_tumor, np.nan, dtype=float)
    for j, c in enumerate(tumor_dataset_indices):
        if c in get_rdr_base_dataset_id:
            depth_normal[:, j] = bb_dp[:, get_rdr_base_dataset_id[c]]

    baf_tumor = baf_mtx_bb[:, tumor_dataset_indices]

    with PdfPages(out_stats_pdf) as pdf:
        plot_segmentation_qc(
            bbs,
            sample_df,
            bb_bases,
            bb_rdcount,
            b_mtx_bb,
            tot_mtx_bb,
            pdf=pdf,
            sample_id=sample_id,
            is_loh=bbs["is_loh"].to_numpy() if detect_loh_tumor_cell_line else None,
        )
    logging.info(f"saved bb-statistics QC PDF to {out_stats_pdf}")

    with PdfPages(out_1d2d_pdf) as pdf:
        plot_rdr_baf(
            bbs,
            bb_rdr,
            baf_tumor,
            depth_tumor,
            depth_normal,
            sample_id,
            tumor_dataset_ids,
            tumor_assays,
            tumor_base_dataset_ids,
            genome_size,
            out_1d2d_pdf,
            feature_label="bb",
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=pdf,
            is_loh=bbs["is_loh"].to_numpy() if detect_loh_tumor_cell_line else None,
        )
        plot_rdr_baf_2d(
            bb_rdr,
            baf_tumor,
            sample_id,
            tumor_dataset_ids,
            tumor_assays,
            pdf=pdf,
        )
    logging.info(f"saved genome-wide QC PDF to {out_1d2d_pdf}")

    write_bb_file(
        bbs if detect_loh_tumor_cell_line else bbs.drop(columns="is_loh"), out_bb
    )
    np.savez_compressed(out_tot, mat=tot_mtx_bb)
    np.savez_compressed(out_a, mat=a_mtx_bb)
    np.savez_compressed(out_b, mat=b_mtx_bb)
    np.savez_compressed(out_dp, mat=bb_dp)
    np.savez_compressed(out_rdr, mat=bb_rdr)
    np.savez_compressed(out_rdcount, mat=bb_rdcount)
    sample_df.to_csv(out_samp, sep="\t", index=False)
    logging.info(f"MSR={msr}: wrote {len(bbs)} bins")

logging.info("finished combine_counts.")
