"""Post-process called SNPs into the run's per-chromosome het/hom-alt VCFs.

Two independent axes. ``params.source`` says which VCFs hold the counts:
- ``bulk``: the bcftools calls, one alignment set, AD/DP per site.
- ``pseudobulk``: cellsnp-lite DP/AD/OTH, summed over the per-modality base calls.

``params.genotyping`` says how the counts become a genotype:
- ``passthrough``: the calls are germline genotypes already, so the outputs symlink to
  them. Bulk only, and only when no genotyped dataset is a tumor.
- ``vaf_cutoff``: threshold depth, minor-allele reads and VAF. The tumor-only default,
  and the only option off bulk since cellsnp-lite emits no GT of its own.

Last update: 2026-10-12

Inputs and outputs: see docs/reference.md#intermediates.

Notes:
- both sources report depth as cellsnp-lite does, ``REF + ALT`` with the other-allele
  reads held apart in OTH (cellsnp-lite ``src/csp.h``: DP is "total counts for ALT and
  REF"), so ``DP - AD`` is the REF count and one genotyping rule serves both
"""

import logging

snakemake_handle = snakemake

from utils import (
    add_chr_prefix,
    maybe_list,
    set_omp_threads,
    setup_logging,
    strip_chr_prefix,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import (
    read_chrom_sizes,
    read_VCF,
    symlink_files,
    write_VCF,
)
from plot_alleles import plot_allele_freqs

GT_DTYPE = "<U3"

REFINE_PLOT_MAX_SNPS = (
    300_000  # QC scatter resolves no more than this; sampled above it
)
# blue vs red is het against hom, the call this step exists to make; grey is no-call
GT_PLOT_COLORS = {"0/0": "red", "0/1": "blue", "1/1": "red", "./.": "gray"}
# VCF fixed fields this step does not compute, and what to write when the input lacks them
VCF_DEFAULT_FALLBACKS = {"ID": ".", "QUAL": ".", "FILTER": "PASS"}


def write_per_chrom(
    snps, chrom_col, ad, dp, oth, gt, chroms, snp_vcfs, input_nochr, genome_size
):
    """Assemble the output VCF records and write one bgzipped, indexed file per contig.

    ``ID``/``QUAL``/``FILTER`` come from *snps* wherever it carries them: this step
    re-genotypes rather than re-calls, so the caller's own fields survive it.

    Args:
        snps: Rows to write, carrying ``POS``, ``REF``, ``ALT`` and *chrom_col*.
        chrom_col: Column of *snps* holding the chr-prefixed contig.
        ad, dp, oth: ALT, total and other-allele counts, aligned to *snps*.
        gt: Genotype string per row, aligned to *snps*.
        chroms: Contigs, chr-prefixed, aligned to *snp_vcfs*.
        snp_vcfs: Output paths.
        input_nochr: Spell output contigs without the ``chr`` prefix.
        genome_size: Chrom sizes TSV, for the ``##contig`` header.
    """
    out = snps[["POS", "REF", "ALT"]].copy()
    out["#CHROM"] = snps[chrom_col].to_numpy()
    for col, placeholder in VCF_DEFAULT_FALLBACKS.items():
        out[col] = snps[col].to_numpy() if col in snps.columns else placeholder
    out["INFO"] = (
        "AD="
        + pd.Series(np.asarray(ad), index=out.index).astype(str)
        + ";DP="
        + pd.Series(np.asarray(dp), index=out.index).astype(str)
        + ";OTH="
        + pd.Series(np.asarray(oth), index=out.index).astype(str)
    )
    out["FORMAT"] = "GT"
    out["SAMPLE"] = gt

    chrom_sizes = {
        f"chr{strip_chr_prefix(k)}": v for k, v in read_chrom_sizes(genome_size).items()
    }
    by_chrom = out.groupby("#CHROM", sort=False)
    for chrom, snp_vcf in zip(chroms, snp_vcfs):
        rows = by_chrom.get_group(chrom) if chrom in by_chrom.groups else out.iloc[:0]
        write_VCF(
            rows,
            snp_vcf,
            strip_chr_prefix(chrom) if input_nochr else chrom,
            chrom_length=chrom_sizes[chrom],
        )


def read_counts_bulk(raw_snp_vcfs):
    """Concatenate the per-chromosome bcftools calls into one counted frame.

    Args:
        raw_snp_vcfs: The caller's per-chromosome VCFs; empty ones are skipped.

    Returns:
        ``(snps, alt, depth, oth)``. *snps* carries the caller's own ``GT`` and a
        chr-prefixed ``#CHR``; *depth* is ``REF + ALT``, matching cellsnp-lite's ``DP``,
        so the other-allele reads live only in *oth*.

    Raises:
        AssertionError: Every input VCF is empty.
    """
    frames = []
    for raw_snp_vcf in raw_snp_vcfs:
        snps = read_VCF(
            raw_snp_vcf,
            addchr=False,
            snps_presorted=True,
            read_AD=True,
            required_cols=["GT", "AD", "DP"],
        )
        if snps is None:
            logging.warning(f"{raw_snp_vcf}: no records")
            continue
        frames.append(snps)
    assert frames, "bulk: every input VCF is empty"
    snps = pd.concat(frames, ignore_index=True)
    snps["#CHR"] = add_chr_prefix(snps["#CHROM"].astype(str))

    ref = snps["REF_COUNT"].to_numpy()
    alt = snps["ALT_COUNT"].to_numpy()
    # NB: FORMAT/DP counts every base, cellsnp-lite's DP only REF+ALT; match the latter
    total = (
        pd.to_numeric(snps["DP"], errors="coerce").fillna(0).to_numpy().astype(np.int64)
    )
    oth = np.maximum(total - ref - alt, 0)
    return snps, alt, ref + alt, oth


def read_counts_pseudobulk(raw_snp_vcfs, modalities, chroms):
    """Sum the per-modality cellsnp-lite DP/AD/OTH onto one row per site.

    Sites are keyed on ``#CHROM_POS``; a modality missing a site contributes zero. Rows
    whose key repeats within or across modalities are dropped, since the sum would double
    count them.

    Args:
        raw_snp_vcfs: One cellSNP.base.vcf.gz per entry of *modalities*, in that order.
        modalities: Modality labels, for the log lines.
        chroms: Contigs to keep; rows on any other contig are dropped.

    Returns:
        ``(snps, alt, depth, oth)`` plus the pre-genotyping site count, as
        ``(snps, alt, depth, oth, n_sites)``.
    """
    KEY = ["#CHROM", "POS", "REF", "ALT"]
    CNT = ["DP", "AD", "OTH"]
    CARRY = list(VCF_DEFAULT_FALLBACKS)
    raw_snps_list = []
    for idx, modality in enumerate(modalities):
        raw_snps = read_VCF(raw_snp_vcfs[idx], addkey=True, required_cols=KEY + CNT)
        for cnt in CNT:
            raw_snps[f"{cnt}{idx}"] = raw_snps[cnt].astype(np.int64)

        raw_snps = raw_snps[raw_snps["#CHROM"].astype(str).isin(chroms)]

        n_dup_rows = raw_snps.duplicated(subset="KEY", keep=False).sum()
        if n_dup_rows > 0:
            n_dup_keys = raw_snps.loc[
                raw_snps.duplicated(subset="KEY", keep=False), "KEY"
            ].nunique()
            logging.warning(
                f"{modality} have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
            )
            logging.warning("drop duplicated rows.")
            raw_snps = raw_snps.drop_duplicates(subset="KEY", keep="first").reset_index(
                drop=True
            )
        raw_snps_list.append(raw_snps)

    base_snps = pd.concat(
        [
            raw_snps[
                ["KEY", "#CHROM", "POS", "REF", "ALT"]
                + [c for c in CARRY if c in raw_snps.columns]
            ]
            for raw_snps in raw_snps_list
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"], keep="first")

    dup_mask = base_snps.duplicated(subset="KEY", keep=False)
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows > 0:
        n_dup_keys = int(base_snps.loc[dup_mask, "KEY"].drop_duplicates().shape[0])
        logging.warning(
            f"merged df have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
        )
        logging.warning("drop duplicated rows.")
        base_snps = base_snps.loc[~dup_mask, :]
    base_snps = base_snps.reset_index(drop=True)
    n_sites = len(base_snps)

    base_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort")
    for cnt in CNT:
        base_snps[cnt] = 0
    for idx, raw_snps in enumerate(raw_snps_list):
        right_cols = ["KEY"] + [f"{cnt}{idx}" for cnt in CNT]
        base_snps = pd.merge(
            left=base_snps,
            right=raw_snps[right_cols],
            on="KEY",
            how="left",
            validate="one_to_one",
            sort=False,
        )
        for cnt in CNT:
            base_snps[cnt] += base_snps[f"{cnt}{idx}"].fillna(0).astype(np.int64)
        base_snps = base_snps.drop(columns=[f"{cnt}{idx}" for cnt in CNT])

    base_snps["#CHR"] = base_snps["#CHROM"].astype(str)
    return (
        base_snps,
        base_snps["AD"].to_numpy(),
        base_snps["DP"].to_numpy(),
        base_snps["OTH"].to_numpy(),
        n_sites,
    )


def call_vaf_cutoff(alt, depth, min_dp, min_hom_dp, min_het_ad, min_vaf_thres):
    """Genotype each site by thresholding its depth, minor-allele reads and VAF.

    The minor-allele side is ``depth - alt``, so a site needs *min_het_ad* on each of
    the two alleles and a VAF inside ``[min_vaf_thres, 1 - min_vaf_thres]`` to be het.
    Hom-alt is the strict ``alt == depth``; hom-ref never survives the caller's own
    filters, so it is recorded for QC but never kept.

    A het call rests on reads of both alleles, which a homozygote shows only through
    sequencing error. A hom call rests on the other allele being absent, which a het at
    depth d shows by chance with probability ``f**d`` (``f`` its ALT fraction, 0.5 when
    balanced, higher under allelic imbalance). Hom calls therefore take their own floor,
    ``max(min_dp, min_hom_dp)``; an all-ALT site below it is a no-call, not a het.

    Args:
        alt: ALT read count per site.
        depth: Total read depth per site.
        min_dp: Depth floor for every call; below it a site is a no-call.
        min_hom_dp: Depth floor for a hom-alt or hom-ref call, on top of *min_dp*.
        min_het_ad: Minimum reads on each allele for a het call.
        min_vaf_thres: Half-width of the het VAF band.

    Returns:
        Array of ``0/1`` / ``1/1`` / ``0/0`` / ``./.`` per site.
    """
    ref_side = depth - alt
    allele_ratio = pd.Series(alt) / pd.Series(depth)
    called = depth >= min_dp
    hom_called = called & (depth >= min_hom_dp)
    is_het = (
        called
        & (np.minimum(alt, ref_side) >= min_het_ad)
        & allele_ratio.between(min_vaf_thres, 1 - min_vaf_thres).to_numpy()
    )
    is_hom_alt = hom_called & (alt == depth)

    gt = np.full(len(depth), "./.", dtype=GT_DTYPE)
    gt[hom_called & (alt == 0)] = "0/0"
    gt[is_hom_alt] = "1/1"
    gt[is_het] = "0/1"
    return gt


params = snakemake_handle.params
source = params["source"]
genotyping = params["genotyping"]
chroms = list(params["chroms"])
input_nochr = params["input_nochr"]
min_dp = int(params["min_dp"])

raw_snp_vcfs = maybe_list(snakemake_handle.input["raw_snp_vcfs"])
snp_vcfs = maybe_list(snakemake_handle.output["snp_vcfs"])
genome_size = snakemake_handle.input["genome_size"]

logging.info(f"start post_genotype_snps, source={source}, genotyping={genotyping}")

# =============================================================================
# read the counts: one frame, one convention, whichever caller produced them
# =============================================================================
if source == "bulk":
    assert len(raw_snp_vcfs) == len(snp_vcfs), (
        f"bulk: {len(raw_snp_vcfs)} input VCF(s) for {len(snp_vcfs)} output(s)"
    )
    snps, alt_count, depth, oth = read_counts_bulk(raw_snp_vcfs)
    nsnps_before_genotyping = len(snps)
else:
    assert genotyping == "vaf_cutoff", (
        f"source={source} genotypes from counts, got genotyping={genotyping!r}"
    )
    snps, alt_count, depth, oth, nsnps_before_genotyping = read_counts_pseudobulk(
        raw_snp_vcfs, list(params["modalities"]), chroms
    )
ref_count = depth - alt_count

# =============================================================================
# genotype: passthrough keeps the caller's GT, vaf_cutoff overwrites it
# =============================================================================
if genotyping == "vaf_cutoff":
    min_hom_dp = int(params["min_hom_dp"])
    min_het_ad = int(params["min_het_ad"])
    min_vaf_thres = float(params["min_vaf_thres"])
    logging.info(
        f"vaf_cutoff: min_dp={min_dp}, min_hom_dp={min_hom_dp}, "
        f"min_het_ad={min_het_ad}, min_vaf_thres={min_vaf_thres}"
    )
    gt = call_vaf_cutoff(
        alt_count, depth, min_dp, min_hom_dp, min_het_ad, min_vaf_thres
    )
else:
    assert genotyping == "passthrough", f"unknown genotyping={genotyping!r}"
    assert source == "bulk", "passthrough needs a caller's GT, which only bulk has"
    gt = snps["GT"].to_numpy()

# =============================================================================
# keep, write, and hand the same points to the QC plot
# =============================================================================
genotyped = snps.assign(
    REF_COUNT=ref_count, ALT_COUNT=alt_count, DP=depth, OTH=oth, GT=gt
)

if genotyping == "passthrough":
    logging.info(f"linking {len(raw_snp_vcfs)} VCF(s) and indexes")
    symlink_files(raw_snp_vcfs, snp_vcfs)
    symlink_files(
        maybe_list(snakemake_handle.input["raw_snp_vcfs_tbi"]),
        maybe_list(snakemake_handle.output["snp_vcfs_tbi"]),
    )
else:
    filter_nz_OTH = bool(params["filter_nz_OTH"])
    filter_hom_ALT = bool(params["filter_hom_ALT"])
    is_het = gt == "0/1"
    # NB: hom-ref never survives the keep below, so it is not genotyped
    keep = is_het if filter_hom_ALT else (is_het | (gt == "1/1"))
    if filter_nz_OTH:
        logging.info("SNPs with nonzero OTHs are filtered.")
        keep &= oth == 0
    logging.info(f"#nz-OTH SNPs={int((oth > 0).sum())}/{nsnps_before_genotyping}")
    logging.info(f"#kept SNPs={int(keep.sum())}/{nsnps_before_genotyping}")
    write_per_chrom(
        genotyped.loc[keep],
        "#CHR",
        alt_count[keep],
        depth[keep],
        oth[keep],
        gt[keep],
        chroms,
        snp_vcfs,
        input_nochr,
        genome_size,
    )

if len(genotyped) > REFINE_PLOT_MAX_SNPS:
    sampled = np.random.default_rng(0).choice(
        len(genotyped), REFINE_PLOT_MAX_SNPS, replace=False
    )
    logging.info(
        f"QC plot: sampled {REFINE_PLOT_MAX_SNPS} of {len(genotyped)} sites, "
        f"the scatter would not resolve more"
    )
    genotyped = genotyped.iloc[np.sort(sampled)].reset_index(drop=True)

ref_counts = genotyped["REF_COUNT"].to_numpy()
depths = ref_counts + genotyped["ALT_COUNT"].to_numpy()

# NB: the pass-through branch carries the caller's GT, which may be phased
genotyped["GT_PLOT"] = genotyped["GT"].astype(str).str.replace("|", "/", regex=False)
with PdfPages(snakemake_handle.output["qc_pdf"]) as qc_pdf:
    plot_allele_freqs(
        genotyped,
        [],
        [],
        [],
        depths.reshape(-1, 1),
        ref_counts.reshape(-1, 1),
        genome_size,
        params["qc_dir"],
        apply_pseudobulk=True,
        cell_dataset_ids=None,
        allele="ref",
        feature_label="SNP",
        snp_groups=genotyped["GT_PLOT"].to_numpy(),
        group_colors=GT_PLOT_COLORS,
        run_id=params["run_id"],
        sample_id=params["sample_id"],
        name_prefix="post_genotype_snps",
        pdf=qc_pdf,
    )
logging.info("finished post_genotype_snps")
