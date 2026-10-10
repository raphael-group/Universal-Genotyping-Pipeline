"""Single-cell: one phased allele matrix over every non-bulk cell, on one shared SNP set.

Last update: 2026-08-11

Inputs and outputs: see docs/reference.md#intermediates.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import hstack, save_npz

from const import ASSAY_TYPE2MODALITY
from io_utils import read_segment_bed, read_VCF, write_sample_ids, write_snp_info
from range_utils import overlaps_any_range
from phase_and_concat_utils import (
    apply_masks_to_df,
    get_mask_by_blacklist,
    get_mask_by_exon,
    get_mask_by_region,
    hstack_replicate_mats,
    interp_pos_ranges,
    map_allele_mat_to_snps,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import plot_allele_freqs, plot_snp_depth
from phasing_utils import apply_phase_to_mat
from feature_utils import annotate_feature_type


##################################################

# inputs, one entry per (assay, replicate) except h5ads (one per RNA assay)
vcf_files = list(snakemake_handle.input["vcfs"])
sample_tsvs = list(snakemake_handle.input["sample_tsvs"])
tot_mtx_files = list(snakemake_handle.input["tot_mtxs"])
ad_mtx_files = list(snakemake_handle.input["ad_mtxs"])
h5ad_files = list(snakemake_handle.input["h5ad_files"])
snp_vcf = snakemake_handle.input["snp_vcf"]
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]
gtf_file = maybe_path(snakemake_handle.input["gtf_file"])
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
assay_types = list(snakemake_handle.params["assay_types"])
dataset_assays = list(snakemake_handle.params["dataset_assays"])
dataset_ids = list(snakemake_handle.params["dataset_ids"])
sample_types = list(snakemake_handle.params["sample_types"])
exon_only = snakemake_handle.params["exon_only"]
run_id = snakemake_handle.params["run_id"]

# outputs; only the QC PDFs are per assay
out_snp_info = snakemake_handle.output["snp_info"]
out_tot_mtx_snp = snakemake_handle.output["tot_mtx_snp"]
out_a_mtx_snp = snakemake_handle.output["a_mtx_snp"]
out_b_mtx_snp = snakemake_handle.output["b_mtx_snp"]
out_all_barcodes = snakemake_handle.output["all_barcodes"]
out_sample_file = snakemake_handle.output["sample_file"]
out_qc_pdf = list(snakemake_handle.output["qc_pdf"])

n_assays = len(assay_types)
rna_assay_types = [at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"]
assert len(h5ad_files) == len(rna_assay_types), (
    f"h5ad_files, {len(h5ad_files)} files for {len(rna_assay_types)} RNA assays"
)
h5ad_by_assay = dict(zip(rna_assay_types, h5ad_files))
assay2dataset_indices = {
    at: [i for i, a in enumerate(dataset_assays) if a == at] for at in assay_types
}
empty = [at for at, ind in assay2dataset_indices.items() if not ind]
assert not empty, f"no replicate for assay(s) {empty}"

logging.info(
    f"phase_and_concat_nonbulk\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#datasets={len(dataset_ids)}\n"
    f"#datasets(per assay)={[len(assay2dataset_indices[at]) for at in assay_types]}"
)

##################################################
# 1. every replicate of every assay onto the shared parent SNP set
snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "SNP VCF, duplicate keys (not bi-allelic)"

barcodes_list = [None] * len(dataset_ids)
tot_by_dataset = [None] * len(dataset_ids)
ad_by_dataset = [None] * len(dataset_ids)
for idx, dataset_id in enumerate(dataset_ids):
    barcodes = pd.read_table(sample_tsvs[idx], sep="\t", header=None, names=["BARCODE"])
    barcodes["BARCODE"] = barcodes["BARCODE"].astype(str)
    barcodes["BARCODE"] += f"_{dataset_id}_{dataset_assays[idx]}"
    barcodes_list[idx] = barcodes
    # cellsnp-lite emits one VCF + two MatrixMarket files per replicate
    tot_by_dataset[idx], ad_by_dataset[idx] = map_allele_mat_to_snps(
        parent_keys,
        read_VCF(vcf_files[idx], addkey=True),
        mmread(tot_mtx_files[idx]).tocsr(),
        mmread(ad_mtx_files[idx]).tocsr(),
        len(barcodes),
    )

tot_list, ref_list, alt_list = [], [], []
barcodes_by_assay, cell_dataset_ids_list = [], []
dataset_ids_by_assay, sample_types_by_assay = [], []
for at in assay_types:
    dataset_indices = assay2dataset_indices[at]
    tot_a, ref_a, alt_a = hstack_replicate_mats(
        [tot_by_dataset[i] for i in dataset_indices],
        [ad_by_dataset[i] for i in dataset_indices],
    )
    tot_list.append(tot_a)
    ref_list.append(ref_a)
    alt_list.append(alt_a)

    all_barcodes = pd.concat(
        [barcodes_list[i] for i in dataset_indices], ignore_index=True
    )
    cell_dataset_ids = np.repeat(
        np.arange(len(dataset_indices), dtype=np.int64),
        [len(barcodes_list[i]) for i in dataset_indices],
    )
    barcodes_by_assay.append(all_barcodes)
    cell_dataset_ids_list.append(cell_dataset_ids)
    dataset_ids_by_assay.append([dataset_ids[i] for i in dataset_indices])
    sample_types_by_assay.append([sample_types[i] for i in dataset_indices])
logging.info(f"#cells(per assay)={[m.shape[1] for m in tot_list]}")

##################################################
# 2. SNP filters, applied once to the shared set
num_snps_before = len(snps)
regions = read_segment_bed(region_bed)

# feature_id (;-joined GTF genes) + feature_type, uniform across all assays
snps = annotate_feature_type(snps, gtf_file)

masks = [
    get_mask_by_region(snps, regions),
    get_mask_by_blacklist(snps, blacklist_bed),
]
if exon_only:
    masks.append(get_mask_by_exon(snps))

# RNA coverage is per assay: RNA reads cover expressed genes, so an RNA assay ignores
# SNPs outside its h5ad feature set. On a shared grid that cannot drop a row for one
# assay, so it zeroes below; only a SNP no assay covers is dropped, which keeps the
# shared set equal to the union of the per-assay sets.
cov_masks = []
for at in assay_types:
    if at not in h5ad_by_assay:
        cov_masks.append(np.ones(len(snps), dtype=bool))
        continue
    adata: sc.AnnData = sc.read_h5ad(h5ad_by_assay[at])
    cov = overlaps_any_range(snps, adata.var)
    logging.info(
        f"{at} feature overlap: {np.sum(cov)}/{len(snps)} "
        f"({np.sum(cov) / len(snps):.3%}); the rest zeroed"
    )
    cov_masks.append(cov)
masks.append(np.logical_or.reduce(cov_masks))

snps, snp_mask = apply_masks_to_df(snps, *masks)
snps = interp_pos_ranges(snps, regions, colname="region_id")
logging.info(f"#SNPs={np.sum(snp_mask)}/{num_snps_before} after filtering")

tot_list = [m[snp_mask, :] for m in tot_list]
ref_list = [m[snp_mask, :] for m in ref_list]
alt_list = [m[snp_mask, :] for m in alt_list]


def zero_uncovered(mat, cov):
    """Zero the rows of a sparse matrix where *cov* is False, keeping its dtype."""
    return mat.multiply(cov[:, None]).tocsr().astype(mat.dtype)


for k, cov in enumerate(c[snp_mask] for c in cov_masks):
    if cov.all():
        continue
    tot_list[k] = zero_uncovered(tot_list[k], cov)
    ref_list[k] = zero_uncovered(ref_list[k], cov)
    alt_list[k] = zero_uncovered(alt_list[k], cov)

phases = snps["PHASE"].to_numpy()
a_list, b_list = [], []
for tot_a, ref_a, alt_a in zip(tot_list, ref_list, alt_list):
    a_a, b_a = apply_phase_to_mat(tot_a, ref_a, alt_a, phases)
    a_list.append(a_a)
    b_list.append(b_a)

##################################################
# 3. QC, one page-set per assay: the cells differ, so the plots cannot merge
for k, assay_type in enumerate(assay_types):
    at_dataset_ids = dataset_ids_by_assay[k]
    at_sample_types = sample_types_by_assay[k]
    at_assay_types = [assay_type] * len(at_dataset_ids)
    cell_dataset_ids = cell_dataset_ids_list[k]

    with PdfPages(out_qc_pdf[k]) as pdf:
        plot_snp_depth(
            tot_list[k],
            at_dataset_ids,
            at_assay_types,
            at_sample_types,
            qc_dir,
            f"{assay_type}.{run_id}",
            ref_mtx=ref_list[k],
            b_mtx=b_list[k],
            is_bulk=False,
            cell_dataset_ids=cell_dataset_ids,
            name_prefix="phase_and_concat",
            pdf=pdf,
            sample_id=sample_id,
        )
        for allele, mat, suffix in (
            ("ref", ref_list[k], ".unphased"),
            ("B", b_list[k], ".phased"),
        ):
            plot_allele_freqs(
                snps,
                at_dataset_ids,
                at_assay_types,
                at_sample_types,
                tot_list[k],
                mat,
                genome_size,
                qc_dir,
                apply_pseudobulk=True,
                allele=allele,
                feature_label="SNP",
                suffix=suffix,
                region_bed=region_bed,
                blacklist_bed=blacklist_bed,
                run_id=run_id,
                sample_id=sample_id,
                pdf=pdf,
                cell_dataset_ids=cell_dataset_ids,
            )

##################################################
# 4. one matrix over every assay's cells; barcodes.tsv.gz is that column axis, while
# sample_ids.tsv is the (dataset_id, assay_type) roster and is NOT column-aligned
logging.info("saving output files")
write_snp_info(snps, out_snp_info)

all_barcodes = pd.concat(barcodes_by_assay, ignore_index=True)
save_npz(out_tot_mtx_snp, hstack(tot_list, format="csr"))
save_npz(out_a_mtx_snp, hstack(a_list, format="csr"))
save_npz(out_b_mtx_snp, hstack(b_list, format="csr"))
all_barcodes.to_csv(out_all_barcodes, sep="\t", header=False, index=False)
write_sample_ids(sample_id, dataset_ids, sample_types, dataset_assays, out_sample_file)
logging.info(
    f"union matrix: {len(snps)} SNPs x {len(all_barcodes)} cells over "
    f"{len(dataset_ids)} (dataset x assay) observations"
)
logging.info("finished joint non-bulk phase_and_concat.")
