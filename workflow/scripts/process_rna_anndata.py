"""One AnnData per RNA-family assay, over all of the sample's datasets.

Last update: 2026-09-02

Inputs:
- barcodes: per-dataset cell barcodes from the sample file
- matrix_h5: per-dataset 10x Ranger feature-barcode matrix
- tissue_positions.csv: spatial spot coordinates
- scalefactors_json.json: spatial scale factors
- tissue_hires_image.png: spatial tissue image, hires
- tissue_lowres_image.png: spatial tissue image, lowres
- gtf_file: gene coordinates, joined on gene_id_colname
- gene_blacklist_file: optional, one gene symbol or id per line
- aux_dir/segment.bed: a gene outside every region is dropped
Outputs:
- bb_dir/{assay}.h5ad: cells x genes, obs_names {raw}_{dataset_id}_{assay_type}
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path, sort_chroms

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
import anndata

from const import SPATIAL_ASSAYS
from io_utils import (
    read_10x_ranger_scRNA,
    read_10x_ranger_spatial,
    read_BED,
    read_barcodes,
    read_GTF,
)
from feature_utils import assign_features_to_ranges

##################################################
# inputs
barcode_files = list(snakemake_handle.input["barcodes"])
matrix_h5_files = list(snakemake_handle.input["matrix_h5"])
spatial_files = list(snakemake_handle.input["spatial_files"])
gtf_file = snakemake_handle.input["gtf_file"]
gene_blacklist_file = maybe_path(snakemake_handle.input["gene_blacklist_file"])
region_bed = snakemake_handle.input["region_bed"]

# parameters
assay_type = snakemake_handle.params["assay_type"]
dataset_ids = list(snakemake_handle.params["dataset_ids"])
# per-dataset spatial/ filenames, aligned with spatial_files
spatial_names = list(snakemake_handle.params["spatial_names"])
gene_id_colname = snakemake_handle.params["gene_id_colname"]
min_frac_barcodes = float(snakemake_handle.params["min_frac_barcodes"])

# outputs
out_h5ad_file = snakemake_handle.output["h5ad_file"]

logging.info(f"prepare rna anndata, assay_type={assay_type}, dataset_ids={dataset_ids}")

# snakemake flattens nested `input:` lists, so regroup per dataset by name count
n_spatial = [len(names) for names in spatial_names]
assert len(spatial_files) == sum(n_spatial), (
    f"spatial_files, {len(spatial_files)} paths for {sum(n_spatial)} names"
)
bounds = np.cumsum([0] + n_spatial)
spatial_paths = [spatial_files[i:j] for i, j in zip(bounds[:-1], bounds[1:])]

adatas = {}
for idx, dataset_id in enumerate(dataset_ids):
    logging.info(f"process {assay_type}-{dataset_id}")
    barcodes = pd.Index(read_barcodes(barcode_files[idx])).astype(str)

    matrix_h5 = matrix_h5_files[idx]
    if assay_type in SPATIAL_ASSAYS:
        names = spatial_names[idx]
        adata = read_10x_ranger_spatial(
            matrix_h5,
            names,
            spatial_paths[idx],
            library_id=dataset_id,
            assay_type=assay_type,
            # squidpy doesn't support load images from 3' data yet.
            load_images=assay_type == "VISIUM",
        )
    else:
        adata = read_10x_ranger_scRNA(matrix_h5)

    adata.obs_names = adata.obs_names.astype(str)
    adata = adata[adata.obs_names.isin(barcodes), :].copy()
    adata.obs_names = adata.obs_names + f"_{dataset_id}_{assay_type}"
    adatas[dataset_id] = adata
    logging.info(f"#barcodes={adata.n_obs}, #features={adata.n_vars}")

if len(adatas) > 1:
    adata = anndata.concat(
        adatas,
        join="outer",  # union of var (genes)
        label="dataset_id",
        merge="same",
        uns_merge="unique",
        fill_value=0,
    )
else:
    adata = adatas[dataset_ids[0]]
adata.X = adata.X.tocsr()
num_total_barcodes = adata.n_obs
logging.info(f"#concat barcodes={num_total_barcodes}, #union features={adata.n_vars}")

genes_gtf = read_GTF(gtf_file, ("gene",), id_col=gene_id_colname)["gene"]
logging.info(f"#genes in the GTF={len(genes_gtf)}")

var_coords = adata.var.merge(
    genes_gtf, how="left", on=gene_id_colname, validate="m:1", sort=False
)
var_coords.index = adata.var.index
var_coords["pseudobulk_counts"] = np.asarray(adata.X.sum(axis=0)).ravel()
num_genes = len(var_coords)
umis = var_coords["pseudobulk_counts"].to_numpy()

# per-gene statistics are independent, so the masks compose and X is subset once
keep = ~var_coords["START"].isna().to_numpy()
if not keep.all():
    logging.warning(
        f"#genes not found in reference GTF file={num_genes - keep.sum()}/{num_genes}"
    )

drop = keep & (umis == 0)
logging.info(f"remove #{drop.sum()}/{num_genes} genes with zero pseudobulk_counts")
keep &= ~drop

if gene_blacklist_file is not None:
    gene_blacklist = (
        pd.read_table(gene_blacklist_file, header=None).iloc[:, 0].to_numpy()
    )
    # NB: var_names carry make_unique suffixes, gene_symbol is the raw 10x symbol
    drop = keep & (
        np.isin(var_coords.index, gene_blacklist)
        | np.isin(var_coords["gene_symbol"], gene_blacklist)
    )
    logging.info(
        f"remove #{drop.sum()}/{num_genes} genes based on {gene_blacklist_file}"
    )
    keep &= ~drop

if assay_type in SPATIAL_ASSAYS:
    min_expressed_barcodes = round(min_frac_barcodes * num_total_barcodes)
    drop = keep & (adata.X.getnnz(axis=0) < min_expressed_barcodes)
    count_ratio = umis[keep & ~drop].sum() / umis[keep].sum()
    logging.info(
        f"remove #{drop.sum()}/{num_genes} genes expressed in "
        f"<{min_expressed_barcodes}/{num_total_barcodes} barcodes "
        f"(min_frac_barcodes={min_frac_barcodes}), "
        f"keeping {100.0 * count_ratio:.2f}% of UMIs"
    )
    keep &= ~drop

adata = adata[:, keep].copy()
adata.var = var_coords.loc[keep, :].copy()
adata.var["#CHR"] = adata.var["#CHR"].astype(str)
adata.var["START"] = adata.var["START"].astype(int)
adata.var["END"] = adata.var["END"].astype(int)
logging.info(f"#genes after filtering={adata.n_vars}/{num_genes}")

regions = read_BED(region_bed, col_id="region_id")
adata = assign_features_to_ranges(adata, regions, assay_type)

chroms = sort_chroms(adata.var["#CHR"].unique().tolist())
adata.var["#CHR"] = pd.Categorical(adata.var["#CHR"], categories=chroms, ordered=True)

assert adata.var_names.is_unique, "AnnData, var_names are not unique"
sort_index = adata.var.sort_values(by=["#CHR", "START"]).index
adata = adata[:, sort_index].copy()

adata.write_h5ad(out_h5ad_file, compression="gzip")

logging.info(f"final processed {assay_type} AnnData")
logging.info(f"final #obs={adata.n_obs}, #vars={adata.n_vars}")
