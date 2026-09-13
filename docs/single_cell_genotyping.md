# Single-Cell Genotyping

This documentation covers input preparation and result interpretation for single-cell and spatial genotyping using **scRNA**, **scATAC** (incl. 10x Epi Multiome), and **Visium** (`VISIUM`/`VISIUM3prime`) data to run [CalicoST](https://github.com/raphael-group/CalicoST). Refer to the [README](../README.md) for Snakemake pipeline installation and execution instructions.

## Table of Contents
1. [Overview](#overview) <br>
2. [Input](#input) <br>
3. [Output](#output) <br>

## Overview

The rule graph below shows the stages of the single-cell genotyping workflow.

<p align="center">
  <img src="imgs/rulegraph.single_cell_genotyping.png" alt="single_cell_genotyping rule graph" width="440">
</p>

## Input

### Sample file
A sample sheet in JSON format records one or more samples representing patients or cell lines identified by `sample_id`. Each sample can have one or more datasets representing multiple sequencing runs collected from same sample identified by `dataset_id`. Detailed JSON format can be found at [sample_sheet.md](./sample_sheet.md). Here is an example for a 10x Epi Multiome dataset `U1` from patient `HT001`.

```json
{
  "version": 1,
  "samples": [
    {
      "sample_id": "HT001",
      "dataset_id": "U1",
      "assay_type": "scRNA",
      "sample_type": "tumor",
      "files": {
        "alignment": "/data/HT001/multiome/outs/gex_possorted_bam.bam",
        "alignment_index": "/data/HT001/multiome/outs/gex_possorted_bam.bam.bai",
        "barcodes": "/data/HT001/multiome/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "matrix_h5": "/data/HT001/multiome/outs/filtered_feature_bc_matrix.h5"
      }
    },
    {
      "sample_id": "HT001",
      "dataset_id": "U1",
      "assay_type": "scATAC",
      "sample_type": "tumor",
      "files": {
        "alignment": "/data/HT001/multiome/outs/atac_possorted_bam.bam",
        "alignment_index": "/data/HT001/multiome/outs/atac_possorted_bam.bam.bai",
        "barcodes": "/data/HT001/multiome/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "fragments": "/data/HT001/multiome/outs/atac_fragments.tsv.gz"
      }
    }
  ]
}
```

> [!IMPORTANT]
> - (`sample_id`, `dataset_id`) must uniquely define a dataset.
> - A multiome paired dataset share the same `dataset_id`, one `scRNA` and one `scATAC` record.
> - alignment file (`alignment`) must be sorted, and its index file (`alignment_index`) must present.

### Config file

A Snakemake config file is required to specify the runtime configurations. Copy the [template](../resources/templates/config.yaml) and adjust following parameters. detailed description can be found at [reference.md](reference.md#configuration).

1. specify workflow mode, assay types (`assay_types`), and patient informations.

```yaml
workflow_mode: "single_cell_genotyping"
assay_types: ["scRNA", "scATAC"]

sample_id: HT001
sample_file: /path/to/samples.json
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
```

2. specify the paths to reference files, see [../resources/README.md](../resources/README.md) for detailed descriptions for pre-built reference files. Only datasets with `reference_version` recorded in config will be processed.

```yaml
species: human
reference_version: hg38
reference: /path/to/reference.fasta
genome_size: resources/data/hg38.chrom.sizes
region_bed: resources/data/hg38.regions.bed
extremity_tsv: null
gtf_file: /path/to/gencode.v38.annotation.gtf.gz
gene_blacklist_file: resources/data/ig_gene_list.txt
```

> [!TIP]
> Set `extremity_tsv` to a TSV of upstream SV breakpoints and no window or bb will span an
> SV junction; see [`extremity_tsv`](reference.md#file-paths).

3. specify the population SNP panel (`snp_panel`) for germline SNP genotyping. Unlike bulk mode, single-cell genotyping piles up a pseudobulk of all datasets of a modality with [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/) over `snp_panel`. See [snp-panels](../resources/README.md#snp-panels) for details.

```yaml
snp_panel: /path/to/snp_panel.vcf.gz
```

> [!TIP]
> - If a set of confident germline (phased) Het SNPs already exist (e.g., from matched-bulk data), user may specify the path via `het_snp_vcf` and set `het_snp_vcf_phased` to indicate if the VCF file is phased or not. This will skip the germline SNP genotyping (and haplotype phasing if `het_snp_vcf_phased=true`).

4. Use Eagle2 to phase germline SNPs. Genetic map file (`gmap_path`, see [genetic-maps](../resources/README.md#genetic-maps)) and population haplotype panel (`phasing_panel`, see [population-haplotype-panels](../resources/README.md#population-haplotype-panels)) are required.

```yaml
phaser: "eagle"
phasing_panel: /path/to/1kGP_3202_hg38/phasing_panel
gmap_path: /path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
```

5. The final step performs adaptive binning over the fixed bins jointly across all tumor datasets and obtain genomic bin by dataset UMI/ATAC-fragment counts, phased B-allele counts, and total-allele counts. Each value in the minimum-SNP-covering reads parameter (`min_snp_reads`) gives one binning result. We recommend user to set `min_snp_reads` to a list of values and inspect the QC plots at `<qc_dir>/combine_counts.{assay_type}.MSR{msr}.pdf` for varying `min_snp_reads`, then pick the lowest value that gives reliable BAF signals.

```yaml
params_combine_counts:
  min_snp_reads: [50, 100, 150, 200]
  nsnp_multi: 2
```

The workflow also produces counts data with genomic bins with at most `nsnp_multi` SNPs per bin, independent of `min_snp_reads`.

## Output

Refer to [Final bins](reference.md#final-bins) for the full specification of each file:

```text
<out_dir>/
  bb/
    {assay_type}.h5ad                          # gene x cell AnnData (scRNA/VISIUM), MSR-independent
    unit/
      {assay_type}/                            # the un-binned grids, MSR-independent
        snp.{tsv.gz,Tallele.npz,Aallele.npz,Ballele.npz}   # per-SNP allele counts, SNPs x cells
        window.{tsv.gz,Xcount.npz}             # scATAC only: fragments per window per cell
        gene.{tsv.gz,Xcount.npz}               # RNA assays only: UMIs per gene per cell
        barcodes.tsv.gz                        # the matrix column axis
        sample_ids.tsv                         # roster: one row per dataset x assay
    multi_snp/
      {assay_type}/                            # multi-SNP diagnostic groups, MSR-independent
        bb.tsv.gz                              # one row per group, same schema as a bb set
        bb.{Tallele,Aallele,Ballele}.npz       # phased allele counts, groups x cells
    MSR{msr}/                                   # one subdir per min_snp_reads value
      {assay_type}/                             # one subdir per assay; the shared bbs duplicated into each
        bb.tsv.gz                              # bb annotations (shared by every matrix below)
        bb.{Tallele,Aallele,Ballele}.npz       # phased allele counts, bins x cells
        bb.Xcount.npz                          # native counts (scATAC fragments / RNA UMIs)
        barcodes.tsv.gz                        # {BARCODE}_{dataset_id}_{assay_type} per row
        sample_ids.tsv                         # roster: one row per dataset x assay (barcodes.tsv.gz is the column axis)
  qc/
    post_genotype_snps.nonbulk.pdf             # SNP allele frequency by called genotype
    phase_and_concat.{assay_type}.pdf          # SNP allele frequency + depth histogram
    combine_counts.{assay_type}.MSR{msr}.pdf   # binning QC, one per min_snp_reads value
```
