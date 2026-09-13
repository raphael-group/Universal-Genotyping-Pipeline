# Copytyping Preprocess

This documentation covers input preparation and result interpretation for the copytyping preprocess mode, which aggregates single-cell / spatial (**scRNA**, **scATAC**, **Visium**) allele and native counts onto a **pre-computed** set of copy-number bbs to run [Copy-typing](https://github.com/raphael-group/Copy-typing). Refer to the [README](../README.md) for Snakemake pipeline installation and execution instructions.

## Table of Contents
1. [Overview](#overview) <br>
2. [Input](#input) <br>
3. [Output](#output) <br>

## Overview

The rule graph below shows the stages of the copytyping preprocess workflow.

<p align="center">
  <img src="imgs/rulegraph.copytyping_preprocess.png" alt="copytyping_preprocess rule graph" width="260">
</p>

## Input

### Sample file
The sample sheet is the same one `single_cell_genotyping` takes: see
[Sample file](single_cell_genotyping.md#sample-file) for the multiome example, and
[sample_sheet.md](./sample_sheet.md) for the schema.

### Config file

A Snakemake config file is required to specify the runtime configurations. Copy the [template](../resources/templates/config.yaml) and adjust following parameters. detailed description can be found at [reference.md](reference.md#configuration).

1. specify workflow mode, assay types (`assay_types`), and patient informations.

```yaml
workflow_mode: "copytyping_preprocess"
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
gtf_file: /path/to/gencode.v38.annotation.gtf.gz
gene_blacklist_file: resources/data/ig_gene_list.txt
```

3. specify the pre-computed **phased** het-SNP VCF (`het_snp_vcf`) and the copy-number bb annotations (`bb_file`). Both are **required** in this mode: genotyping and phasing are skipped, and the counts are aggregated onto the given bbs. A natural source is a prior `bulk_genotyping` run of the same patient (its `phase/phased_het_snps.vcf.gz` and a `bb.tsv.gz`).

```yaml
het_snp_vcf: /path/to/phased_het_snps.vcf.gz
het_snp_vcf_phased: true
bb_file: /path/to/bb.tsv.gz
```

## Output

Refer to [Final bins](reference.md#final-bins) for the full specification of each file:

```text
<out_dir>/
  bb/
    {assay_type}.h5ad                          # gene x cell AnnData (scRNA/VISIUM)
    unit/
      {assay_type}/                            # the un-binned grids, independent of bb_file
        snp.{tsv.gz,Tallele.npz,Aallele.npz,Ballele.npz}   # per-SNP allele counts, SNPs x cells
        window.{tsv.gz,Xcount.npz}             # scATAC only: fragments per window per cell
        gene.{tsv.gz,Xcount.npz}               # RNA assays only: UMIs per gene per cell
        barcodes.tsv.gz                        # the matrix column axis
        sample_ids.tsv                         # roster: one row per dataset x assay
    {assay_type}/                              # per assay, flat (no MSR{msr}/ layer)
      bb.tsv.gz                                # bb annotations (the given bb_file, re-stamped)
      bb.{Xcount,Tallele,Aallele,Ballele}.npz # per-bb native + phased allele counts, bbs x cells
      barcodes.tsv.gz                          # {BARCODE}_{dataset_id}_{assay_type} per row
      sample_ids.tsv                           # roster: one row per dataset x assay (barcodes.tsv.gz is the column axis)
  qc/
    phase_and_concat.{assay_type}.pdf          # SNP allele frequency + depth histogram
    combine_counts_fixed_bins.{assay_type}.pdf # per dataset: SNP BAF, then bb RDR over BAF
```
