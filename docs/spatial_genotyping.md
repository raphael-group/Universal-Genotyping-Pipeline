# Spatial Genotyping

This documentation covers input preparation and result interpretation for spatial genotyping on **Visium HD 3'** data (`assay_type: VISIUM3prime`) to run [CalicoST](https://github.com/raphael-group/CalicoST). Refer to the [README](../README.md) for Snakemake pipeline installation and execution instructions.

## Table of Contents
1. [Overview](#overview) <br>
2. [Preprocessing](#preprocessing) <br>
3. [Input](#input) <br>
4. [Output](#output) <br>

## Overview

```text
FASTQs + images --run_spaceranger.sh--> Space Ranger outs/ (2 um BAM)
                --prepare_visium_hd.py--> X um BAM, matrix, spatial/, sample_sheet.json
                --UGP single_cell_genotyping--> bb/ (X um bins as barcodes)
```

Refer to the rule graph in [single_cell_genotyping.md](single_cell_genotyping.md#overview).

## Preprocessing

### 1. Space Ranger with a BAM

UGP piles up allele counts from the alignment, and Space Ranger writes the BAM only with `--create-bam=true`. See [`resources/scripts/run_spaceranger.sh`](../resources/scripts/run_spaceranger.sh) runs `spaceranger count` (Space Ranger >= 4.0) with that flag:

```sh
bash resources/scripts/run_spaceranger.sh \
    /path/to/spaceranger-4.0.1/spaceranger /path/to/refdata-gex-GRCh38-2024-A \
    <id> <slide>-<area> <fastq_dir> <cytassist_image.tif> <tissue_image.btf> <slide.vlf> <spaceranger_dir>
```

> [!IMPORTANT]
> Use the FASTA and GTF of the transcriptome as UGP's `reference` and `gtf_file`: the BAM coordinates and the matrix gene ids must match them.

Space Ranger writes one BAM, `outs/possorted_genome_bam.bam`, whose `CB` tag is the 2 um barcode `s_002um_<row>_<col>-1`. The 8 um and 16 um binned outputs under `outs/binned_outputs/` have no BAM of their own ([Space Ranger outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)).

### 2. Aggregate 2 um bins into X um bins

[`resources/scripts/prepare_visium_hd.py`](../resources/scripts/prepare_visium_hd.py) builds every UGP input of one slide at an even bin size X from Space Ranger's `outs/`. We recommend `--bin-um 16` to achieve spots with near single-cell resolution.

```sh
python resources/scripts/prepare_visium_hd.py \
    --ranger-outs <spaceranger_dir>/outs --bin-um 16 \
    --out-dir /path/to/ugp_inputs/square_016um \
    --sample-id HT003 --dataset-id U1 --reference-version GRCh38 --threads 16
```

> [!IMPORTANT]
> `barcodes` and `alignment` of a record must be at the same bin size: cellsnp-lite keeps only reads whose `CB` is in `barcodes`.

## Input

### Sample file
`prepare_visium_hd.py` writes the record. Detailed JSON format can be found at [sample_sheet.md](./sample_sheet.md). Several slides of one patient are several records sharing `sample_id`, each with its own `dataset_id`.

```json
{
  "version": 1,
  "samples": [
    {
      "sample_id": "HT003",
      "dataset_id": "U1",
      "assay_type": "VISIUM3prime",
      "sample_type": "tumor",
      "platform": "10x Genomics Visium HD 3prime, square_016um bins",
      "reference_version": "GRCh38",
      "files": {
        "alignment": "/path/to/ugp_inputs/square_016um/possorted_genome_bam.bam",
        "alignment_index": "/path/to/ugp_inputs/square_016um/possorted_genome_bam.bam.bai",
        "barcodes": "/path/to/ugp_inputs/square_016um/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "matrix_h5": "/path/to/ugp_inputs/square_016um/filtered_feature_bc_matrix.h5",
        "tissue_positions": "/path/to/ugp_inputs/square_016um/spatial/tissue_positions.parquet",
        "scalefactors": "/path/to/ugp_inputs/square_016um/spatial/scalefactors_json.json"
      }
    }
  ]
}
```

> [!NOTE]
> `tissue_positions` and `scalefactors` are required for `VISIUM3prime`, but `process_rna_anndata` reads them only together with the tissue images, which it never loads for `VISIUM3prime`. The `.parquet` positions are therefore accepted as is, and `bb/VISIUM3prime.h5ad` carries no spatial coordinates; join them back from `tissue_positions.parquet` by barcode.

### Config file

```yaml
workflow_mode: "single_cell_genotyping"
assay_types: ["VISIUM3prime"]

sample_id: HT003
sample_file: /path/to/ugp_inputs/square_016um/sample_sheet.json
reference_version: hg38
reference: /path/to/refdata-gex-GRCh38-2024-A/fasta/genome.fa
gtf_file: /path/to/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
```

### Rerun UGP with varying spot bin size

SNP genotyping and phasing read every alignment record whatever its `CB`, so `snps/`, `phase/` and `aux/` of one bin size hold for every other bin size of the same slide. A run at a new bin size can copy those directories from an earlier run and mark them complete:

```sh
snakemake --profile profile/ -s workflow/Snakefile --configfile <config.yaml> \
    --directory <new_out_dir> --touch --until concat_and_extract_phased_het_snps build_window_bed
```

Only pileup, allele counting, `process_rna_anndata` and `combine_counts` then rerun.

## Output

The output layout is [single_cell_genotyping.md#output](single_cell_genotyping.md#output) with `{assay_type}` = `VISIUM3prime`. Each barcode is an X um bin, written as `s_<XXX>um_<row>_<col>-1_{dataset_id}_VISIUM3prime` in `barcodes.tsv.gz`.
