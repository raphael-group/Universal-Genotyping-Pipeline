# Sample File Format

Sample file is a JSON-format configuration file contains dataset records over multiple patients and patient-specific datasets. Template can be found at here: [templates](../resources/templates/).


## Schema

```
{
  "version": 1,
  "samples": [
    {
      "sample_id":             <str>,
      "dataset_id":            <str>,
      "rdr_base_dataset_id"?:  <str>,
      "passage_id"?:           <str>,
      "assay_type":            "bulkWGS" | "bulkWGS-lr" | "bulkWES" | "scDNA" |
                               "scRNA" | "scATAC" | "VISIUM" | "VISIUM3prime",
      "sample_type":           "normal" | "tumor",
      "reference_version":     <str>,
      "platform"?:             <str>,
      "cancer_type"?:          <str>,
      "files": {
        "alignment":         <path|url>,
        "alignment_index":   <path|url>,
        "barcodes"?:         <path|url>,
        "fragments"?:        <path|url>,
        "matrix_h5"?:        <path|url>,
        "tissue_positions"?: <path|url>,
        "scalefactors"?:     <path|url>,
        "image_hires"?:      <path|url>,
        "image_lowres"?:     <path|url>,
      }
    },
    ...
  ]
}
```

## Fields

| Key | Required | Description |
|-----|----------|-------------|
| `sample_id` | Yes | Patient ID; `[A-Za-z0-9-]+`. |
| `dataset_id` | Yes | Dataset ID, unique with each patient; `[A-Za-z0-9-]+`. Neither id may hold `_`: both are fields of `_`-joined keys the pipeline splits positionally (`{raw}_{dataset_id}_{assay_type}`, `{sample_id}_{dataset_id}`). |
| `assay_type` | Yes | `bulkWGS` \| `bulkWGS-lr` \| `bulkWES` \| `scDNA` \| `scRNA` \| `scATAC` \| `VISIUM` \| `VISIUM3prime`. `scDNA` runs only under `bulk_genotyping`; see [scDNA](#scdna). |
| `sample_type` | Yes | `normal` \| `tumor`. |
| `reference_version` | Yes | reference version, e.g. `GRCh38`; see [Reference version](#reference-version). |
| `files` | Yes | Input files; see [Files](#files). |
| `rdr_base_dataset_id` | No | matched-normal sample `dataset_id`; see [RDR normalization](reference.md#params_combine_counts). |
| `passage_id` | No | e.g. `p23`. |
| `platform` | No | e.g. `Illumina PCR-free`. |
| `cancer_type` | No | e.g. `PDAC`. |

## Files

| Key | Required for | Example |
|-----|--------------|--------------------|
| `alignment` | every assay type | `outs/*_possorted_bam.bam` |
| `alignment_index` | every assay type | `outs/*_possorted_bam.bam.bai` |
| `barcodes` | `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix/barcodes.tsv.gz` |
| `fragments` | `scATAC` | `outs/atac_fragments.tsv.gz` |
| `matrix_h5` | `scRNA`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix.h5` |
| `tissue_positions` | `VISIUM`, `VISIUM3prime` | `outs/spatial/tissue_positions.csv`; Visium HD: `outs/binned_outputs/square_XXXum/spatial/tissue_positions.parquet` |
| `scalefactors` | `VISIUM`, `VISIUM3prime` | `outs/spatial/scalefactors_json.json` |
| `image_hires` | `VISIUM` | `outs/spatial/tissue_hires_image.png` |
| `image_lowres` | `VISIUM` | `outs/spatial/tissue_lowres_image.png` |

## Features
### scDNA
`scDNA` names a single-cell DNA library (10x CNV, DLP+, ACT). It is accepted **only**
under `workflow_mode: bulk_genotyping`, where it is treated as an ordinary bulk DNA
alignment: every barcode's reads are pooled into one observation, the cell tag is never
read, and the record needs `alignment` + `alignment_index` and nothing else. A `scDNA`
dataset can carry `rdr_base_dataset_id` and be genotyped, phased and binned like any
other bulk dataset, and `bulkWGS` is still preferred as the genotyping source when both
are present.

A per-cell library carries one `@RG` `SM` tag per barcode, which `bcftools mpileup`
would call as one sample each, so both bcftools rules take `--ignore-RG` for a `scDNA`
dataset: one BAM becomes one sample and the barcodes' reads are summed. Because that
flag is per input FILE, a `scDNA` dataset must be genotyped alone; naming it alongside
another dataset in `genotype_dataset_ids` is a parse-time error.

Parsing emits a `WARNING` naming each pooled dataset, so a run that produces no
per-cell output says so up front. `single_cell_genotyping` and `copytyping_preprocess`
reject `scDNA`: their per-cell counting has no source for it.

### Support multiome dataset
A 10x Epi Multiome dataset is two records sharing one (`sample_id`, `dataset_id`), with `assay_type` set to `scRNA` and `scATAC`, respectively.

### Reference version
`reference_version` names the reference genome build the given dataset was aligned to, and matches config's `reference`. The workflow only selects the datasets whose `reference_version` exactly matches the config's [`reference_version`](./reference.md#input-data) under aliasing:

| Canonical | Also accepted (case-insensitive) |
|-----------|----------------------------------|
| `hg19` | `GRCh37`, `b37`, `hs37`, `hs37d5` |
| `hg38` | `GRCh38`, `hs38`, `hs38DH`, `GRCh38.p13`, `GRCh38_no_alt` |
| `chm13v2` | `CHM13v2.0`, `CHM13`, `T2T`, `T2T-CHM13v2`, `T2T-CHM13v2.0` |
| `mm10` | `GRCm38`, `MGSCv38` |

> [!IMPORTANT]
> If CRAM format is used, `reference_version` and config's `reference` must match exactly
> to the input CRAM, otherwise CRAM's MD5 checksum will fail.

### Support remote files
A input file could be a local disk path or an `http(s)` URL to a FTP server, see [`remote_mode`](./reference.md#input-data) for more details. Here is an example for the GIAB HG008 tumor/normal cell line (bulkWGS) served from the NCBI GIAB FTP:
```json
{
  "version": 1,
  "samples": [
    {
      "sample_id": "HG008",
      "dataset_id": "N",
      "assay_type": "bulkWGS",
      "sample_type": "normal",
      "reference_version": "GRCh38",
      "files": {
        "alignment": "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/HG008-N.bam",
        "alignment_index": "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/HG008-N.bam.bai"
      }
    },
    {
      "sample_id": "HG008",
      "dataset_id": "T",
      "rdr_base_dataset_id": "N",
      "assay_type": "bulkWGS",
      "sample_type": "tumor",
      "reference_version": "GRCh38",
      "files": {
        "alignment": "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/HG008-T.bam",
        "alignment_index": "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/HG008-T.bam.bai"
      }
    }
  ]
}
```
> [!IMPORTANT]
> In `stream` mode, `alignment_index` must also be reachable via URL.

## Supporting TSV-formatted sample sheet

User may also use TSV-formatted sample sheet with same schema as JSON file but strongly unrecommended. All files column must have prefix `files.`. Here is an example:

```
sample_id  dataset_id  assay_type  sample_type  reference_version  files.alignment  files.alignment_index
T1         N1          bulkWGS     normal       hg38               /d/normal.bam    /d/normal.bam.bai
T1         D1          bulkWGS     tumor        hg38               /d/tumor.bam     /d/tumor.bam.bai
```
