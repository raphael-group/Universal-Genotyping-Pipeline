# Universal Genotyping Pipeline

[![CI](https://github.com/raphael-group/Universal-Genotyping-Pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/raphael-group/Universal-Genotyping-Pipeline/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.1.0b1-blue.svg)](VERSION)
[![Snakemake](https://img.shields.io/badge/snakemake->=9.0-brightgreen.svg)](https://snakemake.readthedocs.io)

Universal Genotyping Pipeline is a Snakemake preprocessing pipeline for downstream allele-specific CNA inference softwares including:
- [HATCHet](https://github.com/raphael-group/hatchet) - bulk **short-read** WGS/WES, **long-read** PacBio HiFi/Oxford Nanopore,
- [Copy-typing](https://github.com/raphael-group/Copy-typing) - scRNA-seq, scATAC-seq, Epi Multiome,
- [CalicoST](https://github.com/raphael-group/CalicoST) - Visium ST.

---

## Installation

`Universal-Genotyping-Pipeline` requires a 64-bit Linux system and [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/) (version 9 or newer). Create the runner environment from [`environment.yaml`](environment.yaml) and activate it:

```sh
conda env create -f environment.yaml
conda activate genotyping-env
```

Build the conda environments for pipeline dependencies **once** for all future runs:

```sh
snakemake --profile profile/ \
    --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile
```

All pipeline dependencies can be found at [`workflow/envs/`](./workflow/envs/), see [dependencies](./docs/reference.md#dependencies) for details.

---

## Configuration

---
The profile at [`profile/config.yaml`](profile/config.yaml) holds run-wide settings. The main keys to check:

| Key | Usage |
|-----|--------------|
| `cores` | Max CPU cores the whole pipeline may use for job scheduling. |
| `use-conda` | Keep `true`. |
| `conda-prefix` | Where the pipeline dependencies environments live; set to an **absolute path**. |
| `resources: downloads` | How many remote (URL) input files download at once; default `2`. |
| `local-storage-prefix` | Where remote input files are downloaded; set to a large scratch disk. |

> [!IMPORTANT]
> The default `conda-prefix` is the relative path `.snakemake/conda` w.r.t. current working directory. After the environments are built, change it to an absolute path to avoid re-building the environment for every new run.

---

The sample sheet ([template](resources/templates/samples.json)) is a JSON file lists the input datasets & configurations. Copy the template and modify from it according to [schema](docs/sample_sheet.md). Validate the format and check the file existence via:

```sh
python resources/scripts/validate_sample_file.py /path/to/samples.json --check-files
```

---

The config file ([template](resources/templates/config.yaml)) is a YAML file lists the detailed workflow parameters and paths. See [reference.md](docs/reference.md#configuration) for detailed guidelines for each specific workflow: [bulk-genotyping](docs/bulk_genotyping.md), [single-cell-genotyping](docs/single_cell_genotyping.md), [copytyping-preprocess](docs/copytyping_preprocess.md).

---

## Running the pipeline

Run the snakemake pipeline using following commands. `--config` allows user to overwrite config parameters against the config file.

```sh
snakemake --profile profile/ \
    -s workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.json sample_id=<PATIENT_ID>
```


> [!TIP]
> - For the first run, use CMD argument `--dry-run` (`-n`). It lists the jobs Snakemake would run without actual executions, so you can confirm the plan and catch potential configuration errors.
> - If a run failed at intermediate jobs or user changed downstream parameters, use CMD argument 
> `--rerun-incomplete` to resume execution.

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/bulk_genotyping.md](docs/bulk_genotyping.md) | Bulk genotyping workflow |
| [docs/single_cell_genotyping.md](docs/single_cell_genotyping.md) | Single-cell/spatial genotyping workflow |
| [docs/copytyping_preprocess.md](docs/copytyping_preprocess.md) | Copy-typing preprocessing workflow |
| [docs/sample_sheet.md](docs/sample_sheet.md) | Sample sheet specification |
| [docs/reference.md](docs/reference.md) | Reference manual |
| [resources/README.md](resources/README.md) | External data resources |
| [CHANGELOG.md](CHANGELOG.md) | Release notes and version history |
