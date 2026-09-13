#!/usr/bin/env python3
"""Parse and validate the workflow configfile and sample sheet into the rules' globals.

Last update: 2026-08-07

Inputs:
- configfile: the Snakemake config
- sample_file: JSON records, or the same schema flattened to TSV
Outputs:
- parse_workflow(config): dict of globals, unpacked by the Snakefile
References:
- docs/reference.md: config keys and output layout
- docs/sample_sheet.md: the sample-file schema
"""

import csv
import json
import os
import re

from const import (
    ALLOWED_ASSAY_TYPES,
    GT_ASSAY_ORD,
    ASSAY_TYPE2MODALITY,
    BULK_ASSAYS,
    BULK_LR_ASSAYS,
    BULK_TARGETS,
    COPYTYPING_TARGETS,
    FILES_COLUMN_PREFIX,
    LONGREAD_PHASER,
    MULTIOME_ASSAYS,
    NONBULK_ASSAYS,
    PSEUDOBULK_ASSAYS,
    PANEL_PHASER,
    RECORD_ID_KEYS,
    RECORD_ID_PATTERN,
    REMOTE_MODES,
    REPLISEQ_REFVERS,
    RDR_NORMALIZATIONS,
    REFVERS,
    REQUIRED_FILES,
    REQUIRED_RECORD_KEYS,
    SAMPLE_FILE_EXTS,
    SAMPLE_TYPES,
    SCALAR_RECORD_KEYS,
    SINGLE_CELL_TARGETS,
    SNP_PANEL_EXTS,
    SPECIES,
    WORKFLOW_MODES,
    canonical_refver,
    get_genetic_map_path,
    is_known_refver,
    get_phasing_panel_path,
)
from io_utils import read_chrom_sizes
from utils import check_local_path, is_url, logging_snakemake, strip_chr_prefix

_RECORD_ID_RE = re.compile(RECORD_ID_PATTERN)


def read_sample_sheet(path):
    """Read the sample file, spec: docs/sample_sheet.md.

    Args:
        path: Path to the sample file.

    Returns:
        List of record dicts, with scalar fields coerced to str.

    Raises:
        AssertionError: The file has the wrong shape, or a record is malformed,
            missing a required key, or carries an id outside RECORD_ID_PATTERN.
    """
    ext = os.path.splitext(path)[1].lower()
    assert ext in SAMPLE_FILE_EXTS, (
        f"{path}: extension must be one of {SAMPLE_FILE_EXTS}, got {ext!r}"
    )
    records = []
    if ext == ".json":
        with open(path) as fh:
            doc = json.load(fh)
        assert isinstance(doc, dict), f"{path}: not a JSON object"
        samples = doc.get("samples")
        assert isinstance(samples, list), f"{path}: `samples` is not a list"
        for idx, rec in enumerate(samples):
            assert isinstance(rec, dict), f"{path}: record {idx} is not an object"
            norm = dict(rec)
            for key in SCALAR_RECORD_KEYS:
                if key in norm and norm[key] is not None:
                    norm[key] = str(norm[key])
            files = norm.get("files")
            assert isinstance(norm["files"], dict), (
                f"{path}: record {idx}, `files` is not a dict"
            )
            if isinstance(files, dict):
                norm["files"] = {k: str(v) for k, v in files.items() if v is not None}
            records.append(norm)
    else:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for row in rows:
            rec, files = {}, {}
            for col, val in row.items():
                val = (val or "").strip()
                if not val or col is None:
                    continue
                if col.startswith(FILES_COLUMN_PREFIX):
                    files[col[len(FILES_COLUMN_PREFIX) :]] = val
                else:
                    rec[col] = val
            rec["files"] = files
            records.append(rec)

    assert records, f"{path}: no records"
    for idx, rec in enumerate(records):
        missing = [k for k in REQUIRED_RECORD_KEYS if rec.get(k) in (None, "")]
        assert not missing, f"{idx}th dataset: missing required key(s) {missing}"
        for key in RECORD_ID_KEYS:
            assert _RECORD_ID_RE.fullmatch(rec[key]), (
                f"{idx}th dataset: {key} {rec[key]!r} must match {RECORD_ID_PATTERN}"
            )
    return records


def parse_records(
    records: list,
    sample_id: str,
    reference_version: str,
    config_assay_types: list,
    select_dataset_ids: list = None,
):
    """Select and parse sample records based on sample_id, assay types, and reference build.

    Filters run in the order sample_id -> assay_type -> select_dataset_ids ->
    reference_version, so a requested dataset_id whose records all sit on another
    build is reported as a build mismatch rather than as a missing dataset. A
    dataset_id carrying rows for several builds keeps only the matching rows. With
    no selection, records of another build are dropped silently.

    A multiome pair shares one dataset_id, so selecting that id keeps both of its
    records.

    Args:
        records: Records from read_sample_sheet.
        sample_id: The sample_id being processed.
        reference_version: canonical reference version.
        config_assay_types: Assay types enabled for this run.
        select_dataset_ids: Restrict to these dataset_ids; empty/None keeps all.

    Returns:
        The selected datasets.

    Raises:
        AssertionError: The selection is empty, a requested dataset_id does not exist
            or has no record on reference_version, or a record violates the spec or a
            replicate rule.
    """
    select = set(select_dataset_ids or [])
    parsed_records = []
    available_ids = set()
    id2refvers = {}
    for rec in records:
        if rec["sample_id"] != sample_id:
            continue
        assay_type = rec["assay_type"]
        if assay_type not in config_assay_types:
            continue
        dataset_id = rec["dataset_id"]
        available_ids.add(dataset_id)
        if select and dataset_id not in select:
            continue
        rec_refver = rec["reference_version"]
        id2refvers.setdefault(dataset_id, set()).add(rec_refver)
        if canonical_refver(rec_refver) != reference_version:
            continue
        sample_type = rec["sample_type"]
        assert rec["sample_type"] in SAMPLE_TYPES, (
            f"{dataset_id}: sample_type must be one of {SAMPLE_TYPES}, "
            f"got {sample_type!r}"
        )
        files = rec["files"]
        readable = REQUIRED_FILES[assay_type]
        files = {k: v for k, v in files.items() if k in readable}
        for key in sorted(REQUIRED_FILES[assay_type]):
            assert files.get(key), (
                f"{dataset_id}: files.{key} is required for {assay_type}"
            )
        for key, path in files.items():
            check_local_path(path, f"{dataset_id}: files.{key}")
        rec["files"] = files
        rec["modality"] = ASSAY_TYPE2MODALITY[assay_type]
        parsed_records.append(rec)

    if select:
        missing = sorted(select - available_ids)
        assert not missing, (
            f"dataset_ids not found for sample_id={sample_id!r}: {missing}; "
            f"available: {sorted(available_ids)}"
        )
        wrong_build = sorted(select - {rec["dataset_id"] for rec in parsed_records})
        assert not wrong_build, (
            "dataset_ids have no record on "
            f"reference_version={reference_version!r}: "
            + "; ".join(f"{d} is on {sorted(id2refvers[d])}" for d in wrong_build)
        )
    assert parsed_records, "no datasets exist after selection"

    parsed_records.sort(
        key=lambda r: (
            GT_ASSAY_ORD.get(r["assay_type"], len(GT_ASSAY_ORD)),
            r["sample_type"] != "normal",
            r["dataset_id"],
        )
    )

    dataset2assays = {}
    for rec in parsed_records:
        dataset2assays.setdefault(rec["dataset_id"], []).append(rec["assay_type"])
    for dataset_id, assays in dataset2assays.items():
        if len(assays) == 1:
            continue
        assert len(assays) == 2 and set(assays) == MULTIOME_ASSAYS, (
            f"dataset_id {dataset_id!r}, assays {assays} may not share one "
            f"({' + '.join(sorted(MULTIOME_ASSAYS))} only)"
        )

    dataset_ids = {rec["dataset_id"] for rec in parsed_records}
    for rec in parsed_records:
        rdr_base_id = rec.get("rdr_base_dataset_id")
        if rdr_base_id:
            dataset_id = rec["dataset_id"]
            assert rec["sample_type"] == "tumor", (
                f"{dataset_id}: rdr_base_dataset_id set on a non-tumor record"
            )
            assert rdr_base_id != dataset_id, (
                f"{dataset_id}: rdr_base_dataset_id is itself"
            )
            assert rdr_base_id in dataset_ids, (
                f"{dataset_id}: rdr_base_dataset_id {rdr_base_id!r} is not among the "
                "selected datasets; add it to dataset_ids or drop the key"
            )
    return parsed_records


def parse_workflow(config):
    """Parse and validate the config + sample file into the workflow's globals.

    Args:
        config: The Snakemake config dict.

    Returns:
        Dict: parsed configurations

    Raises:
        AssertionError: Invalid configuration.
    """

    def require_per_chrom(get_path, label):
        """Every configured chromosome has its file."""
        for chrname in nochr_chromosomes:
            check_local_path(get_path(chrname), f"{label} (chromosome {chrname})")

    # === workflow mode ===
    workflow_mode = config["workflow_mode"]
    assert workflow_mode in WORKFLOW_MODES, (
        f"workflow_mode must be one of {list(WORKFLOW_MODES)}, got {workflow_mode!r}"
    )

    # === output directories ===
    snp_dir = config["snp_dir"]
    phase_dir = config["phase_dir"]
    pileup_dir = config["pileup_dir"]
    allele_dir = config["allele_dir"]
    bb_dir = config["bb_dir"]
    qc_dir = config["qc_dir"]
    log_dir = config["log_dir"]
    aux_dir = config["aux_dir"]
    bench_dir = config["bench_dir"]

    # === assay_types in configfile ===
    config_assay_types = config["assay_types"]
    invalid_assay_types = [
        a for a in config_assay_types if a not in ALLOWED_ASSAY_TYPES
    ]
    assert not invalid_assay_types, (
        f"assay_types {invalid_assay_types} not in {sorted(ALLOWED_ASSAY_TYPES)}"
    )
    allowed = BULK_ASSAYS if workflow_mode == "bulk_genotyping" else NONBULK_ASSAYS
    config_assay_types = [a for a in config_assay_types if a in allowed]
    assert len(config_assay_types) > 0, (
        f"assay_types, none valid for workflow_mode={workflow_mode}"
    )

    # === sample_file ===
    sample_id = config["sample_id"]
    sample_file = config["sample_file"]
    records = read_sample_sheet(sample_file)

    # === reference version ===
    raw_refver = config["reference_version"]
    assert raw_refver, f"reference_version is required, one of {REFVERS} (or an alias)"
    reference_version = canonical_refver(raw_refver)
    if not is_known_refver(raw_refver):
        logging_snakemake(
            f"WARNING: reference_version={raw_refver!r} is not natively supported "
            f"({REFVERS})."
        )
    reference = config["reference"]
    assert reference, "reference is required (genome FASTA)"
    check_local_path(reference, "reference")

    select_dataset_ids = config["dataset_ids"]
    assert len(set(select_dataset_ids)) == len(select_dataset_ids), (
        f"dataset_ids has duplicate dataset_ids: {select_dataset_ids}"
    )
    records = parse_records(
        records, sample_id, reference_version, config_assay_types, select_dataset_ids
    )
    dataset_ids = {r["dataset_id"]: r for r in records}
    assay_types = list(dict.fromkeys(rec["assay_type"] for rec in records))

    if select_dataset_ids:
        logging_snakemake(f"dataset_ids restricts the run to: {select_dataset_ids}")
    logging_snakemake("The following datasets will be processed:")
    for rec in records:
        logging_snakemake(
            f"dataset: {rec['dataset_id']}\t{rec['reference_version']}\t"
            f"{rec['assay_type']}\t{rec['sample_type']}"
        )

    pooled = [r for r in records if r["assay_type"] in PSEUDOBULK_ASSAYS]
    if pooled:
        logging_snakemake(
            "WARNING: single-cell dataset(s) "
            + ", ".join(f"{r['dataset_id']} ({r['assay_type']})" for r in pooled)
            + " are treated as pseudobulks ignoring cell barcodes."
        )

    # === chromosomes: the run's list, spelled both ways ===
    nochr_chromosomes = [strip_chr_prefix(c) for c in config["chromosomes"]]
    assert nochr_chromosomes, "chromosomes is empty"

    genome_size = config["genome_size"]
    assert genome_size, "genome_size is required (two-column chrom<TAB>size file)"
    ref_chroms_nochr = {}
    for name in read_chrom_sizes(genome_size):
        ref_chroms_nochr.setdefault(strip_chr_prefix(name), name)

    absent_chroms = [c for c in nochr_chromosomes if c not in ref_chroms_nochr]
    assert not absent_chroms, (
        f"chromosomes {absent_chroms} are not found in {genome_size}"
    )
    chr_chromosomes = [f"chr{c}" for c in nochr_chromosomes]
    logging_snakemake(f"chromosomes: {chr_chromosomes}")

    input_nochr = True
    if ref_chroms_nochr[nochr_chromosomes[0]].lower().startswith("chr"):
        for chrname in ref_chroms_nochr.values():
            assert chrname.lower().startswith("chr"), (
                f"chrom.sizes having mixing chr and nochr-prefix chroms: {chrname}"
            )
        input_nochr = False
    logging_snakemake(f"chrom.sizes chr-prefix={not input_nochr}")

    # === remote input mode: whole-file storage() download vs direct URL streaming ===
    remote_mode = config["remote_mode"]
    assert remote_mode in REMOTE_MODES, (
        f"remote_mode must be one of {REMOTE_MODES}, got {remote_mode!r}"
    )
    if remote_mode == "stream":
        assert workflow_mode == "bulk_genotyping", (
            "remote_mode='stream' is only supported for bulk_genotyping"
        )

    # === species ===
    species = config["species"]
    assert species, f"species is required, one of {list(SPECIES)}"
    if species not in SPECIES:
        logging_snakemake(
            f"WARNING: species={species!r} is not natively supported ({list(SPECIES)})."
        )

    # === gtf_file ===
    gtf_file = config["gtf_file"]
    assert gtf_file, "gtf_file is required"
    check_local_path(gtf_file, "gtf_file")

    # === segment BED: the arms cut at the SV extremities, blacklist subtracted ===
    region_bed = config["region_bed"]
    assert region_bed, "region_bed is required (chromosome arms)"
    check_local_path(region_bed, "region_bed")
    extremity_tsv = config["extremity_tsv"]
    if extremity_tsv:
        check_local_path(extremity_tsv, "extremity_tsv")
        logging_snakemake(f"cut region_bed arms at SV extremities: {extremity_tsv}")
    else:
        extremity_tsv = []
        logging_snakemake("extremity_tsv unset, one segment per region_bed arm")
    segment_bed = aux_dir + "/segment.bed"

    # === blacklist files ===
    blacklist_bed = config["blacklist_bed"]
    if blacklist_bed:
        check_local_path(blacklist_bed, "blacklist_bed")
    else:
        blacklist_bed = []

    gene_blacklist_file = config["gene_blacklist_file"]
    if gene_blacklist_file:
        check_local_path(gene_blacklist_file, "gene_blacklist_file")
    else:
        gene_blacklist_file = []

    mappability_bed = config["mappability_bed"]
    if mappability_bed:
        check_local_path(mappability_bed, "mappability_bed")
    else:
        mappability_bed = []

    # === window BED ===
    window_size = int(config["params_build_windows"]["window_size"])
    window_bed = config["window_bed"]
    if window_bed and extremity_tsv:
        logging_snakemake(
            f"NOTE: extremity_tsv is set, so window_bed={window_bed} is ignored; "
            "windows are re-tiled from the cut segments"
        )
        window_bed = None
    if window_bed:
        check_local_path(window_bed, "window_bed")
        logging_snakemake(f"use pre-built window BED: {window_bed}")
        build_windows = False
    else:
        build_windows = True
        logging_snakemake(
            f"build window BED from {segment_bed}, window_size={window_size}"
        )
        window_bed = aux_dir + "/windows.bed.gz"

    # === capture targets ===
    target_bed = config.get("target_bed")
    window_target = []
    if target_bed:
        check_local_path(target_bed, "target_bed")
        window_target = aux_dir + "/window.target.npz"
        if "bulkWES" not in assay_types:
            logging_snakemake(
                f"NOTE: target_bed={target_bed} is set but no bulkWES dataset is in the "
                "run; the target sites are annotated and then unused"
            )
        else:
            logging_snakemake(
                f"split bulkWES windows on/off capture target: {target_bed}"
            )
    else:
        target_bed = []

    # === pre-built files ===
    het_snp_vcf = config["het_snp_vcf"]
    het_snp_vcf_phased = bool(config["het_snp_vcf_phased"])
    if het_snp_vcf is not None:
        check_local_path(het_snp_vcf, "het_snp_vcf")

    bb_file = config["bb_file"]
    if bb_file:
        check_local_path(bb_file, "bb_file")
    else:
        bb_file = []

    # === copytyping_preprocess requirements ===
    if workflow_mode == "copytyping_preprocess":
        assert het_snp_vcf is not None, (
            "het_snp_vcf is required for copytyping_preprocess"
        )
        assert het_snp_vcf_phased, (
            "het_snp_vcf must be phased for copytyping_preprocess"
        )
        assert bb_file, "bb_file is required for copytyping_preprocess"

    # === genotyping check ===
    run_genotyping = True
    if het_snp_vcf is not None:
        logging_snakemake(f"run_genotyping is skipped, use input {het_snp_vcf}")
        run_genotyping = False

    snp_panel = config["snp_panel"]

    # the strategy that actually runs; None passes the caller's germline calls through
    tumor_genotyping_mode = None
    genotype_files = None
    genotype_ignore_rg = False
    # NB: cellsnp-lite -R already genotypes at the panel's alleles
    panel_allele_only = bool(config["panel_allele_only"]) and (
        run_genotyping and workflow_mode == "bulk_genotyping"
    )
    if run_genotyping:
        assert snp_panel, f"snp_panel is required for {workflow_mode}"
        check_local_path(snp_panel, "snp_panel")
        assert snp_panel.endswith(SNP_PANEL_EXTS), (
            f"snp_panel: extension must be one of {SNP_PANEL_EXTS} "
            f"(bcftools --targets-file cannot read BCF), got {snp_panel!r}"
        )
        if not is_url(snp_panel):
            assert os.path.exists(snp_panel + ".tbi") or os.path.exists(
                snp_panel + ".csi"
            ), f"snp_panel is not indexed, run: tabix -p vcf {snp_panel}"
        genotype_dataset_ids = config["genotype_dataset_ids"]
        genotype_records = []
        if genotype_dataset_ids:
            assert len(set(genotype_dataset_ids)) == len(genotype_dataset_ids), (
                f"genotype_dataset_ids has duplicate dataset_ids: {genotype_dataset_ids}"
            )
            missing = [d for d in genotype_dataset_ids if d not in dataset_ids]
            assert not missing, (
                f"genotype_dataset_ids not found in the records: {missing}"
            )
            genotype_records = [
                rec for rec in records if rec["dataset_id"] in set(genotype_dataset_ids)
            ]
        if len(genotype_records) == 0:
            # the first bulkWGS normal, falling through to a tumor when none exists
            genotype_records = sorted(
                records,
                key=lambda r: (
                    r["sample_type"] != "normal",
                    GT_ASSAY_ORD.get(r["assay_type"], len(GT_ASSAY_ORD)),
                ),
            )
            assert len(genotype_records) > 0, (
                "genotype_dataset_ids, no dataset to genotype"
            )
            genotype_records = genotype_records[:1]
        genotype_dataset_ids = [rec["dataset_id"] for rec in genotype_records]
        logging_snakemake(f"genotype_dataset_ids: {genotype_dataset_ids}")
        if workflow_mode == "bulk_genotyping":
            # only bulk can fall back on the caller's germline GT; cellsnp-lite emits none
            genotype_tumor = any(
                rec["sample_type"] == "tumor" for rec in genotype_records
            )
            tumor_genotyping_mode = "vaf_cutoff" if genotype_tumor else None
            logging_snakemake(
                f"tumor_genotyping_mode={tumor_genotyping_mode or 'none'} "
                f"(genotyped tumor={genotype_tumor})"
            )
            if genotype_tumor and not panel_allele_only:
                logging_snakemake(
                    "genotyping a tumor: forcing panel_allele_only, so ALT comes from "
                    "snp_panel and a somatic allele cannot become the called ALT"
                )
            panel_allele_only = panel_allele_only or genotype_tumor
        else:
            tumor_genotyping_mode = "vaf_cutoff"
            logging_snakemake(
                f"tumor_genotyping_mode={tumor_genotyping_mode} "
                "(single-cell genotypes from the pseudobulk counts)"
            )
        genotype_files = [rec["files"] for rec in genotype_records]
        # a per-cell library carries one @RG SM per barcode, which mpileup would call as
        # one sample each; --ignore-RG pools them, but then one sample per input FILE
        genotype_ignore_rg = any(
            rec["assay_type"] in PSEUDOBULK_ASSAYS for rec in genotype_records
        )
        if genotype_ignore_rg:
            assert len(genotype_records) == 1, (
                f"genotype_dataset_ids pools {genotype_dataset_ids}, including a "
                f"per-cell library {sorted(PSEUDOBULK_ASSAYS)}; bcftools --ignore-RG "
                "yields one sample per input file, so genotype exactly one dataset"
            )
            logging_snakemake(
                f"genotyping {genotype_dataset_ids[0]} with bcftools --ignore-RG: its "
                "per-barcode @RG SM tags collapse into one pooled sample"
            )

    # === phasing check ===
    run_phasing = True
    phased_snp_vcf = phase_dir + "/phased_het_snps.vcf.gz"
    if het_snp_vcf is not None and het_snp_vcf_phased:
        phased_snp_vcf = het_snp_vcf
        logging_snakemake(f"run_phasing is skipped, use input {het_snp_vcf}")
        run_phasing = False

    phaser = config["phaser"] if run_phasing else None
    phase_files = None
    gmap_file = []
    get_genetic_map = None
    get_phasing_panel = None
    if run_phasing:
        assert phaser in PANEL_PHASER | LONGREAD_PHASER, f"unknown phaser: {phaser}"
        if phaser in PANEL_PHASER:
            gmap_path = config["gmap_path"]
            assert gmap_path, f"gmap_path is required for {phaser}"
            get_genetic_map = get_genetic_map_path(gmap_path)
            require_per_chrom(get_genetic_map, "gmap file")

            phasing_panel = config["phasing_panel"]
            assert phasing_panel and os.path.isdir(phasing_panel), (
                f"phasing_panel is not a directory: {phasing_panel}"
            )
            get_phasing_panel = get_phasing_panel_path(phasing_panel)
            require_per_chrom(get_phasing_panel, "panel file")
            logging_snakemake(
                f"phaser={phaser}, gmap={gmap_path}, panel={phasing_panel}"
            )
            gmap_file = phase_dir + "/genetic_map.tsv.gz"

        if phaser in LONGREAD_PHASER:
            phase_dataset_ids = config["phase_dataset_ids"]
            phase_records = []
            if phase_dataset_ids:
                assert len(set(phase_dataset_ids)) == len(phase_dataset_ids), (
                    f"phase_dataset_ids has duplicate dataset_ids: {phase_dataset_ids}"
                )
                missing = [d for d in phase_dataset_ids if d not in dataset_ids]
                assert not missing, (
                    f"phase_dataset_ids not found in the records: {missing}"
                )
                phase_records = [
                    rec
                    for rec in records
                    if rec["dataset_id"] in set(phase_dataset_ids)
                ]
            if len(phase_records) == 0:
                # co-phases all long-read normals, else all long-read tumors
                lr_records = [r for r in records if r["assay_type"] in BULK_LR_ASSAYS]
                normals = [r for r in lr_records if r["sample_type"] == "normal"]
                phase_records = normals or lr_records
            assert phase_records, (
                f"{phaser} requires a long-read assay {sorted(BULK_LR_ASSAYS)}"
            )
            short_read = [
                rec["dataset_id"]
                for rec in phase_records
                if rec["assay_type"] not in BULK_LR_ASSAYS
            ]
            assert not short_read, f"{phaser} needs long reads, got: {short_read}"
            phase_dataset_ids = [rec["dataset_id"] for rec in phase_records]
            logging_snakemake(f"phase_dataset_ids: {phase_dataset_ids}")
            non_normal = [
                rec["dataset_id"]
                for rec in phase_records
                if rec["sample_type"] != "normal"
            ]
            if non_normal:
                logging_snakemake(
                    f"WARN: phase_dataset_ids includes non-normal dataset(s) "
                    f"{non_normal}; phasing from tumor reads"
                )
            phase_files = [rec["files"] for rec in phase_records]

    # === RDR normalization mode ===
    assert "params_mosdepth" not in config, (
        "params_mosdepth was folded into params_count_reads; move read_quality there and "
        "rename extra_params to mosdepth_extra_params"
    )
    do_repliseq = False
    if workflow_mode == "bulk_genotyping":
        rdr_normalization = config["params_combine_counts"]["rdr_normalization"]
        assert rdr_normalization in RDR_NORMALIZATIONS, (
            f"rdr_normalization must be one of {list(RDR_NORMALIZATIONS)}, "
            f"got {rdr_normalization!r}"
        )
        do_repliseq = (reference_version in REPLISEQ_REFVERS) and bool(
            config["params_count_reads"]["rt_correct"]
        )
        for rec in records:
            if rec["sample_type"] != "tumor":
                continue
            dataset_id = rec["dataset_id"]
            base = rec.get("rdr_base_dataset_id")
            if rdr_normalization == "median" and base:
                rec.pop("rdr_base_dataset_id")
                logging_snakemake(
                    f"NOTE: {dataset_id}: rdr_normalization='median' -> "
                    f"ignoring rdr_base_dataset_id={base}"
                )
            elif rdr_normalization == "normal":
                assert base, (
                    f"{dataset_id}: rdr_normalization='normal' requires rdr_base_dataset_id"
                )
            elif rdr_normalization == "auto" and not base:
                logging_snakemake(
                    f"NOTE: {dataset_id}: no rdr_base_dataset_id, RDR median-normalized"
                )

    # === min_snp_reads sweep (one MSR{msr}/ subdir per value) ===
    assert "max_blocksize" not in config["params_combine_counts"], (
        "params_combine_counts.max_blocksize was removed: a span cap can only fire by "
        "cutting a bb that has not met its read thresholds. Use min_total_reads, the "
        "read-start floor every bb must clear"
    )
    msr_raw = config["params_combine_counts"]["min_snp_reads"]
    msr_list = [int(m) for m in (msr_raw if isinstance(msr_raw, list) else [msr_raw])]

    # === Snakemake rule's lookup tables ===
    modalities = list(dict.fromkeys(rec["modality"] for rec in records))
    modality2files = {}
    for rec in records:
        modality2files.setdefault(rec["modality"], []).append(rec["files"])
    by_assay = {}
    for rec in records:
        by_assay.setdefault(rec["assay_type"], []).append(rec)
    assay2dataset_ids, assay2sample_types, assay2base_dataset_ids = {}, {}, {}
    for assay_type in ALLOWED_ASSAY_TYPES:
        rows = by_assay.get(assay_type, [])
        assay2dataset_ids[assay_type] = [rec["dataset_id"] for rec in rows]
        assay2sample_types[assay_type] = [rec["sample_type"] for rec in rows]
        assay2base_dataset_ids[assay_type] = [
            rec.get("rdr_base_dataset_id", "") for rec in rows
        ]
    get_data = {(rec["assay_type"], rec["dataset_id"]): rec["files"] for rec in records}

    # === final targets  ===
    if workflow_mode == "bulk_genotyping":
        final_targets = [
            f"{bb_dir}/MSR{m}/bulk/{f}" for m in msr_list for f in BULK_TARGETS
        ]
    elif workflow_mode == "single_cell_genotyping":
        final_targets = [
            f"{bb_dir}/MSR{m}/{at}/{f}"
            for at in assay_types
            for m in msr_list
            for f in SINGLE_CELL_TARGETS
        ]
    else:  # copytyping_preprocess
        final_targets = [
            f"{bb_dir}/{at}/{f}" for at in assay_types for f in COPYTYPING_TARGETS
        ]

    # === assemble the globals the Snakefile unpacks ===
    return {
        "workflow_mode": workflow_mode,
        "sample_id": sample_id,
        "remote_mode": remote_mode,
        "reference_version": reference_version,
        "species": species,
        "nochr_chromosomes": nochr_chromosomes,
        "chr_chromosomes": chr_chromosomes,
        "input_nochr": input_nochr,
        "assay_types": assay_types,
        "modalities": modalities,
        "msr_list": msr_list,
        "phaser": phaser,
        "run_genotyping": run_genotyping,
        "run_phasing": run_phasing,
        "het_snp_vcf": het_snp_vcf,
        "phased_snp_vcf": phased_snp_vcf,
        "gmap_file": gmap_file,
        "final_targets": final_targets,
        "get_data": get_data,
        "modality2files": modality2files,
        "assay2dataset_ids": assay2dataset_ids,
        "assay2sample_types": assay2sample_types,
        "assay2base_dataset_ids": assay2base_dataset_ids,
        "genotype_files": genotype_files,
        "genotype_ignore_rg": genotype_ignore_rg,
        "phase_files": phase_files,
        "get_genetic_map": get_genetic_map,
        "get_phasing_panel": get_phasing_panel,
        "snp_dir": snp_dir,
        "phase_dir": phase_dir,
        "pileup_dir": pileup_dir,
        "allele_dir": allele_dir,
        "bb_dir": bb_dir,
        "qc_dir": qc_dir,
        "log_dir": log_dir,
        "aux_dir": aux_dir,
        "bench_dir": bench_dir,
        "extremity_tsv": extremity_tsv,
        "snp_panel": snp_panel,
        "tumor_genotyping_mode": tumor_genotyping_mode,
        "panel_allele_only": panel_allele_only,
        "reference": reference,
        "genome_size": genome_size,
        "gtf_file": gtf_file,
        "bb_file": bb_file,
        "mappability_bed": mappability_bed,
        "region_bed": region_bed,
        "blacklist_bed": blacklist_bed,
        "gene_blacklist_file": gene_blacklist_file,
        "segment_bed": segment_bed,
        "build_windows": build_windows,
        "do_repliseq": do_repliseq,
        "window_bed": window_bed,
        "window_size": window_size,
        "target_bed": target_bed,
        "window_target": window_target,
    }
