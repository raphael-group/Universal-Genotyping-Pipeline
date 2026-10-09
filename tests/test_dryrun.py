#!/usr/bin/env python3
"""Dry-run the DAG for every workflow mode, from both sample-file formats.

Last update: 2026-08-06

Covers:
- modes: each of the three builds a non-empty DAG
- formats: a TSV sheet plans exactly the same jobs as JSON
- validation: bad sample sheets and configs fail at DAG build
- wiring: window build, repliseq, streaming, and the short-circuit VCFs
Notes:
- cost: no rule runs and no real data is needed
"""

import json
import os

import pytest

from conftest import dryrun, job_counts

# (id, sample-file key, sample_id, workflow_mode, assay_types)
BULK_CASES = [
    ("bulk", "bulk", "T1", "bulk_genotyping", ["bulkWGS"]),
]
SC_CASES = [
    ("multiome", "sc", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ("visium", "sc", "V1", "single_cell_genotyping", ["VISIUM"]),
    ("copytyping", "sc", "S1", "copytyping_preprocess", ["scATAC"]),
]
ALL_CASES = BULK_CASES + SC_CASES


@pytest.mark.parametrize("fmt", ["json", "tsv"])
@pytest.mark.parametrize(
    "case_id,sheet,sample_id,mode,assays",
    ALL_CASES,
    ids=[c[0] for c in ALL_CASES],
)
def test_dag_builds(workspace, fmt, case_id, sheet, sample_id, mode, assays):
    """Every mode builds a non-empty DAG from either sample-file format."""
    proc = dryrun(workspace, workspace[f"{sheet}_{fmt}"], sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts, f"no jobs planned:\n{proc.stdout[-2000:]}"
    assert "all" in counts


@pytest.mark.parametrize(
    "case_id,sheet,sample_id,mode,assays",
    ALL_CASES,
    ids=[c[0] for c in ALL_CASES],
)
def test_json_and_tsv_agree(workspace, case_id, sheet, sample_id, mode, assays):
    """A TSV sheet plans exactly the same jobs as the equivalent JSON."""
    a = dryrun(workspace, workspace[f"{sheet}_json"], sample_id, mode, assays)
    b = dryrun(workspace, workspace[f"{sheet}_tsv"], sample_id, mode, assays)
    assert a.returncode == 0 and b.returncode == 0
    assert job_counts(a.stdout) == job_counts(b.stdout)


def test_bulk_rules(workspace):
    """Bulk genotyping plans the genotype/phase/pileup/depth chain."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    counts = job_counts(proc.stdout)
    assert {"genotype_snps_bulk", "phase_snps_eagle", "combine_counts"} <= set(counts)
    # one pileup chunk per (dataset, chromosome) + one merge per dataset; chromosomes=[22]
    assert counts["pileup_snps_bulk_bcftools_chrom"] == 2
    assert counts["merge_pileup_counts"] == 2
    assert counts["run_mosdepth"] == 2
    # read starts fan out per (dataset, chromosome), sharing one per-chrom window BED
    assert counts["count_read_starts_chrom"] == 2
    assert counts["merge_read_starts"] == 2
    assert counts["window_bed_to_3bed_chrom"] == 1
    assert "pileup/bulkWGS/D1.rdcount.bed.gz" in proc.stdout
    # combine_counts aggregates them onto every level's rows
    assert "unit/bulk/window.rdcount.npz" in proc.stdout
    assert "multi_snp/bulk/bb.rdcount.npz" in proc.stdout
    assert "bulk/bb.rdcount.npz" in proc.stdout
    # the caller writes snps/raw/, post-processing writes the file phasing reads
    assert counts["post_genotype_snps_bulk"] == 1
    assert "raw/chr22.vcf.gz" in proc.stdout
    # genotyping targets the panel VCF and scopes the chromosome itself
    assert "--targets-file" in proc.stdout and "snp_panel.vcf.gz" in proc.stdout
    assert "--regions chr22" in proc.stdout


def test_bulk_stream_mode(workspace):
    """remote_mode=stream: mosdepth is per-chrom + merged, bcftools keeps its regions."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=("remote_mode=stream",),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # depth becomes per-chrom mosdepth + a merge; the whole-file rule is gone
    assert "run_mosdepth_chrom" in counts and "merge_mosdepth" in counts
    assert "run_mosdepth" not in counts
    # read starts are already per-chrom, so stream mode reuses the same two rules
    assert counts["count_read_starts_chrom"] == 2
    assert counts["merge_read_starts"] == 2
    # genotype/pileup restrict to the config chroms via index jumps
    assert "--regions chr22" in proc.stdout


def test_read_filters_are_shared_with_mosdepth(workspace):
    """One read_quality/exclude_flags pair reaches both mosdepth and samtools view."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=('params_count_reads={"read_quality": 30, "exclude_flags": 3844}',),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "-Q 30" in proc.stdout and "-q 30" in proc.stdout
    assert proc.stdout.count("-F 3844") >= 2


def test_snp_panel_must_be_vcf_gz(workspace):
    """A BCF panel is refused: bcftools --targets-file matches nothing on BCF."""
    ref = workspace["ref"]
    panel = os.path.join(ref, "snp_panel.bcf")
    open(panel, "wb").close()
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=(f"snp_panel={panel}",),
    )
    assert proc.returncode != 0
    assert "extension must be one of" in (proc.stdout + proc.stderr)


def test_snp_panel_must_be_indexed(workspace):
    """A panel with no .tbi/.csi beside it is refused at parse time."""
    ref = workspace["ref"]
    panel = os.path.join(ref, "unindexed.vcf.gz")
    open(panel, "wb").close()
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=(f"snp_panel={panel}",),
    )
    assert proc.returncode != 0
    assert "not indexed" in (proc.stdout + proc.stderr)


def test_stream_rejected_for_single_cell(workspace):
    """remote_mode=stream errors for non-bulk modes (cellsnp-lite cannot read URLs)."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
        extra=("remote_mode=stream",),
    )
    assert proc.returncode != 0
    assert "only supported for bulk_genotyping" in (proc.stdout + proc.stderr)


@pytest.mark.parametrize(
    "sheet,sample_id,mode,assays",
    [
        ("bulk_json", "T1", "bulk_genotyping", ["bulkWGS"]),
        ("sc_json", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ],
    ids=["bulk", "single_cell"],
)
def test_windows_are_built_from_the_segments(workspace, sheet, sample_id, mode, assays):
    """With no window_bed, every mode builds the segment BED and tiles it."""
    ref = workspace["ref"]
    proc = dryrun(workspace, workspace[sheet], sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_segment_bed" in counts
    assert "build_window_bed" in counts
    # the arms feed build_segment_bed, which feeds the tiling
    assert f"{ref}/region.bed" in proc.stdout
    assert "/windows.bed.gz" in proc.stdout


@pytest.mark.parametrize(
    "mode,assays",
    [
        ("bulk_genotyping", ["bulkWGS"]),
        ("single_cell_genotyping", ["scRNA", "scATAC"]),
        ("copytyping_preprocess", ["scATAC"]),
    ],
)
def test_segments_default_to_region_bed_arms(workspace, mode, assays):
    """An unset extremity_tsv leaves the arms uncut: one segment per arm."""
    sheet = (
        workspace["bulk_json"] if mode == "bulk_genotyping" else workspace["sc_json"]
    )
    sample_id = "T1" if mode == "bulk_genotyping" else "S1"
    proc = dryrun(workspace, sheet, sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "extremity_tsv unset, one segment per region_bed arm" in proc.stdout
    assert job_counts(proc.stdout)["build_segment_bed"] == 1
    assert f"{workspace['ref']}/region.bed" in proc.stdout


@pytest.mark.parametrize(
    "mode,assays",
    [
        ("bulk_genotyping", ["bulkWGS"]),
        ("single_cell_genotyping", ["scRNA", "scATAC"]),
        ("copytyping_preprocess", ["scATAC"]),
    ],
)
def test_extremity_tsv_feeds_build_segment_bed(workspace, mode, assays):
    """A configured extremity_tsv is an input of build_segment_bed, in every mode."""
    ref = workspace["ref"]
    sheet = (
        workspace["bulk_json"] if mode == "bulk_genotyping" else workspace["sc_json"]
    )
    sample_id = "T1" if mode == "bulk_genotyping" else "S1"
    proc = dryrun(
        workspace,
        sheet,
        sample_id,
        mode,
        assays,
        extra=[f"extremity_tsv={ref}/extremity.tsv"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "cut region_bed arms at SV extremities" in proc.stdout
    assert job_counts(proc.stdout)["build_segment_bed"] == 1
    assert f"{ref}/extremity.tsv" in proc.stdout


def test_extremity_tsv_overrides_a_prebuilt_window_bed(workspace):
    """A pre-built grid cannot honor the cuts, so it is ignored and the windows re-tiled."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"extremity_tsv={ref}/extremity.tsv", f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "is ignored" in proc.stdout
    assert job_counts(proc.stdout).get("build_window_bed", 0) == 1, proc.stdout[-2000:]
    assert "/windows.bed.gz" in proc.stdout


@pytest.mark.parametrize(
    "mode,assays",
    [
        ("bulk_genotyping", ["bulkWGS"]),
        ("single_cell_genotyping", ["scRNA", "scATAC"]),
        ("copytyping_preprocess", ["scATAC"]),
    ],
)
def test_segment_bed_is_built_in_every_mode(workspace, mode, assays):
    """Every mode reads aux/segment.bed, so build_segment_bed is always planned."""
    sheet = (
        workspace["bulk_json"] if mode == "bulk_genotyping" else workspace["sc_json"]
    )
    sample_id = "T1" if mode == "bulk_genotyping" else "S1"
    proc = dryrun(workspace, sheet, sample_id, mode, assays)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "build_segment_bed" in job_counts(proc.stdout)


@pytest.mark.network
def test_repliseq_is_fetched_when_rt_correct(workspace):
    """rt_correct on a Repli-seq build pulls the ENCODE tracks (16 + one liftOver).

    The suite runs with rt_correct off (conftest), so this is the only test that
    resolves UCSC URLs; hence the `network` marker.
    """
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['params_count_reads={"rt_correct": True}'],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts.get("repliseq_bigwig_to_bedgraph") == 16, proc.stdout[-2000:]
    assert counts.get("repliseq_liftover") == 16
    assert ".chm13v2.bedGraph" in proc.stdout


def test_mixed_wgs_wes(workspace):
    """bulkWGS + bulkWES share one window grid and one joint bulk bb dir (WES == WGS)."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # segment BED + exactly one shared window BED (no per-stream / wes_targets)
    assert "build_segment_bed" in counts
    assert counts.get("build_window_bed", 0) == 1, proc.stdout[-2000:]
    assert "/windows.bed.gz" in proc.stdout
    assert "wes_targets" not in proc.stdout
    assert "wgs_windows.bed.gz" not in proc.stdout
    assert "wes_windows.bed.gz" not in proc.stdout
    assert counts.get("rd_correct", 0) == 1, proc.stdout[-2000:]
    assert "/bulk/window.dp.npz" in proc.stdout
    assert "/bulk/window.raw.dp.npz" in proc.stdout
    assert "bulkWGS/window.tsv.gz" not in proc.stdout
    assert "bulkWES/window.tsv.gz" not in proc.stdout
    # one joint binning into a single bb/bulk dir (no per-stream subdir)
    assert "combine_counts" in counts
    assert "/bulk/bb.tsv.gz" in proc.stdout
    assert "bulkWES/bb.tsv.gz" not in proc.stdout


def test_target_bed_adds_target_site_annotation(workspace):
    """target_bed adds one annotation job feeding both rd_correct and combine_counts."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=[f"target_bed={workspace['ref']}/targets.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts.get("annotate_window_targets", 0) == 1, proc.stdout[-2000:]
    assert "/window.target.npz" in proc.stdout


def test_no_target_bed_skips_target_site_annotation(workspace):
    """Without target_bed no annotation job exists, and the RDR SE is still emitted."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "annotate_window_targets" not in job_counts(proc.stdout)
    assert "window.target.npz" not in proc.stdout


@pytest.mark.parametrize(
    "sheet,sample_id,mode,assays",
    [
        ("bulk_json", "T1", "bulk_genotyping", ["bulkWGS"]),
        ("sc_json", "S1", "single_cell_genotyping", ["scRNA", "scATAC"]),
    ],
    ids=["bulk", "single_cell"],
)
def test_prebuilt_windows_are_used_not_built(workspace, sheet, sample_id, mode, assays):
    """A pre-built window_bed is consumed as-is; nothing is tiled."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace[sheet],
        sample_id,
        mode,
        assays,
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    # nothing is tiled; the supplied grid is read straight from its configured path
    assert "build_window_bed" not in counts
    # build_segment_bed still runs: the SNP side and the RD QC overlay both read it
    assert "build_segment_bed" in counts
    assert f"{ref}/window.bed" in proc.stdout


def test_prebuilt_windows_cover_all_assays(workspace):
    """A pre-built window_bed serves WGS + WES alike, so nothing is built."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=[f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" not in counts
    assert f"{ref}/window.bed" in proc.stdout


def test_prebuilt_windows_skip_repliseq(workspace):
    """A pre-built window_bed skips the Repli-seq fetch (do_repliseq active on hg38).

    Network-free: the Repli-seq rules are the only URL-storage inputs here, and they
    are gated inside `if build_windows`, so nothing queries a remote host.
    """
    ref = workspace["ref"]
    sheet = _sheet_with_refvers(workspace, "repliseq.json", ["hg38", "hg38"])
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=hg38", f"window_bed={ref}/window.bed"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    counts = job_counts(proc.stdout)
    assert "repliseq_bigwig_to_bedgraph" not in counts
    assert "build_window_bed" not in counts


def test_scatac_fragments_are_tracked(workspace):
    """The ATAC fragments file is a tracked input, not resolved inside a script."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert proc.returncode == 0
    assert "atac_fragments.tsv.gz" in proc.stdout


def test_visium_spatial_files_are_tracked(workspace):
    """Each Space Ranger spatial/ file is a tracked input of process_rna_anndata."""
    proc = dryrun(
        workspace, workspace["sc_json"], "V1", "single_cell_genotyping", ["VISIUM"]
    )
    assert proc.returncode == 0
    for name in (
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png",
    ):
        assert name in proc.stdout, f"{name} is not a tracked input"


def test_unread_files_key_is_dropped(workspace):
    """A files key the assay never reads is dropped, null or not, and plans the same DAG."""
    sheet = os.path.join(workspace["root"], "null_file.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    doc["samples"][1]["files"]["fragments"] = (
        None  # scATAC-only key on a bulkWGS record
    )
    doc["samples"][0]["files"]["barcodes"] = "/path/to/nowhere.tsv.gz"
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    base = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert job_counts(proc.stdout) == job_counts(base.stdout)


def test_tsv_missing_required_column_fails(workspace):
    """TSV columns are the record keys; a missing one is named, not silently empty."""
    sheet = os.path.join(workspace["root"], "no_col.tsv")
    lines = open(workspace["bulk_tsv"]).read().splitlines()
    header = lines[0].split("\t")
    drop = header.index("reference_version")
    rows = [
        "\t".join(c for i, c in enumerate(ln.split("\t")) if i != drop) for ln in lines
    ]
    with open(sheet, "w") as fh:
        fh.write("\n".join(rows) + "\n")
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "missing required key(s)" in out and "reference_version" in out


def test_record_without_reference_version_fails(workspace):
    """reference_version is required on every record, not optional provenance."""
    sheet = os.path.join(workspace["root"], "no_refver.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    del doc["samples"][1]["reference_version"]
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "missing required key(s)" in out and "reference_version" in out


def _sheet_with_refvers(workspace, name, refvers, source="bulk_json"):
    """Copy a sample sheet, setting each record's reference_version from *refvers*."""
    sheet = os.path.join(workspace["root"], name)
    doc = json.loads(open(workspace[source]).read())
    for rec, refver in zip(doc["samples"], refvers):
        rec["reference_version"] = refver
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    return sheet


def test_chromosome_absent_from_genome_size_fails(workspace):
    """Config chromosomes are checked against the genome once, at DAG build."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["chromosomes=[22,99]"],
    )
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "are not found in" in out and "'99'" in out


def test_config_species_required(workspace):
    """An unset config species is an error, not a silent fallback to human."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species="],
    )
    assert proc.returncode != 0
    assert "species is required" in proc.stdout + proc.stderr


def test_species_reaches_parse_genetic_map(workspace):
    """The shipped config says human, and that is what parse_genetic_map is given."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "parse_genetic_map" in job_counts(proc.stdout)
    mouse = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species=mouse"],
    )
    assert mouse.returncode == 0, mouse.stderr[-1500:]


def test_unsupported_species_warns(workspace):
    """An unknown species warns but runs; it only matters if a gmap needs relabeling."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["species=zebrafish"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "species='zebrafish' is not natively supported" in proc.stdout


def test_config_reference_version_required(workspace):
    """An unset config reference_version is an error, not a silent no-filter."""
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version="],
    )
    assert proc.returncode != 0
    assert "reference_version is required" in proc.stdout + proc.stderr


def test_record_reference_version_alias_matches(workspace):
    """Config and records may spell the build differently; both are canonicalized."""
    sheet = _sheet_with_refvers(workspace, "alias.json", ["T2T-CHM13v2.0", "chm13v2.0"])
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=CHM13"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]


def test_records_of_another_build_are_dropped(workspace):
    """A record on a different build is not selected, leaving nothing to run."""
    sheet = _sheet_with_refvers(workspace, "otherbuild.json", ["hg19", "hg19"])
    proc = dryrun(workspace, sheet, "T1", "bulk_genotyping", ["bulkWGS"])
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "no datasets exist after selection" in out


def test_unrecognized_reference_version_still_selects(workspace):
    """An unsupported build warns but still runs, so long as records match it."""
    sheet = _sheet_with_refvers(
        workspace, "giabv3.json", ["GRCh38-GIABv3", "GRCh38-GIABv3"]
    )
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=["reference_version=GRCh38-GIABv3"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "is not natively supported" in proc.stdout


def test_unknown_sample_id_fails(workspace):
    """An unmatched sample_id is an error, not an empty DAG."""
    proc = dryrun(
        workspace, workspace["bulk_json"], "nope", "bulk_genotyping", ["bulkWGS"]
    )
    assert proc.returncode != 0
    assert "no datasets exist after selection" in proc.stdout + proc.stderr


def test_dataset_ids_selects_a_subset(workspace):
    """dataset_ids drops the sample_id's other datasets from the DAG."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=(f"dataset_ids={json.dumps(['N1', 'D1'])}",),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts["run_mosdepth"] == 2, proc.stdout[-2000:]
    assert counts["merge_pileup_counts"] == 2
    assert "bulkWES_E1" not in proc.stdout


def test_dataset_ids_empty_runs_every_dataset(workspace):
    """The default [] keeps all three datasets."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert job_counts(proc.stdout)["run_mosdepth"] == 3


def test_dataset_ids_unknown_fails(workspace):
    """A dataset_id absent from the selection is an error, not a silent drop."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=(f"dataset_ids={json.dumps(['N1', 'nope'])}",),
    )
    assert proc.returncode != 0
    assert "dataset_ids not found" in proc.stdout + proc.stderr


def test_dataset_ids_on_another_build_fails(workspace):
    """A selected dataset_id with no record on the run's build names the build."""
    sheet = _sheet_with_refvers(workspace, "selbuild.json", ["chm13v2", "hg19"])
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=(f"dataset_ids={json.dumps(['N1', 'D1'])}",),
    )
    assert proc.returncode != 0
    out = proc.stdout + proc.stderr
    assert "dataset_ids have no record on reference_version" in out
    assert "D1 is on ['hg19']" in out


def test_dataset_ids_keeps_the_matching_build_row(workspace):
    """One dataset_id with a row per build keeps only the row on the run's build."""
    sheet = os.path.join(workspace["root"], "dupbuild.json")
    doc = json.loads(open(workspace["bulk_json"]).read())
    other = dict(doc["samples"][1])
    other["reference_version"] = "hg19"
    doc["samples"].append(other)
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=(f"dataset_ids={json.dumps(['N1', 'D1'])}",),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert job_counts(proc.stdout)["run_mosdepth"] == 2, proc.stdout[-2000:]


def test_dataset_ids_excluding_an_rdr_base_fails(workspace):
    """Dropping the normal a tumor normalizes against is caught at parse time."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS", "bulkWES"],
        extra=(f"dataset_ids={json.dumps(['D1'])}",),
    )
    assert proc.returncode != 0
    assert "is not among the selected datasets" in proc.stdout + proc.stderr


@pytest.mark.parametrize("mode", ["bulk_genotyping", "single_cell_genotyping"])
def test_het_snp_vcf_skips_genotyping(workspace, mode):
    """A supplied het_snp_vcf replaces genotyping in every mode."""
    sheet, sample_id, assays = (
        ("bulk_json", "T1", ["bulkWGS"])
        if mode == "bulk_genotyping"
        else ("sc_json", "S1", ["scRNA", "scATAC"])
    )
    ref = workspace["ref"]
    base = dryrun(workspace, workspace[sheet], sample_id, mode, assays)
    with_vcf = dryrun(
        workspace,
        workspace[sheet],
        sample_id,
        mode,
        assays,
        extra=[f"het_snp_vcf={ref}/het_snps.vcf.gz", "het_snp_vcf_phased=False"],
    )
    assert base.returncode == 0 and with_vcf.returncode == 0, with_vcf.stderr[-1500:]
    genotype_rules = {
        r for r in job_counts(base.stdout) if r.startswith("genotype_snps")
    }
    assert genotype_rules, "baseline should genotype"
    assert not genotype_rules & set(job_counts(with_vcf.stdout))
    # het_snp_vcf_phased=false: split per chromosome, then phased
    counts = job_counts(with_vcf.stdout)
    assert "split_het_snp_vcf" in counts
    assert any(r.startswith("phase_snps") for r in counts)


def test_phased_het_snp_vcf_skips_phasing(workspace):
    """het_snp_vcf_phased=true replaces genotyping and phasing."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["bulk_json"],
        "T1",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=[f"het_snp_vcf={ref}/het_snps.vcf.gz", "het_snp_vcf_phased=True"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    counts = job_counts(proc.stdout)
    for rule in counts:
        assert not rule.startswith(("genotype_snps", "phase_snps", "split_het_snp_vcf"))
    assert "concat_and_extract_phased_het_snps" not in counts
    assert "pileup_snps_bulk_bcftools_chrom" in counts


def test_copytyping_requires_phased_vcf(workspace):
    """copytyping_preprocess never phases, so its het_snp_vcf must be declared phased."""
    ref = workspace["ref"]
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "copytyping_preprocess",
        ["scATAC"],
        extra=[
            f"het_snp_vcf={ref}/het_snps.vcf.gz",
            "het_snp_vcf_phased=False",
            f"bb_file={ref}/bb.tsv.gz",
        ],
    )
    assert proc.returncode != 0
    assert "het_snp_vcf must be phased" in proc.stdout + proc.stderr


def _longread_sheet(workspace, name):
    """Bulk sheet with a short-read normal, a long-read normal and a long-read tumor."""
    ref = workspace["ref"]

    def files(stem):
        return {
            "alignment": f"{ref}/{stem}.bam",
            "alignment_index": f"{ref}/{stem}.bam.bai",
        }

    doc = {
        "version": 1,
        "samples": [
            {
                "sample_id": "T1",
                "dataset_id": "N1",
                "assay_type": "bulkWGS",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": files("normal"),
            },
            {
                "sample_id": "T1",
                "dataset_id": "L1",
                "assay_type": "bulkWGS-lr",
                "sample_type": "normal",
                "reference_version": "chm13v2",
                "files": files("normal"),
            },
            {
                "sample_id": "T1",
                "dataset_id": "D1",
                "assay_type": "bulkWGS-lr",
                "sample_type": "tumor",
                "reference_version": "chm13v2",
                "files": files("tumor"),
            },
        ],
    }
    sheet = os.path.join(workspace["root"], name)
    with open(sheet, "w") as fh:
        json.dump(doc, fh)
    return sheet


def test_longphase_auto_picks_the_long_read_normal(workspace):
    """phaser=longphase with no phase_dataset_ids co-phases the long-read normal."""
    sheet = _longread_sheet(workspace, "longread.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS", "bulkWGS-lr"],
        extra=["phaser=longphase"],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "phase_dataset_ids: ['L1']" in proc.stdout
    counts = job_counts(proc.stdout)
    assert "phase_snps_longphase" in counts
    assert "phase_snps_eagle" not in counts and "parse_genetic_map" not in counts


def test_phase_dataset_ids_must_be_long_read(workspace):
    """A short-read dataset named in phase_dataset_ids is rejected, not silently kept."""
    sheet = _longread_sheet(workspace, "longread_shortread.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS", "bulkWGS-lr"],
        extra=["phaser=longphase", 'phase_dataset_ids=["N1"]'],
    )
    assert proc.returncode != 0
    assert "needs long reads" in proc.stdout + proc.stderr


def test_phase_dataset_ids_outside_selection_fails(workspace):
    """An id in the sheet but dropped by assay_types is an error, not an auto-pick."""
    sheet = _longread_sheet(workspace, "longread_selection.json")
    proc = dryrun(
        workspace,
        sheet,
        "T1",
        "bulk_genotyping",
        ["bulkWGS-lr"],
        extra=["phaser=longphase", 'phase_dataset_ids=["N1"]'],
    )
    assert proc.returncode != 0
    assert "phase_dataset_ids not found in the records" in proc.stdout + proc.stderr


def test_genotype_dataset_ids_outside_selection_fails(workspace):
    """bulkWES E1 is in the sheet but dropped by assay_types, so naming it fails."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['genotype_dataset_ids=["E1"]'],
    )
    assert proc.returncode != 0
    assert "genotype_dataset_ids not found in the records" in proc.stdout + proc.stderr


def test_genotype_dataset_ids_in_selection_is_used(workspace):
    """The same id check accepts a dataset the selection kept."""
    proc = dryrun(
        workspace,
        workspace["bulk_mixed_json"],
        "MX",
        "bulk_genotyping",
        ["bulkWGS"],
        extra=['genotype_dataset_ids=["N1"]'],
    )
    assert proc.returncode == 0, proc.stderr[-1500:]
    assert "genotype_dataset_ids: ['N1']" in proc.stdout


def test_single_cell_skips_repliseq(workspace):
    """`do_repliseq` is the one thing build_window_bed still gates on the mode.

    GC and MAP are annotated identically in every mode, but the ENCODE bigWig fetch +
    liftOver must stay out of a single-cell DAG: they are 15 UCSC downloads feeding a
    covariate only rd_correct (bulk) reads.
    """
    sheet = _sheet_with_refvers(
        workspace, "sc_hg38.json", ["hg38"] * 4, source="sc_json"
    )
    proc = dryrun(
        workspace,
        sheet,
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
        extra=["reference_version=hg38"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "build_window_bed" in counts
    assert "repliseq_bigwig_to_bedgraph" not in counts
    assert "repliseq_liftover" not in counts


def test_scdna_runs_as_bulk(workspace):
    """scDNA takes the bulk path: mosdepth, bcftools pileup, one joint bb set."""
    proc = dryrun(
        workspace,
        workspace["bulk_scdna_json"],
        "SD",
        "bulk_genotyping",
        ["bulkWGS", "scDNA"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts.get("rd_correct", 0) == 1, proc.stdout[-2000:]
    assert "scDNA/out_mosdepth/C1.regions.bed.gz" in proc.stdout, proc.stdout[-2000:]
    assert "scDNA_C1/bcftools.counts.tsv.gz" in proc.stdout, proc.stdout[-2000:]
    assert "combine_counts" in counts
    # no single-cell rule is reachable: the reads are pooled, cell tags never read
    assert "pileup_snps_nonbulk_mode1a" not in counts
    assert "/bb.Xcount.npz" not in proc.stdout

    # the pooling is warned about at parse time; a dry run has no run-log file
    # handler, so logging_snakemake falls back to stdout
    assert "WARNING: single-cell dataset(s) C1 (scDNA)" in proc.stdout + proc.stderr
    # the scDNA pileup collapses the per-barcode @RG SM tags; the bulkWGS one does not
    assert "--ignore-RG" in proc.stdout, proc.stdout[-2000:]


def test_scdna_genotyping_ignores_read_groups(workspace):
    """Genotyping off the scDNA BAM pools its per-barcode @RG SM tags into one sample."""
    proc = dryrun(
        workspace,
        workspace["bulk_scdna_json"],
        "SD",
        "bulk_genotyping",
        ["bulkWGS", "scDNA"],
        extra=(f"genotype_dataset_ids={json.dumps(['C1'])}",),
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "--ignore-RG" in proc.stdout, proc.stdout[-2000:]
    assert "per-barcode @RG SM tags collapse" in proc.stdout + proc.stderr


def test_scdna_genotyping_cannot_pool_with_another_dataset(workspace):
    """--ignore-RG is per input file, so a per-cell library must be genotyped alone."""
    proc = dryrun(
        workspace,
        workspace["bulk_scdna_json"],
        "SD",
        "bulk_genotyping",
        ["bulkWGS", "scDNA"],
        extra=(f"genotype_dataset_ids={json.dumps(['N1', 'C1'])}",),
    )
    assert proc.returncode != 0
    assert "genotype exactly one dataset" in proc.stdout + proc.stderr


def test_scdna_rejected_outside_bulk(workspace):
    """scDNA is not a non-bulk assay, so the other two modes have nothing to run."""
    for mode in ("single_cell_genotyping", "copytyping_preprocess"):
        proc = dryrun(workspace, workspace["bulk_scdna_json"], "SD", mode, ["scDNA"])
        assert proc.returncode != 0, mode
        assert "none valid for workflow_mode" in proc.stdout + proc.stderr


def test_copytyping_builds_the_window_grid(workspace):
    """copytyping_preprocess bins onto its own bb_file, but the unit level needs windows."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "copytyping_preprocess",
        ["scRNA", "scATAC"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert counts.get("build_window_bed", 0) == 1, proc.stdout[-2000:]
    # the two rules cover disjoint assays, so both run in one DAG without ambiguity
    assert counts.get("combine_counts_fixed_bins", 0) == 1, proc.stdout[-2000:]
    assert counts.get("combine_counts_fixed_bins_rna", 0) == 1, proc.stdout[-2000:]
    assert "/unit/scATAC/window.Xcount.npz" in proc.stdout, proc.stdout[-2000:]
    assert "/unit/scATAC/snp.tsv.gz" in proc.stdout, proc.stdout[-2000:]
    assert "/unit/scRNA/gene.Xcount.npz" in proc.stdout, proc.stdout[-2000:]


def test_copytyping_rna_units_are_genes(workspace):
    """An RNA assay's unit is the gene, so it takes the _rna rule and writes no windows."""
    proc = dryrun(
        workspace, workspace["sc_json"], "S1", "copytyping_preprocess", ["scRNA"]
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    counts = job_counts(proc.stdout)
    assert "combine_counts_fixed_bins_rna" in counts, proc.stdout[-2000:]
    assert "combine_counts_fixed_bins" not in counts, proc.stdout[-2000:]
    assert "/unit/scRNA/gene.Xcount.npz" in proc.stdout, proc.stdout[-2000:]
    assert "/unit/scRNA/window.tsv.gz" not in proc.stdout


def test_single_cell_binning_reads_the_window_grid(workspace):
    """combine_counts_nonbulk takes the window BED as its fixed bins."""
    proc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "combine_counts_nonbulk" in job_counts(proc.stdout)
    assert "/windows.bed.gz" in proc.stdout


def test_unit_level_outputs(workspace):
    """Every mode that bins writes the un-binned unit level next to the bbs."""
    bulk = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert bulk.returncode == 0, bulk.stderr[-2000:]
    for suffix in (
        "snp.tsv.gz",
        "snp.Tallele.npz",
        "window.tsv.gz",
        "window.depth.npz",
        "sample_ids.tsv",
    ):
        assert f"/unit/bulk/{suffix}" in bulk.stdout, bulk.stdout[-2000:]

    sc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert sc.returncode == 0, sc.stderr[-2000:]
    assert "/unit/scATAC/window.Xcount.npz" in sc.stdout, sc.stdout[-2000:]
    assert "/unit/scRNA/gene.Xcount.npz" in sc.stdout, sc.stdout[-2000:]
    # a gene is indivisible, so RNA has no window level and ATAC no gene level
    assert "/unit/scRNA/window.tsv.gz" not in sc.stdout
    assert "/unit/scATAC/gene.tsv.gz" not in sc.stdout


def test_qc_pdf_names(workspace):
    """The QC PDFs each rule declares, and the bulk binning split into two files."""
    bulk = dryrun(
        workspace, workspace["bulk_json"], "T1", "bulk_genotyping", ["bulkWGS"]
    )
    assert bulk.returncode == 0, bulk.stderr[-2000:]
    for name in (
        "qc/genotype_snps.pdf",
        "qc/rd_correction.pdf",
        "qc/phase_and_concat.bulk.pdf",
        "qc/combine_counts.stats.bulk.MSR100.pdf",
        "qc/combine_counts.1d2d.bulk.MSR100.pdf",
    ):
        assert name in bulk.stdout, bulk.stdout[-2000:]

    sc = dryrun(
        workspace,
        workspace["sc_json"],
        "S1",
        "single_cell_genotyping",
        ["scRNA", "scATAC"],
    )
    assert sc.returncode == 0, sc.stderr[-2000:]
    assert "qc/genotype_snps.pdf" in sc.stdout, sc.stdout[-2000:]
    for name in (
        "qc/combine_counts.scATAC.MSR100.pdf",
        "qc/combine_counts.stats.scATAC.pdf",
        "qc/combine_counts.stats.scATAC.tsv",
        "qc/combine_counts.stats.scRNA.pdf",
    ):
        assert name in sc.stdout, sc.stdout[-2000:]

    ct = dryrun(
        workspace, workspace["sc_json"], "S1", "copytyping_preprocess", ["scATAC"]
    )
    assert ct.returncode == 0, ct.stderr[-2000:]
    # the rule suffix is dropped: the bb page is combine_counts.{assay}.pdf
    assert "qc/combine_counts.scATAC.pdf" in ct.stdout, ct.stdout[-2000:]
    assert "qc/combine_counts_fixed_bins." not in ct.stdout
    for name in (
        "qc/combine_counts.stats.scATAC.pdf",
        "qc/combine_counts.stats.scATAC.tsv",
    ):
        assert name in ct.stdout, ct.stdout[-2000:]
