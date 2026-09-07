"""Read the pipeline's input files and write its tabular outputs.

Last update: 2026-09-02

Grouped by format, each block reading before writing:
- Universal: read_chrom_sizes, symlink_files
- GTF: read_GTF
- VCF: read_VCF, write_VCF
- BED: read_BED, read_segment_bed, read_window_bed, read_mosdepth_bed, read_bedgraph,
  read_extremity_tsv
- Pipeline: read_bcftools_pileup_counts, read_allele_mat, read_snp_mats, read_barcodes,
  read_barcodes_by_dataset, read_chunks_from_atac_fragments, uniquify_var_names,
  read_10x_ranger_scRNA, read_10x_ranger_spatial, write_snp_info, write_bb_file,
  write_sample_ids
"""

import logging
import os
import subprocess
import tempfile
import warnings
from collections import OrderedDict

import pandas as pd
import numpy as np

from const import (
    GTF_COLUMNS,
    RANGER_MATRIX_H5,
    RANGER_SPATIAL_DIR,
    VCF_COLUMNS,
    VCF_SAMPLE_COLUMNS,
)
from utils import add_chr_prefix, log_ratios, sort_chroms, sort_df_chr


# --------------------------------------------------------------------------
# Universal
# --------------------------------------------------------------------------


def read_chrom_sizes(sz_file: str):
    """Read a two-column chromosome-sizes file.

    Args:
        sz_file: Tab-separated ``chrom<TAB>size``.

    Returns:
        OrderedDict of chromosome name to length, in file order.
    """
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def symlink_files(srcs, dsts):
    """Point each destination at its source, as a relative symlink.

    For a rule whose output is its input unchanged: linking keeps one copy of the bytes
    and records the provenance in the link. Relative so the run directory stays movable,
    and an existing destination is replaced so a rerun does not fail on it. The source
    must outlive the link, so it may not be a ``temp()`` output.

    Args:
        srcs: Existing files, in destination order.
        dsts: Paths to create.
    """
    assert len(srcs) == len(dsts), f"symlink_files: {len(srcs)} src for {len(dsts)} dst"
    for src, dst in zip(srcs, dsts):
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.relpath(src, os.path.dirname(dst) or "."), dst)


# --------------------------------------------------------------------------
# GTF
# --------------------------------------------------------------------------


def read_GTF(gtf_file: str, feature_types=("gene",), id_col="gene_id"):
    """Parse a GTF once and split it by feature type.

    Contigs are chr-normalized and coordinates become 0-based half-open; genes are
    deduplicated by their id.

    Args:
        gtf_file: Path to a GTF annotation file (optionally gzipped).
        feature_types: Feature types to extract, e.g. ``("gene", "exon")``.
        id_col: Name of the gene-id column in the output.

    Returns:
        ``{feature_type: DataFrame}`` with ``#CHR``, ``START``, ``END``, *id_col*.

    Notes/References:
        Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format4
    """
    wanted = list(feature_types)
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLUMNS,
        dtype={"seqname": str},
        low_memory=False,
    )
    gtf = gtf.loc[
        gtf["feature"].isin(wanted),
        ["feature", "seqname", "start", "end", "attributes"],
    ]
    flat = pd.DataFrame(
        {
            "feature": gtf["feature"].values,
            "#CHR": add_chr_prefix(gtf["seqname"]).values,
            "START": gtf["start"].values - 1,  # GTF is 1-based -> 0-based
            "END": gtf["end"].values,  # GTF end is inclusive -> half-open
            id_col: gtf["attributes"]
            .str.extract(r'gene_id "([^"]+)"', expand=False)
            .values,
        }
    )
    out = {}
    for feature_type in wanted:
        sub = flat.loc[flat["feature"] == feature_type].drop(columns="feature")
        if feature_type == "gene":
            sub = sub.drop_duplicates(id_col, keep="first")
        out[feature_type] = sub.reset_index(drop=True)
    return out


# --------------------------------------------------------------------------
# VCF
# --------------------------------------------------------------------------


def read_VCF(
    vcf_file: str,
    addchr=True,
    addkey=False,
    snps_presorted=False,
    add_pos0=False,
    add_phase1=False,
    read_AD=False,
    required_cols=[],
):
    """Read a VCF into a DataFrame, exploding its INFO and FORMAT fields into columns.

    Args:
        vcf_file: Path to the (optionally gzipped) VCF.
        addchr: Prepend ``chr`` to contigs named without it.
        addkey: Add ``KEY`` (``#CHROM_POS``).
        snps_presorted: Skip the genomic sort.
        add_pos0: Add ``POS0`` (0-based).
        add_phase1: Add ``PHASE``, the second GT allele.
        read_AD: Split the comma-joined ``AD`` into ``REF_COUNT`` and ``ALT_COUNT``, which
            requires a bi-allelic single-sample VCF.
        required_cols: Columns the caller needs, asserted once the INFO/FORMAT keys are
            exploded; the frame may carry more.

    Returns:
        DataFrame with the 8 fixed VCF columns, ``#CHR``, ``RAW_SNP_DF_IDX`` and one
        column per INFO/FORMAT key; None when the file has no records.
    """
    snps = pd.read_csv(
        vcf_file, comment="#", sep="\t", header=None, dtype={0: "string"}
    )
    if snps.empty:
        return None
    ncols = snps.shape[1]
    assert ncols == 8 or ncols >= 10, (
        f"VCF file, expected 8 or >=10 columns, got {ncols}"
    )
    colnames = list(VCF_COLUMNS)
    if ncols >= 10:
        colnames += VCF_SAMPLE_COLUMNS
        snps = snps.iloc[:, :10].copy()
    snps.columns = colnames
    snps["POS"] = snps["POS"].astype(np.int64)
    snps["RAW_SNP_DF_IDX"] = np.arange(len(snps))
    if addchr and not str(snps["#CHROM"].iloc[0]).startswith("chr"):
        snps["#CHROM"] = "chr" + snps["#CHROM"].astype(str)

    snps["#CHROM"] = snps["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    chrom_order = sort_chroms(snps["#CHROM"].unique().tolist())
    snps["#CHROM"] = pd.Categorical(
        snps["#CHROM"], categories=chrom_order, ordered=True
    )
    if not snps_presorted:
        snps = snps.sort_values(["#CHROM", "POS"], kind="mergesort")

    snps["#CHR"] = snps["#CHROM"]

    # parse INFO column
    info_kvs = (
        snps["INFO"]
        .fillna("")
        .str.split(";")
        .explode()
        .loc[lambda s: s.ne("")]
        .to_frame("kv")
    )
    kv = info_kvs["kv"].str.split("=", n=1, expand=True)
    info_kvs["key"] = kv[0]
    info_kvs["val"] = kv[1] if kv.shape[1] > 1 else None
    info_kvs["val"] = info_kvs["val"].fillna(True)  # INFO flags (no '=') -> True
    info_kvs["row"] = info_kvs.index

    info_wide = info_kvs.pivot_table(
        index="row", columns="key", values="val", aggfunc="first"
    )
    snps = snps.join(info_wide)

    # parse FORMAT column
    if "FORMAT" in snps.columns:
        fmt_keys = snps["FORMAT"].fillna("").str.split(":")
        samp_vals = snps["SAMPLE"].fillna("").str.split(":")

        fmt_long = pd.DataFrame(
            {"key": fmt_keys.explode(), "val": samp_vals.explode()}
        ).dropna(subset=["key"])
        fmt_long["row"] = fmt_keys.explode().index  # original variant row index

        fmt_wide = fmt_long.pivot_table(
            index="row", columns="key", values="val", aggfunc="first"
        )
        snps = snps.drop(columns=["FORMAT", "SAMPLE"])
        # A key can appear in both INFO and FORMAT (e.g. AD/DP when a phaser
        # such as longphase preserves the input INFO and FORMAT fields). Detect
        # any such overlap and let the per-sample FORMAT value take precedence,
        # dropping the INFO-derived duplicate to avoid a "columns overlap" join
        # error. No overlap (eagle/shapeit emit only GT) -> no-op.
        dup_cols = snps.columns.intersection(fmt_wide.columns)
        if len(dup_cols):
            snps = snps.drop(columns=dup_cols)
        snps = snps.join(fmt_wide)

    if addkey:
        snps["KEY"] = snps["#CHROM"].astype(str) + "_" + snps["POS"].astype(str)
    if add_pos0:
        snps["POS0"] = snps["POS"] - 1
    if add_phase1:
        snps["PHASE"] = snps["GT"].str[2].astype(np.int8).to_numpy()
    if read_AD:
        allele_depths = snps["AD"].astype(str).str.split(",", expand=True)
        assert allele_depths.shape[1] >= 2, (
            f"{vcf_file}: FORMAT/AD must hold at least REF,ALT; got "
            f"{allele_depths.shape[1]} field(s). bcftools call needs --keep-alts."
        )
        for name, field in (("REF_COUNT", 0), ("ALT_COUNT", 1)):
            snps[name] = (
                pd.to_numeric(allele_depths[field], errors="coerce")
                .fillna(0)
                .astype(np.int64)
            )
    snps = snps.reset_index(drop=True)

    missing = [c for c in required_cols if c not in snps.columns]
    assert not missing, f"{vcf_file}: missing required VCF column(s) {missing}"
    return snps


VCF_HEADER_LINES = [
    "##fileformat=VCFv4.2\n",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth, REF+ALT">\n',
    '##INFO=<ID=AD,Number=1,Type=Integer,Description="Allele Depth for ALT allele">\n',
    '##INFO=<ID=OTH,Number=1,Type=Integer,Description="Allele Depth other than REF and ALT">\n',
]


def write_VCF(snps, out_file, chrname, chrom_length=None, create_index=True):
    """Write one chromosome's genotyped SNPs as a bgzipped VCF.

    The counterpart of :func:`read_VCF` for the VCFs this pipeline writes itself: one
    contig, one sample, ``GT`` in FORMAT and the ``DP``/``AD``/``OTH`` counts in INFO.
    ``#CHROM`` is stamped with *chrname*, so the caller decides the contig spelling once
    and the frame it passes need not carry it.

    Args:
        snps: Rows for this contig, carrying VCF_COLUMNS + VCF_SAMPLE_COLUMNS; may be
            empty, which writes a header-only VCF.
        out_file: Output path, ending ``.vcf.gz``.
        chrname: Contig name to write, in the run's input spelling.
        chrom_length: Contig length for the ``##contig`` header; omitted when None.
        create_index: Also write the tabix index.

    Returns:
        The DataFrame written.
    """
    assert out_file.endswith(".vcf.gz"), f"write_VCF: {out_file} is not a .vcf.gz"
    cols = list(VCF_COLUMNS) + VCF_SAMPLE_COLUMNS
    missing = [c for c in cols if c not in snps.columns]
    assert not missing, f"VCF table, missing column(s) {missing}"

    contig = f"##contig=<ID={chrname}"
    contig += ">\n" if chrom_length is None else f",length={chrom_length}>\n"
    out = snps[cols].copy()
    out["#CHROM"] = chrname

    plain_file = out_file[:-3]
    with open(plain_file, "w") as fd:
        fd.writelines(VCF_HEADER_LINES[:1] + [contig] + VCF_HEADER_LINES[1:])
        fd.write("\t".join(cols) + "\n")
        out.to_csv(fd, sep="\t", index=False, header=False)

    # bgzip removes plain_file
    subprocess.run(["bgzip", "-f", plain_file], check=True)
    if create_index:
        subprocess.run(["tabix", "-f", "-p", "vcf", out_file], check=True)
    return out


# --------------------------------------------------------------------------
# BED
# --------------------------------------------------------------------------


def read_BED(bed_file: str, addchr=True, col_id="region_id"):
    """Read BED file.

    Args:
        bed_file: Path to a BED file, at least 3 columns.
        addchr: Prepend ``chr`` to contigs if contigs are not chr-prefix.
        col_id: Name of the id column, read from column 4 or derived.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END`` and *col_id*.

    Notes/References:
        Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    df = pd.read_table(bed_file, sep="\t", header=None, dtype={0: "string"})
    assert len(df.columns) >= 3, (
        f"{bed_file}: BED file, expected >=3 columns, got {len(df.columns)}"
    )
    has_name = len(df.columns) >= 4
    df = df.iloc[:, : 4 if has_name else 3].copy()
    df.columns = ["#CHR", "START", "END"] + ([col_id] if has_name else [])
    if not str(df["#CHR"].iloc[0]).startswith("chr") and addchr:
        df["#CHR"] = "chr" + df["#CHR"].astype(str)
    if not has_name:
        df[col_id] = (
            df["#CHR"] + ":" + df["START"].astype(str) + "-" + df["END"].astype(str)
        )
    return df


def read_segment_bed(bed_file: str, addchr=True):
    """Read custom BED file with 5th column defines segment ID.

    Args:
        bed_file: Path to the custom BED, 5 columns.
        addchr: Prepend ``chr`` to contigs named without it.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, ``region_id``, ``seg_id``.
    """
    df = read_BED(bed_file, addchr=addchr, col_id="region_id")
    df["seg_id"] = pd.read_table(
        bed_file, sep="\t", header=None, usecols=[4], dtype="string"
    ).iloc[:, 0]
    return df


def read_window_bed(bed_file, chroms=None, keep_covariates=False):
    """Read the window BED and sort by genomic positions.

    Args:
        bed_file: Path to the window BED (headered TSV).
        chroms: Keep only these contigs; ``None`` keeps every row.
        keep_covariates: Also carry optional ``GC``/``MAP``/``REPLI`` columns.

    Returns:
        Pandas DataFrame with ``#CHR``, ``START``, ``END``, ``region_id``, ``seg_id``, ``bin_id``.
    """
    bin_df = pd.read_table(bed_file, sep="\t", dtype={"#CHR": str})
    cols = ["#CHR", "START", "END", "region_id"]
    cols += [c for c in ("seg_id",) if c in bin_df.columns]
    if keep_covariates:
        cols += [c for c in ("GC", "MAP", "REPLI") if c in bin_df.columns]
    bin_df = bin_df[cols].copy()
    bin_df["#CHR"] = add_chr_prefix(bin_df["#CHR"])
    if chroms is not None:
        n_raw = len(bin_df)
        keep = bin_df["#CHR"].isin(chroms)
        dropped = (bin_df["END"] - bin_df["START"])[~keep]
        bin_df = bin_df[keep]
        log_ratios(
            f"{bed_file}, bins off the run's {len(chroms)} contigs",
            n_raw - len(bin_df),
            n_raw,
            dropped,
            "bin",
        )
    bin_df = sort_df_chr(bin_df, ch="#CHR", pos="START").reset_index(drop=True)
    if "seg_id" not in bin_df.columns:
        bin_df["seg_id"] = bin_df["region_id"]
    bin_df["bin_id"] = np.arange(len(bin_df))
    return bin_df


def read_mosdepth_bed(mosdepth_bed: str, addchr=True, value_col="DEPTH"):
    """Read a mosdepth ``--by`` regions BED: headerless ``#CHR START END <value>``.

    One row per window of the BED mosdepth was given, in BAM ``@SQ`` order.
    ``count_read_starts_chrom`` writes the same four-column shape, so *value_col* names
    the fourth column: ``DEPTH`` for mosdepth, ``COUNT`` for the read-start counts.

    Args:
        mosdepth_bed: Path to ``{dataset_id}.regions.bed.gz`` or ``.rdcount.bed.gz``.
        addchr: Prepend ``chr`` to contigs named without it.
        value_col: Name to give the fourth column.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, *value_col*.
    """
    df = pd.read_table(
        mosdepth_bed,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", value_col],
        dtype={"#CHR": str},
    )
    if addchr:
        df["#CHR"] = add_chr_prefix(df["#CHR"])
    return df


def read_bedgraph(bg_file: str, chroms=None):
    """Read a bedGraph track: ``chrom start end value``, 0-based half-open.

    Args:
        bg_file: Path to the (optionally gzipped) bedGraph.
        chroms: Keep only these contigs; ``None`` keeps every row.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, ``signal``, reindexed from 0.

    Notes/References:
        Format: https://genome.ucsc.edu/goldenPath/help/bedgraph.html
    """
    df = pd.read_csv(
        bg_file,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "signal"],
        dtype={"#CHR": str, "START": np.int64, "END": np.int64, "signal": np.float64},
    )
    df["#CHR"] = add_chr_prefix(df["#CHR"])
    if chroms is not None:
        df = df[df["#CHR"].isin(chroms)]
    return df.reset_index(drop=True)


def read_extremity_tsv(ext_file: str, chroms=None):
    """Read an SV-extremity file: a headered TSV, one breakpoint per row.

    Only the contig and the breakpoint position are read; any other column an upstream
    caller writes (event id, support counts, VAF) is ignored. The position is taken from
    ``POS0`` when the file has one, else derived from a 1-based ``POS``. Rows are
    deduplicated, since two events sharing a breakpoint cut a segment once.

    Args:
        ext_file: Path to the (optionally gzipped) extremity TSV.
        chroms: Keep only these contigs; ``None`` keeps every row.

    Returns:
        DataFrame with ``#CHR`` and ``POS0``, sorted, reindexed from 0.
    """
    df = pd.read_table(ext_file, sep="\t", dtype={"#CHR": str})
    assert "#CHR" in df.columns, (
        f"extremity file, no `#CHR` column: {ext_file}, got {df.columns.tolist()}"
    )
    assert "POS0" in df.columns or "POS" in df.columns, (
        f"extremity file, no `POS0` or `POS` column: {ext_file}, "
        f"got {df.columns.tolist()}"
    )
    pos = df["POS0"] if "POS0" in df.columns else df["POS"] - 1
    df = pd.DataFrame(
        {"#CHR": add_chr_prefix(df["#CHR"]), "POS0": pos.astype(np.int64)}
    )
    if chroms is not None:
        df = df[df["#CHR"].isin(chroms)]
    df = df.drop_duplicates()
    return sort_df_chr(df, ch="#CHR", pos="POS0")


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------


def read_bcftools_pileup_counts(tsv_file: str, parent_alt_by_key: dict):
    """Read a bcftools per-locus AD table as pseudobulk depth/alt count matrices.

    Counts the parent ALT allele only, looked up in the locus's own ALT list, so
    ``DP = ref + alt`` and a downstream ``REF = DP - ALT`` is exact.

    Args:
        tsv_file: Path to the (optionally gzipped) counts TSV.
        parent_alt_by_key: Parent ALT allele keyed by ``#CHROM_POS``.

    Returns:
        ``(snps, tot_mtx, ad_mtx)``: *snps* carries ``KEY`` (matching ``read_VCF``) and
        ``RAW_SNP_DF_IDX`` (file row order); the matrices are ``(len(snps), 1)`` csr,
        shaped for ``map_allele_mat_to_snps``.
    """
    from scipy.sparse import csr_matrix  # scipy is not a runner-env dependency

    df = pd.read_csv(
        tsv_file,
        sep="\t",
        header=None,
        names=["#CHROM", "POS", "REF", "ALT", "AD"],
        dtype={"#CHROM": "string", "REF": "string", "ALT": "string", "AD": "string"},
    )
    chrom = add_chr_prefix(df["#CHROM"]).str.replace("^chrMT$", "chrM", regex=True)
    keys = chrom + "_" + df["POS"].astype(np.int64).astype(str)
    alt_lists = df["ALT"].str.split(",")
    ad_lists = df["AD"].str.split(",").apply(lambda xs: [int(x) for x in xs])
    parent_alts = keys.map(parent_alt_by_key)

    ref = np.array([ad[0] if ad else 0 for ad in ad_lists], dtype=np.int64)
    alt = np.array(
        [
            ad[1 + alts.index(pa)] if pa in alts else 0
            for ad, alts, pa in zip(ad_lists, alt_lists, parent_alts)
        ],
        dtype=np.int64,
    )
    snps = pd.DataFrame({"KEY": keys.to_numpy(), "RAW_SNP_DF_IDX": np.arange(len(df))})
    # a csr built from a dense column drops the zeros itself
    tot_mtx = csr_matrix((ref + alt).reshape(-1, 1))
    ad_mtx = csr_matrix(alt.reshape(-1, 1))
    return snps, tot_mtx, ad_mtx


def read_allele_mat(npz_file, mat_dtype=None):
    """Read one SNP-level allele matrix, dense or sparse, as it was written.

    Bulk writes a dense ``mat`` key, single-cell a ``scipy.sparse`` archive; the keys
    tell them apart.

    Args:
        npz_file: Path to the ``.npz``.
        mat_dtype: Cast the matrix to this dtype; ``None`` keeps the stored one.

    Returns:
        ``np.ndarray`` for a dense file, ``scipy.sparse.csr_matrix`` for a sparse one.
    """
    from scipy.sparse import load_npz  # scipy is not a runner-env dependency

    with np.load(npz_file) as npz:
        is_dense = "mat" in npz.files
        mat = npz["mat"] if is_dense else load_npz(npz_file)
    return mat if mat_dtype is None else mat.astype(mat_dtype)


def read_snp_mats(snp_info_file, tot_file, a_file, b_file, mat_dtype=None):
    """Read a SNP table and its T/A/B allele matrices, in file row order.

    Args:
        snp_info_file: ``snps.tsv.gz`` from phase_and_concat.
        tot_file, a_file, b_file: the total / A-allele / B-allele ``.npz``.
        mat_dtype: Cast the three matrices to this dtype; ``None`` keeps the stored one,
            which for a sparse file avoids copying its ``data`` array.

    Returns:
        ``(snps, tot_mtx, a_mtx, b_mtx)``; the matrices are dense or sparse per
        ``read_allele_mat``.
    """
    snps = pd.read_table(snp_info_file, sep="\t")
    return (
        snps,
        read_allele_mat(tot_file, mat_dtype),
        read_allele_mat(a_file, mat_dtype),
        read_allele_mat(b_file, mat_dtype),
    )


def read_barcodes(bc_file: str):
    """Read a barcode file, one barcode per line.

    Args:
        bc_file: Path to the barcode list.

    Returns:
        List of barcodes, in file order.
    """
    barcodes = (
        pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
    )
    return barcodes


def read_barcodes_by_dataset(bc_file: str):
    """Read ``barcodes.tsv.gz`` and split each key back into its three fields.

    The parse is positional: no assay_type holds ``_`` and no raw barcode does either
    (asserted at write), so the last and first ``_`` bound the dataset_id.

    Args:
        bc_file: One ``{raw}_{dataset_id}_{assay_type}`` per line, no header, in
            matrix-column order.

    Returns:
        DataFrame with ``raw``, ``dataset_id``, ``assay_type`` and ``BARCODE`` (the full
        key), one row per matrix column, in file order.
    """
    full = pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0]
    head, assay_type = _split_once(full, bc_file, from_right=True)
    raw, dataset_id = _split_once(head, bc_file, from_right=False)
    return pd.DataFrame(
        {
            "raw": raw,
            "dataset_id": dataset_id,
            "assay_type": assay_type,
            "BARCODE": full,
        }
    )


def _split_once(values: pd.Series, bc_file: str, from_right: bool):
    """Split every barcode field on one ``_``, asserting each side is non-empty."""
    parts = (
        values.str.rsplit("_", n=1, expand=True)
        if from_right
        else values.str.split("_", n=1, expand=True)
    )
    assert parts.shape[1] == 2 and parts.notna().all().all(), (
        f"{bc_file}, every barcode must be '{{raw}}_{{dataset_id}}_{{assay_type}}'"
    )
    return parts[0], parts[1]


def read_chunks_from_atac_fragments(frag_file: str, chunksize=5_000_000):
    """Read a 10x ATAC fragment file in chunks, keeping its first four columns.

    Contig names are left as the file spells them: one sample runs to hundreds of
    millions of records, so a caller chr-normalizes the rows it keeps.

    Args:
        frag_file: Path to the (optionally gzipped) fragment TSV.
        chunksize: Records per chunk.

    Returns:
        Iterator of DataFrames with ``#CHR``, ``start``, ``end``, ``BC``.

    Notes/References:
        Format: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/fragments-file
    """
    return pd.read_csv(
        frag_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["#CHR", "start", "end", "BC"],
        dtype={0: str, 1: np.int64, 2: np.int64, 3: str},
        chunksize=chunksize,
    )


def uniquify_var_names(adata, label: str):
    """Keep the raw 10x gene symbols, then suffix the repeated ones.

    A 10x reference repeats a gene symbol across loci while its gene id stays unique, so
    ``var_names`` (the symbols) are not unique on read. ``var_names_make_unique`` appends
    ``-1``, ``-2`` to the second and later copies; the raw symbol is preserved in
    ``var["gene_symbol"]`` so a caller matching on symbols (the gene blacklist) still sees
    every copy. The GTF join keys on the gene id and is unaffected.

    Args:
        adata: AnnData from a 10x reader, var_names being gene symbols.
        label: Dataset label for the log line.

    Returns:
        The same AnnData, var_names unique.

    Notes/References:
        Suffixing rule (the first occurrence keeps the bare symbol):
        https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.var_names_make_unique.html
    """
    adata.var["gene_symbol"] = adata.var_names.astype(str)
    num_dup = int(adata.var_names.duplicated().sum())
    if num_dup > 0:
        logging.warning(
            f"{label}, #gene symbols repeated={num_dup}/{adata.n_vars}, "
            "suffixed to make var_names unique"
        )
    adata.var_names_make_unique()
    return adata


def read_10x_ranger_scRNA(matrix_h5):
    """Read one Cell Ranger gene-expression matrix into an AnnData.

    Args:
        matrix_h5: Path to ``filtered_feature_bc_matrix.h5``.

    Returns:
        AnnData of the gene-expression features only, with unique var_names.

    Notes/References:
        Format: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
    """
    import scanpy as sc

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r".*names are not unique", UserWarning)
        adata = sc.read_10x_h5(matrix_h5, gex_only=True)
    return uniquify_var_names(adata, os.path.basename(matrix_h5))


def read_10x_ranger_spatial(
    matrix_h5, names, paths, library_id, assay_type, load_images=True
):
    """Read one Space Ranger spatial dataset into an AnnData.

    squidpy takes a directory, so the Space Ranger layout is rebuilt as symlinks in a
    temporary directory for the read; names come from ``RANGER_*`` in const.py.

    Args:
        matrix_h5: Path to this dataset's feature-barcode matrix.
        names: Space Ranger filenames under spatial/, for this dataset.
        paths: Paths supplying those files, in the same order.
        library_id: Library id squidpy records in ``uns``; the dataset_id.
        assay_type: VISIUM | VISIUM3prime.
        load_images: Read the tissue images; VISIUM3prime must pass False.

    Returns:
        AnnData with unique var_names.

    Raises:
        AssertionError: load_images is set for VISIUM3prime.

    Notes/References:
        spatial/ layout: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs
    """
    import squidpy as sq

    if assay_type == "VISIUM3prime":
        assert not load_images, "VISIUM3prime, squidpy cannot load its tissue images"

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.symlink(
            os.path.abspath(matrix_h5),
            os.path.join(tmp_dir, RANGER_MATRIX_H5[0]),
        )
        spatial_dir = os.path.join(tmp_dir, RANGER_SPATIAL_DIR)
        os.makedirs(spatial_dir)
        for name, path in zip(names, paths):
            os.symlink(os.path.abspath(path), os.path.join(spatial_dir, name))
        logging.info(f"staged {len(names) + 1} files for squidpy: {names}")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r".*names are not unique", UserWarning)
            adata = sq.read.visium(
                tmp_dir, load_images=load_images, library_id=library_id
            )
    return uniquify_var_names(adata, library_id)


def write_snp_info(
    snps: pd.DataFrame,
    out_file: str,
):
    """Write the SNP feature axis of the allele matrices.

    ``PS`` (the phaser's phase-set label, which becomes the binning phase clusters) and
    ``seg_id`` are carried only when present.

    Args:
        snps: Filtered SNPs, in matrix-feature order.
        out_file: Output TSV path.

    Returns:
        The DataFrame written.
    """
    snp_cols = ["#CHR", "POS", "POS0", "START", "END", "GT", "PHASE"]
    # upstream phaser's phaseset label
    if "PS" in snps.columns:
        snp_cols.append("PS")
    logging.info(f"phase set (PS) column carried: {'PS' in snps.columns}")

    snp_cols += ["region_id"]
    if "seg_id" in snps.columns:
        snp_cols.append("seg_id")
    snp_cols += ["feature_id", "feature_type"]
    snp_info = snps[snp_cols]
    snp_info.to_csv(out_file, sep="\t", header=True, index=False)
    return snp_info


def write_bb_file(bbs: pd.DataFrame, out_file: str):
    """Write ``bb.tsv.gz``, the feature axis of every bb matrix.

    The three coordinate columns are required and the optional ones are written when
    present, so one schema covers all three modes; the binning internals (``bb_id``,
    ``BLOCKSIZE``, ``seg_id``, ``PS``, the cluster keys) are dropped.

    Args:
        bbs: bbs carrying at least ``#CHR``, ``START``, ``END``.
        out_file: Output TSV path; ``.gz`` is compressed by pandas.

    Returns:
        The DataFrame written.
    """
    bb_cols = ["#CHR", "START", "END"]
    missing = [c for c in bb_cols if c not in bbs.columns]
    assert not missing, f"bb table, missing column(s) {missing}"
    bb_cols += [
        c
        for c in (
            "#SNPS",
            "region_id",
            "switchprobs",
            "feature_id",
            "#feature",
            "is_loh",
        )
        if c in bbs.columns
    ]
    bb_out = bbs[bb_cols]
    bb_out.to_csv(out_file, sep="\t", header=True, index=False)
    return bb_out


def write_sample_ids(
    sample_id: str,
    dataset_ids: list,
    sample_types: list,
    assay_types: list,
    out_file: str,
    rdr_base_dataset_ids=None,
):
    """Write ``sample_ids.tsv``, one row per observation of the bb matrices.

    Every column is a sample-file record key spelled the same way, except the derived
    leading ``SAMPLE``.

    Args:
        sample_id: Sample (patient) id, one per file.
        dataset_ids: Dataset id per observation, in matrix-observation order.
        sample_types: ``tumor``/``normal`` per observation.
        assay_types: Assay type per observation; the column is omitted when None.
        out_file: Output TSV path.
        rdr_base_dataset_ids: RDR baseline dataset id per observation; the column is
            omitted when None.

    Returns:
        The DataFrame written.
    """
    record = {
        "sample_id": sample_id,
        "dataset_id": dataset_ids,
        "sample_type": sample_types,
        "assay_type": assay_types,
        "rdr_base_dataset_id": rdr_base_dataset_ids,
    }
    # a multiome pair shares one dataset_id, so the assay is what keeps SAMPLE unique
    repeated = len(set(dataset_ids)) != len(dataset_ids)
    samples = [
        f"{sample_id}_{d}_{a}" if repeated else f"{sample_id}_{d}"
        for d, a in zip(dataset_ids, assay_types)
    ]
    assert len(set(samples)) == len(samples), (
        f"sample_ids.tsv, duplicate SAMPLE: {sorted(samples)}"
    )
    sample_dict = {"SAMPLE": samples}
    sample_dict.update({key: val for key, val in record.items() if val is not None})
    sample_df = pd.DataFrame(sample_dict)
    sample_df.to_csv(out_file, sep="\t", header=True, index=False)
    return sample_df
