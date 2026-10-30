#!/usr/bin/env python
"""Build UGP inputs for one Visium HD 3' slide at any even bin size from Space Ranger outs/.

Runpeng Luo (2026-10-29)

Dependencies:
    numpy, pandas, scipy, h5py, pyarrow, pysam

Usage:
    python resources/scripts/prepare_visium_hd.py --ranger-outs <outs> --bin-um <X> --out-dir <dir>
                                --sample-id <id> [--dataset-id U1] [--sample-type tumor]
                                [--reference-version GRCh38] [--from-2um] [--skip-bam]
                                [--contigs chr1 ...] [--threads 16] [-v]
            ranger-outs  # Space Ranger count outs/ (default layout, --create-bam=true)
            bin-um       # target square bin edge in um; even, >= 2
            out-dir      # output dir for this bin size
            sample-id    # UGP sample_id of the slide
            from-2um     # build 8/16um from the 2um grid instead of copying Space Ranger's
            skip-bam     # do not write the retagged BAM
            contigs      # retag only these contigs (default: every contig in the header)

Inputs
    outs/possorted_genome_bam.bam{,.bai}: BAM with CB = square_002um barcode
    outs/barcode_mappings.parquet: square_002um -> square_008um/square_016um/cell
    outs/binned_outputs/square_002um/{raw,filtered}_feature_bc_matrix.h5: 2um counts
    outs/binned_outputs/square_002um/spatial/{tissue_positions.parquet,scalefactors_json.json}
    outs/binned_outputs/square_{008,016}um/: Space Ranger's own 8/16um outputs
    outs/spatial/*.{png,jpg,tiff}: slide images
Parameters
    bin_um: X; a bin holds the k x k 2um squares (r // k, c // k), k = X / 2
Outputs:
    <out_dir>/possorted_genome_bam.bam{,.bai}: CB rewritten to s_<XXX>um_<row>_<col>-1 (X > 2)
    <out_dir>/filtered_feature_bc_matrix.h5: 10x h5, genes x in-tissue bins
    <out_dir>/filtered_feature_bc_matrix/barcodes.tsv.gz: in-tissue bins, h5 order
    <out_dir>/spatial/tissue_positions.parquet: Space Ranger schema, full bin grid
    <out_dir>/spatial/scalefactors_json.json, spatial/*.{png,jpg,tiff}
    <out_dir>/barcode_mappings.parquet: Space Ranger mapping plus square_<XXX>um
    <out_dir>/sample_sheet.json: UGP sample sheet with one VISIUM3prime record
Notes/References:
    Space Ranger nests its 8um and 16um bins in the 2um grid as (r // k, c // k); the script
    asserts this against barcode_mappings.parquet whenever X is 8 or 16, so other sizes follow
    the same rule. For X = 2 no BAM is written and the sheet points at the Space Ranger BAM.
    For X = 8/16 without --from-2um, the matrix, barcodes and spatial/ are Space Ranger's files.
    Otherwise:
    - counts of a bin are the sum over its 2um squares in square_002um/raw_feature_bc_matrix.h5
      (Space Ranger default output); without it, the filtered 2um matrix is summed, which
      omits out-of-tissue 2um squares inside in-tissue bins;
    - in_tissue of a bin is 1 if any member 2um bin is in tissue; the h5 holds those bins;
    - pixel centres are the mean of the member 2um centres, shifted to the full-bin centre by
      an affine (array_row, array_col) -> pixel fit on the 2um grid (moves edge bins only);
    - scalefactors are the 2um ones with bin_size_um = X and
      spot_diameter_fullres = X / microns_per_pixel.
    BAM: CB is mapped from its own 2um row/col, so no lookup table is held. Reads without CB, or
    with a CB outside the 2um grid, are dropped, as are unplaced reads (contigs are fetched one
    by one in parallel, then concatenated in header order).
    Binned outputs: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview
    barcode_mappings.parquet: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/segmented-outputs
    10x h5: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
    UGP sample sheet: docs/sample_sheet.md; walkthrough: docs/spatial_genotyping.md
"""

import argparse
import gzip
import json
import logging
import multiprocessing as mp
import os
import shutil
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import csc_matrix, csr_matrix

NATIVE_UM = (2, 8, 16)
POS_COLS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
log = logging.getLogger(__name__)
IMAGE_EXT = (".png", ".jpg", ".tiff", ".tif")


def parse_args():
    """Parse command-line arguments.

    Returns:
        argparse.Namespace with the options in the module docstring.
    """
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--ranger-outs", required=True)
    p.add_argument("--bin-um", type=int, required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--dataset-id", default="U1")
    p.add_argument("--sample-type", default="tumor", choices=["tumor", "normal"])
    p.add_argument("--reference-version", default="GRCh38")
    p.add_argument("--from-2um", action="store_true")
    p.add_argument("--skip-bam", action="store_true")
    p.add_argument("--contigs", nargs="+", default=None)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args()
    if a.bin_um < 2 or a.bin_um % 2:
        p.error(f"--bin-um must be even and >= 2, got {a.bin_um}")
    return a


def bin_name(x):
    """Return the Space Ranger bin label, e.g. square_016um."""
    return f"square_{x:03d}um"


def bin_barcodes(x, rows, cols):
    """Format Space Ranger barcodes s_<XXX>um_<row:05d>_<col:05d>-1.

    Args:
        x: Bin size in um.
        rows: Integer array of bin rows.
        cols: Integer array of bin columns.

    Returns:
        numpy array of str.
    """
    r = pd.Series(rows).astype(str).str.zfill(5)
    c = pd.Series(cols).astype(str).str.zfill(5)
    return (f"s_{x:03d}um_" + r + "_" + c + "-1").to_numpy()


def load_2um_grid(outs):
    """Read the 2um tissue positions.

    Args:
        outs: Space Ranger outs/ dir.

    Returns:
        DataFrame with POS_COLS, one row per 2um square.
    """
    path = os.path.join(
        outs, "binned_outputs", bin_name(2), "spatial", "tissue_positions.parquet"
    )
    pos = pd.read_parquet(path, columns=POS_COLS)
    log.info(f"{path}: {len(pos)} 2um squares, {int(pos.in_tissue.sum())} in tissue")
    return pos


def write_mapping(outs, out_dir, x, pos2, k):
    """Write barcode_mappings.parquet with a square_<XXX>um column; check 8/16um nesting.

    Args:
        outs: Space Ranger outs/ dir.
        out_dir: Output dir.
        x: Bin size in um.
        pos2: 2um positions from load_2um_grid.
        k: 2um squares per bin edge.

    Raises:
        AssertionError: the (r // k, c // k) rule disagrees with Space Ranger's column.
    """
    bm = pd.read_parquet(os.path.join(outs, "barcode_mappings.parquet"))
    rc = pos2.set_index("barcode").loc[bm["square_002um"], ["array_row", "array_col"]]
    tgt = bin_barcodes(x, rc.array_row.to_numpy() // k, rc.array_col.to_numpy() // k)
    col = bin_name(x)
    if col in bm.columns:
        n_bad = int((bm[col].astype(str).to_numpy() != tgt).sum())
        assert n_bad == 0, f"{n_bad} 2um squares disagree with Space Ranger's {col}"
        log.info(
            f"{col}: (r // {k}, c // {k}) matches Space Ranger on {len(bm)} squares"
        )
    else:
        bm[col] = tgt
    bm.to_parquet(os.path.join(out_dir, "barcode_mappings.parquet"), index=False)


def build_positions(pos2, x, k):
    """Aggregate 2um positions into bins.

    Args:
        pos2: 2um positions from load_2um_grid.
        x: Bin size in um.
        k: 2um squares per bin edge.

    Returns:
        (positions DataFrame over the full bin grid sorted by row/col, Series mapping each
        2um barcode to its bin barcode).
    """
    r2 = pos2.array_row.to_numpy()
    c2 = pos2.array_col.to_numpy()
    rows, cols = r2 // k, c2 // k
    g = (
        pd.DataFrame(
            {
                "array_row": rows,
                "array_col": cols,
                "in_tissue": pos2.in_tissue.to_numpy(),
                "r2": r2,
                "c2": c2,
                "pxl_row_in_fullres": pos2.pxl_row_in_fullres.to_numpy(),
                "pxl_col_in_fullres": pos2.pxl_col_in_fullres.to_numpy(),
            }
        )
        .groupby(["array_row", "array_col"], sort=True)
        .agg(
            in_tissue=("in_tissue", "max"),
            r2=("r2", "mean"),
            c2=("c2", "mean"),
            pxl_row_in_fullres=("pxl_row_in_fullres", "mean"),
            pxl_col_in_fullres=("pxl_col_in_fullres", "mean"),
        )
        .reset_index()
    )
    design = np.column_stack([np.ones(len(r2)), r2, c2])
    pxl = pos2[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
    coef = np.linalg.lstsq(design, pxl, rcond=None)[0][1:]
    shift = np.column_stack(
        [
            g.array_row * k + (k - 1) / 2 - g.r2,
            g.array_col * k + (k - 1) / 2 - g.c2,
        ]
    )
    g[["pxl_row_in_fullres", "pxl_col_in_fullres"]] += shift @ coef
    g["barcode"] = bin_barcodes(x, g.array_row.to_numpy(), g.array_col.to_numpy())
    g["in_tissue"] = g.in_tissue.astype(np.uint8)
    g["array_row"] = g.array_row.astype(np.uint32)
    g["array_col"] = g.array_col.astype(np.uint32)
    to_bin = pd.Series(bin_barcodes(x, rows, cols), index=pos2.barcode.to_numpy())
    log.info(f"{bin_name(x)}: {len(g)} bins, {int(g.in_tissue.sum())} in tissue")
    return g[POS_COLS], to_bin


def write_scalefactors(outs, dst, x):
    """Write the 2um scalefactors with bin_size_um and spot_diameter_fullres set for x."""
    src = os.path.join(
        outs, "binned_outputs", bin_name(2), "spatial", "scalefactors_json.json"
    )
    with open(src) as f:
        sf = json.load(f)
    sf["bin_size_um"] = float(x)
    sf["spot_diameter_fullres"] = x / sf["microns_per_pixel"]
    with open(dst, "w") as f:
        json.dump(sf, f, indent=2)


def copy_images(outs, sp_dir):
    """Copy the slide images of outs/spatial/ into sp_dir."""
    src = os.path.join(outs, "spatial")
    names = [n for n in sorted(os.listdir(src)) if n.endswith(IMAGE_EXT)]
    for n in names:
        shutil.copy2(os.path.join(src, n), sp_dir)
    log.info(f"copied {len(names)} images from {src}")


def read_10x_h5(path):
    """Read a 10x h5 matrix.

    Args:
        path: filtered_feature_bc_matrix.h5.

    Returns:
        (genes x barcodes csc_matrix, barcodes as str array).
    """
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        X = csc_matrix(
            (m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])
        )
        bc = m["barcodes"][:].astype(str)
    log.info(f"{path}: {X.shape[0]} genes x {X.shape[1]} barcodes, nnz {X.nnz}")
    return X, bc


def write_10x_h5(src, dst, X, barcodes):
    """Write a 10x h5 with features and file attributes copied from src.

    Args:
        src: 10x h5 whose matrix/features group and root attrs are copied.
        dst: Output path.
        X: genes x barcodes count matrix.
        barcodes: Column barcodes.
    """
    X = csc_matrix(X)
    X.sort_indices()
    with h5py.File(src, "r") as fs, h5py.File(dst, "w") as fd:
        for key, val in fs.attrs.items():
            fd.attrs[key] = val
        m = fd.create_group("matrix")
        m.create_dataset("barcodes", data=np.asarray(barcodes, dtype="S"))
        m.create_dataset("data", data=X.data.astype(np.int32), compression="gzip")
        m.create_dataset("indices", data=X.indices.astype(np.int64), compression="gzip")
        m.create_dataset("indptr", data=X.indptr.astype(np.int64), compression="gzip")
        m.create_dataset("shape", data=np.array(X.shape, dtype=np.int32))
        fs.copy(fs["matrix/features"], m, name="features")


def write_barcodes(path, barcodes):
    """Write barcodes.tsv.gz, one barcode per line."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write("\n".join(barcodes) + "\n")


def build_from_2um(outs, out_dir, x, pos2, k):
    """Write the matrix, barcodes and spatial/ of bin size x from the 2um outputs.

    Args:
        outs: Space Ranger outs/ dir.
        out_dir: Output dir.
        x: Bin size in um.
        pos2: 2um positions from load_2um_grid.
        k: 2um squares per bin edge.
    """
    sp_dir = os.path.join(out_dir, "spatial")
    os.makedirs(sp_dir, exist_ok=True)
    pos, to_bin = build_positions(pos2, x, k)
    pos.to_parquet(os.path.join(sp_dir, "tissue_positions.parquet"), index=False)
    write_scalefactors(outs, os.path.join(sp_dir, "scalefactors_json.json"), x)
    copy_images(outs, sp_dir)

    d2 = os.path.join(outs, "binned_outputs", bin_name(2))
    h5_2um = os.path.join(d2, "raw_feature_bc_matrix.h5")
    if not os.path.exists(h5_2um):
        h5_2um = os.path.join(d2, "filtered_feature_bc_matrix.h5")
        log.warning(
            f"no 2um raw_feature_bc_matrix.h5; using {h5_2um}, which omits UMIs of "
            "out-of-tissue 2um squares inside in-tissue bins"
        )
    X2, bc2 = read_10x_h5(h5_2um)
    keep = pos.barcode[pos.in_tissue == 1].to_numpy()
    col = pd.Index(keep).get_indexer(to_bin.loc[bc2].to_numpy())
    used = np.flatnonzero(col >= 0)
    A = csr_matrix(
        (np.ones(len(used), dtype=np.int64), (used, col[used])),
        shape=(len(bc2), len(keep)),
    )
    X = X2 @ A
    log.info(f"{bin_name(x)}: {X.shape[0]} genes x {X.shape[1]} bins, nnz {X.nnz}")
    write_10x_h5(
        h5_2um, os.path.join(out_dir, "filtered_feature_bc_matrix.h5"), X, keep
    )
    write_barcodes(
        os.path.join(out_dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), keep
    )


def copy_native(outs, out_dir, x):
    """Copy Space Ranger's own binned outputs of x (8/16um) and the slide images."""
    src = os.path.join(outs, "binned_outputs", bin_name(x))
    shutil.copy2(os.path.join(src, "filtered_feature_bc_matrix.h5"), out_dir)
    shutil.copytree(
        os.path.join(src, "filtered_feature_bc_matrix"),
        os.path.join(out_dir, "filtered_feature_bc_matrix"),
        dirs_exist_ok=True,
    )
    sp_dir = os.path.join(out_dir, "spatial")
    shutil.copytree(os.path.join(src, "spatial"), sp_dir, dirs_exist_ok=True)
    copy_images(outs, sp_dir)
    log.info(f"copied Space Ranger {src}")


def _retag_contig(job):
    """Retag the CB of one contig's reads; worker of retag_bam.

    Args:
        job: (bam_in, part_out, contig, prefix, k).

    Returns:
        (contig, reads written, reads dropped).
    """
    bam_in, part, contig, prefix, k = job
    kept = dropped = 0
    with (
        pysam.AlignmentFile(bam_in, "rb") as fin,
        pysam.AlignmentFile(part, "wb", template=fin) as fout,
    ):
        for r in fin.fetch(contig):
            if not r.has_tag("CB"):
                dropped += 1
                continue
            cb = r.get_tag("CB")
            if not (cb.startswith("s_002um_") and len(cb) == 21):
                dropped += 1
                continue
            row, col = int(cb[8:13]), int(cb[14:19])
            r.set_tag("CB", f"{prefix}{row // k:05d}_{col // k:05d}-1", value_type="Z")
            fout.write(r)
            kept += 1
    return contig, kept, dropped


def retag_bam(bam_in, bam_out, x, k, threads, contigs=None):
    """Write bam_in with CB mapped from 2um to x um bins, then index it.

    Args:
        bam_in: Space Ranger position-sorted BAM with 2um CB tags.
        bam_out: Output BAM; written via a temporary name and renamed when indexed.
        x: Bin size in um.
        k: 2um squares per bin edge.
        threads: Worker processes, one contig each.
        contigs: Contigs to retag; default every contig in the header.
    """
    with pysam.AlignmentFile(bam_in, "rb") as f:
        refs = list(f.references)
    contigs = refs if contigs is None else [c for c in refs if c in set(contigs)]
    prefix = f"s_{x:03d}um_"
    tmp_dir = tempfile.mkdtemp(prefix=".retag_", dir=os.path.dirname(bam_out))
    jobs = [
        (bam_in, os.path.join(tmp_dir, f"{i:04d}.bam"), c, prefix, k)
        for i, c in enumerate(contigs)
    ]
    log.info(f"retag {bam_in} -> {bam_out}: {len(contigs)} contigs, {threads} workers")
    kept = dropped = 0
    with mp.Pool(threads) as pool:
        for contig, n_kept, n_drop in pool.imap_unordered(_retag_contig, jobs):
            kept += n_kept
            dropped += n_drop
            log.debug(f"{contig}: kept {n_kept}, dropped {n_drop}")
    tmp_bam = bam_out + ".tmp.bam"
    pysam.cat("-o", tmp_bam, *[j[1] for j in jobs])
    pysam.index("-@", str(threads), tmp_bam)
    os.replace(tmp_bam + ".bai", bam_out + ".bai")
    os.replace(tmp_bam, bam_out)
    shutil.rmtree(tmp_dir)
    log.info(f"wrote {bam_out}: kept {kept} reads, dropped {dropped} without a 2um CB")


def write_sample_sheet(a, out_dir, bam):
    """Write sample_sheet.json with one VISIUM3prime record pointing at out_dir.

    Args:
        a: Parsed arguments.
        out_dir: Output dir (absolute).
        bam: Alignment path for the record.
    """
    rec = {
        "sample_id": a.sample_id,
        "dataset_id": a.dataset_id,
        "assay_type": "VISIUM3prime",
        "sample_type": a.sample_type,
        "platform": f"10x Genomics Visium HD 3prime, {bin_name(a.bin_um)} bins",
        "reference_version": a.reference_version,
        "files": {
            "alignment": bam,
            "alignment_index": bam + ".bai",
            "barcodes": os.path.join(
                out_dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz"
            ),
            "matrix_h5": os.path.join(out_dir, "filtered_feature_bc_matrix.h5"),
            "tissue_positions": os.path.join(
                out_dir, "spatial", "tissue_positions.parquet"
            ),
            "scalefactors": os.path.join(out_dir, "spatial", "scalefactors_json.json"),
        },
        "meta": {"ranger_outs": os.path.abspath(a.ranger_outs)},
    }
    path = os.path.join(out_dir, "sample_sheet.json")
    with open(path, "w") as f:
        json.dump({"version": 1, "samples": [rec]}, f, indent=2)
        f.write("\n")
    log.info(f"wrote {path}")


def main():
    """Run the steps in the module docstring."""
    a = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if a.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )
    outs = os.path.abspath(a.ranger_outs)
    out_dir = os.path.abspath(a.out_dir)
    x, k = a.bin_um, a.bin_um // 2
    os.makedirs(out_dir, exist_ok=True)
    log.info(f"{outs} -> {out_dir}: {bin_name(x)}, k = {k}")

    pos2 = load_2um_grid(outs)
    write_mapping(outs, out_dir, x, pos2, k)
    if x in NATIVE_UM and not a.from_2um:
        copy_native(outs, out_dir, x)
    else:
        build_from_2um(outs, out_dir, x, pos2, k)

    bam_src = os.path.join(outs, "possorted_genome_bam.bam")
    if x == 2:
        bam = bam_src
    else:
        bam = os.path.join(out_dir, "possorted_genome_bam.bam")
        if a.skip_bam:
            log.info("skip BAM")
        elif os.path.exists(bam + ".bai"):
            log.info(f"skip BAM: {bam} exists")
        else:
            retag_bam(bam_src, bam, x, k, a.threads, a.contigs)
    write_sample_sheet(a, out_dir, bam)


if __name__ == "__main__":
    main()
