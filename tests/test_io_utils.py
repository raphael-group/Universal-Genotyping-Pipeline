#!/usr/bin/env python3
"""Unit tests for the chromosome-naming helpers and the size-file reader.

Last update: 2026-08-06

Covers:
- naming: strip and add the chr prefix, per row
- sizes: file order is kept, both namings are read
- bundled: every shipped size file covers its chromosomes
"""

import os
import sys

import pytest

_REPO = os.path.join(os.path.dirname(__file__), "..")
for _sub in ("config", "workflow/scripts/script_utils"):
    sys.path.insert(0, os.path.join(_REPO, _sub))
io_utils = pytest.importorskip("io_utils")
utils = pytest.importorskip("utils")

HG38 = os.path.join(_REPO, "resources", "data", "hg38.chrom.sizes")
MM10 = os.path.join(_REPO, "resources", "data", "mm10.chrom.sizes")


@pytest.mark.parametrize(
    "name,core",
    [("chr1", "1"), ("1", "1"), ("chrX", "X"), ("X", "X"), ("CHR22", "22"), (22, "22")],
)
def test_strip_chr_prefix(name, core):
    """Either spelling, any case, str or int, reduces to the same core name."""
    assert utils.strip_chr_prefix(name) == core


def test_get_chr_sizes_keeps_file_names_and_order(tmp_path):
    """Names come back as the file spells them, in file order."""
    path = tmp_path / "genome_size.txt"
    path.write_text("chr2\t20\nchr1\t10\nchrX\t5\n")
    sizes = io_utils.read_chrom_sizes(str(path))
    assert list(sizes) == ["chr2", "chr1", "chrX"]
    assert sizes["chr1"] == 10


def test_bare_contig_size_file(tmp_path):
    """An Ensembl-style file has no chr prefix; the cores still match."""
    path = tmp_path / "genome_size.txt"
    path.write_text("1\t10\n2\t20\nX\t5\n")
    sizes = io_utils.read_chrom_sizes(str(path))
    assert list(sizes) == ["1", "2", "X"]
    assert [utils.strip_chr_prefix(c) for c in sizes] == ["1", "2", "X"]


def test_bundled_size_files_cover_their_chromosomes():
    """The shipped size files hold the chromosome set each genome has."""
    hg38 = io_utils.read_chrom_sizes(HG38)
    mm10 = io_utils.read_chrom_sizes(MM10)
    assert {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]} <= set(hg38)
    assert {f"chr{c}" for c in list(range(1, 20)) + ["X", "Y"]} <= set(mm10)
    assert "chr20" not in mm10


@pytest.mark.parametrize(
    "names,want",
    [
        (["22", "X"], ["chr22", "chrX"]),
        (["chr22", "chrX"], ["chr22", "chrX"]),
        ([], []),
    ],
)
def test_add_chr_prefix(names, want):
    """Readers that take a file in the genome's naming normalize on ingest."""
    pd = pytest.importorskip("pandas")
    assert utils.add_chr_prefix(pd.Series(names, dtype=str)).tolist() == want


def test_add_chr_prefix_is_per_row():
    """A GTF can start on a scaffold, so the decision cannot be made per file."""
    pd = pytest.importorskip("pandas")
    got = utils.add_chr_prefix(pd.Series(["GL000009.2", "chr1", "2"]))
    assert got.tolist() == ["chrGL000009.2", "chr1", "chr2"]


##################################################
# RNA AnnData concatenation


def _var(gene_ids, symbols, genome="ref"):
    """A per-dataset ``var`` as the 10x readers leave it: gene-id index, symbols kept."""
    pd = pytest.importorskip("pandas")
    return pd.DataFrame(
        {"gene_ids": gene_ids, "gene_symbol": symbols, "genome": genome},
        index=pd.Index(gene_ids, name=None),
    )


def test_combine_var_frames_keeps_columns_across_disjoint_references():
    """A gene missing from one reference must not delete the GTF join key."""
    pd = pytest.importorskip("pandas")
    frames = {
        "d0": _var(["g1", "g2"], ["A", "B"]),
        "d1": _var(["g2", "g3"], ["B", "C"]),
    }
    got = io_utils.combine_var_frames(frames, pd.Index(["g1", "g2", "g3"]))
    assert list(got.columns) == ["gene_ids", "gene_symbol", "genome"]
    assert got["gene_ids"].tolist() == ["g1", "g2", "g3"]
    assert got["gene_symbol"].tolist() == ["A", "B", "C"]


def test_combine_var_frames_first_dataset_wins_on_conflict(caplog):
    """A symbol disagreement is reported, not resolved; the id stays the key."""
    pd = pytest.importorskip("pandas")
    frames = {
        "d0": _var(["g1"], ["A"]),
        "d1": _var(["g1"], ["A_alias"]),
    }
    with caplog.at_level("WARNING"):
        got = io_utils.combine_var_frames(frames, pd.Index(["g1"]))
    assert got["gene_symbol"].tolist() == ["A"]
    assert "symbol differs between datasets" in caplog.text


def test_combine_var_frames_reindexes_to_the_given_order():
    """The output row order is the concatenated object's, not the frames'."""
    pd = pytest.importorskip("pandas")
    frames = {"d0": _var(["g2", "g1"], ["B", "A"])}
    got = io_utils.combine_var_frames(frames, pd.Index(["g1", "g2"]))
    assert got["gene_symbol"].tolist() == ["A", "B"]


def test_concat_rna_adatas_joins_on_gene_id_not_symbol():
    """Suffixing is per file, so symbols cannot key the join; gene ids can."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    # the same locus is the first 'DUP' copy in d0 and the second in d1
    a0 = anndata.AnnData(np.ones((2, 2)), var=_var(["g1", "g2"], ["DUP", "X"]))
    a1 = anndata.AnnData(np.ones((3, 2)), var=_var(["g2", "g1"], ["X", "DUP"]))
    a0.obs_names = ["c0_d0", "c1_d0"]
    a1.obs_names = ["c0_d1", "c1_d1", "c2_d1"]
    a0.var_names = ["DUP", "X"]
    a1.var_names = ["X", "DUP"]
    got = io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")
    assert got.n_obs == 5
    assert got.n_vars == 2, "a shared gene must not split on its symbol suffix"
    assert set(got.var["gene_ids"]) == {"g1", "g2"}
    assert got.var_names.is_unique


def test_concat_rna_adatas_rejects_repeated_gene_ids():
    """A probe-barcode matrix repeats the gene id, so it cannot key the join."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    var = _var(["g1", "g1"], ["A", "A"])
    var.index = ["A", "A-1"]
    a0 = anndata.AnnData(np.ones((2, 2)), var=var)
    a1 = anndata.AnnData(np.ones((2, 1)), var=_var(["g2"], ["B"]))
    a1.var_names = ["B"]
    with pytest.raises(AssertionError, match="probe-barcode"):
        io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")


def test_concat_rna_adatas_requires_the_id_column():
    """A missing id column fails naming the dataset, not with a bare KeyError."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    a0 = anndata.AnnData(np.ones((2, 1)), var=pd.DataFrame(index=["A"]))
    a1 = anndata.AnnData(np.ones((2, 1)), var=_var(["g2"], ["B"]))
    with pytest.raises(AssertionError, match="d0, var has no 'gene_ids'"):
        io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")


def test_concat_rna_adatas_keeps_the_union_and_reports_the_gap(caplog):
    """A gene one reference lacks is kept, zero-filled there, and reported."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    a0 = anndata.AnnData(np.ones((2, 2)), var=_var(["g1", "g2"], ["A", "B"]))
    a1 = anndata.AnnData(np.ones((2, 3)), var=_var(["g1", "g2", "g3"], ["A", "B", "C"]))
    a0.obs_names = ["c0_d0", "c1_d0"]
    a1.obs_names = ["c0_d1", "c1_d1"]
    a0.var_names = ["A", "B"]
    a1.var_names = ["A", "B", "C"]
    with caplog.at_level("WARNING"):
        got = io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")
    assert set(got.var["gene_ids"]) == {"g1", "g2", "g3"}
    assert got.n_obs == 4
    assert got[["c0_d0", "c1_d0"], got.var["gene_ids"] == "g3"].X.sum() == 0
    assert "d1, #genes absent from another dataset's reference=1/3" in caplog.text


def test_concat_rna_adatas_keeps_everything_when_references_agree(caplog):
    """Identical gene sets make the intersection the whole set and log nothing."""
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    a0 = anndata.AnnData(np.ones((2, 2)), var=_var(["g1", "g2"], ["A", "B"]))
    a1 = anndata.AnnData(np.ones((2, 2)), var=_var(["g2", "g1"], ["B", "A"]))
    a0.obs_names = ["c0_d0", "c1_d0"]
    a1.obs_names = ["c0_d1", "c1_d1"]
    a0.var_names = ["A", "B"]
    a1.var_names = ["B", "A"]
    with caplog.at_level("WARNING"):
        got = io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")
    assert got.n_vars == 2
    assert "absent from another dataset's reference" not in caplog.text


def test_concat_rna_adatas_result_writes_to_h5ad(tmp_path):
    """The result must survive write_h5ad: repeated symbols are the failing case.

    ``write_h5ad`` rejects a var index whose name is also a column holding different
    values, which is what assigning the named ``gene_symbol`` Series to ``var_names``
    produces once ``var_names_make_unique`` suffixes the repeats.
    """
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    a0 = anndata.AnnData(np.ones((2, 2)), var=_var(["g1", "g2"], ["DUP", "DUP"]))
    a1 = anndata.AnnData(
        np.ones((3, 3)), var=_var(["g1", "g2", "g3"], ["DUP", "DUP", "C"])
    )
    a0.obs_names = ["c0_d0", "c1_d0"]
    a1.obs_names = ["c0_d1", "c1_d1", "c2_d1"]
    a0.var_names = ["DUP", "DUP-1"]
    a1.var_names = ["DUP", "DUP-1", "C"]

    got = io_utils.concat_rna_adatas({"d0": a0, "d1": a1}, "gene_ids")
    assert got.var.index.name is None
    path = tmp_path / "concat.h5ad"
    got.write_h5ad(path, compression="gzip")
    back = anndata.read_h5ad(path)
    assert back.var_names.is_unique
    assert list(back.var["gene_ids"]) == ["g1", "g2", "g3"]


def test_read_barcodes_by_dataset_keeps_underscored_raw_barcodes(tmp_path):
    """A Visium HD barcode carries its own "_", so only the last two may be split on."""
    path = tmp_path / "barcodes.tsv"
    path.write_text(
        "s_002um_00000_00000-1_D1_VISIUM\n"
        "s_008um_00123_00456-1_D1_VISIUM\n"
        "AAACCTGAGAAACCAT-1_D2_scATAC\n"
    )
    cells = io_utils.read_barcodes_by_dataset(str(path))
    assert cells["raw"].tolist() == [
        "s_002um_00000_00000-1",
        "s_008um_00123_00456-1",
        "AAACCTGAGAAACCAT-1",
    ]
    assert cells["dataset_id"].tolist() == ["D1", "D1", "D2"]
    assert cells["assay_type"].tolist() == ["VISIUM", "VISIUM", "scATAC"]
    assert cells["BARCODE"].tolist() == path.read_text().split()


def test_read_barcodes_by_dataset_rejects_a_key_missing_a_field(tmp_path):
    """Two fields cannot be recovered from one "_"."""
    path = tmp_path / "barcodes.tsv"
    path.write_text("AAACCTGAGAAACCAT-1_scATAC\n")
    with pytest.raises(AssertionError, match="every barcode must be"):
        io_utils.read_barcodes_by_dataset(str(path))
