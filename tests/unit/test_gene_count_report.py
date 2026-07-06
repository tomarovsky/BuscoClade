"""Unit tests for workflow/scripts/gene_count_report.py.

Covers per-species single-copy vs common-set accounting and the source
classification used to surface gene loss / bottleneck species.
"""
import gene_count_report as gcr


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_classify_source():
    assert gcr.classify_source("Sample1.ref1.AltRef") == "reconstruct:vcf"
    assert gcr.classify_source("Sample1.ref1.Consensus") == "reconstruct:consensus"
    assert gcr.classify_source("Homo_sapiens") == "genome"


def test_read_ids_strips_and_skips_blank_lines(tmp_path):
    f = tmp_path / "x.ids"
    f.write_text("gene1\n\n  gene2  \ngene3\n\n")
    assert gcr.read_ids(f) == {"gene1", "gene2", "gene3"}


def test_read_ids_missing_file_is_empty_set(tmp_path):
    assert gcr.read_ids(tmp_path / "does_not_exist.ids") == set()


# ---------------------------------------------------------------------------
# End-to-end report
# ---------------------------------------------------------------------------

def test_main_writes_sorted_report(tmp_path, monkeypatch):
    sc = tmp_path / "single_copy"
    mc = tmp_path / "multi_copy"
    sc.mkdir()
    mc.mkdir()

    (sc / "speciesA.ids").write_text("g1\ng2\ng3\ng4\n")  # 4 single-copy
    (sc / "speciesB.ids").write_text("g1\ng2\n")           # 2 single-copy
    (mc / "speciesA.ids").write_text("g5\n")               # 1 multi-copy
    # speciesB has no multi-copy file -> counts as 0

    common = tmp_path / "common.ids"
    common.write_text("g1\ng2\ng3\n")  # 3 common genes

    out = tmp_path / "gene_counts.tsv"
    monkeypatch.setattr(
        "sys.argv",
        [
            "gene_count_report.py",
            "--single_copy_dir", str(sc),
            "--multi_copy_dir", str(mc),
            "--common", str(common),
            "--output", str(out),
        ],
    )
    gcr.main()

    lines = out.read_text().splitlines()
    assert lines[0] == "# common_genes\t3"
    assert lines[1] == "species\tsource\tsingle_copy\tmulti_copy\tin_common\tsingle_copy_not_in_common"
    # Sorted by single_copy ascending: speciesB (2) before speciesA (4).
    assert lines[2] == "speciesB\tgenome\t2\t0\t2\t0"
    assert lines[3] == "speciesA\tgenome\t4\t1\t3\t1"
