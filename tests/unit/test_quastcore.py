"""Unit tests for workflow/scripts/quastcore.py.

Covers the assembly-statistics primitives: transparent (de)compression, FASTA
parsing, per-threshold contig/length counts, N content, GC%, and N/L statistics.
"""
import gzip

import quastcore as qc


FASTA = ">c1\nATGCATGCAT\n>c2\nGGGCCCNNNN\n>c3\nAT\n"
# c1: len 10, GC 4, N 0
# c2: len 10, GC 6, N 4
# c3: len 2,  GC 0, N 0


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def test_metaopen_plain_and_gzip(tmp_path):
    plain = tmp_path / "seq.fasta"
    plain.write_text(FASTA)
    with qc.metaopen(str(plain), "rt") as fh:
        assert fh.read() == FASTA

    gz = tmp_path / "seq.fasta.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(FASTA)
    with qc.metaopen(str(gz), "rt") as fh:
        assert fh.read() == FASTA


def test_parse_fasta(tmp_path):
    p = tmp_path / "seq.fasta"
    p.write_text(FASTA)
    seqs, lengths_df = qc.parse_fasta(str(p))
    assert seqs["c1"] == "ATGCATGCAT"
    assert lengths_df.loc["c1", "lengths"] == 10
    assert lengths_df.loc["c3", "lengths"] == 2


# ---------------------------------------------------------------------------
# Statistics (driven off the parsed lengths DataFrame)
# ---------------------------------------------------------------------------

def _fixtures(tmp_path):
    p = tmp_path / "seq.fasta"
    p.write_text(FASTA)
    return qc.parse_fasta(str(p))


def test_contig_count_respects_threshold(tmp_path):
    _, df = _fixtures(tmp_path)
    assert qc.contig_count(df, 0) == 3
    assert qc.contig_count(df, 3) == 2  # c3 (len 2) excluded


def test_largest_contig_length(tmp_path):
    _, df = _fixtures(tmp_path)
    assert qc.largest_contig_length(df) == 10


def test_total_length_respects_threshold(tmp_path):
    _, df = _fixtures(tmp_path)
    assert qc.total_length(df, 0) == 22
    assert qc.total_length(df, 3) == 20


def test_n_amount(tmp_path):
    seqs, df = _fixtures(tmp_path)
    assert qc.n_amount(df, seqs, 0) == 4
    assert qc.n_amount(df, seqs, 3) == 4  # the N-bearing contig is length 10


def test_gc_content(tmp_path):
    seqs, df = _fixtures(tmp_path)
    assert qc.gc_content(df, seqs, 0) == 45.45   # 10 / 22 * 100
    assert qc.gc_content(df, seqs, 3) == 50.0     # 10 / 20 * 100


def test_n_l_statistics(tmp_path):
    _, df = _fixtures(tmp_path)
    # lengths sorted desc [10, 10, 2], total 22; 50% target 11 -> 2nd contig.
    assert qc.n_l_statistics(df, 50, 0) == [10, 2]


def test_n_l_statistics_empty_when_threshold_excludes_all(tmp_path):
    _, df = _fixtures(tmp_path)
    assert qc.n_l_statistics(df, 50, 1000) == [None, None]
