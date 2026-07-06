"""Unit tests for workflow/scripts/apply_consensus_to_busco.py.

Covers the consensus-genome-facing logic: collecting required contigs/extents
from BUSCO FNA headers, the coordinate-identity guard, and exon slicing
(including minus-strand reverse-complement) from the consensus genome.
"""
import pytest

pysam = pytest.importorskip("pysam")

import apply_consensus_to_busco as acb


def _make_fasta(tmp_path):
    """Write and faidx a small consensus genome; return an open FastaFile."""
    p = tmp_path / "consensus.fasta"
    p.write_text(">chr1\nACGTACGTAC\n>chr2\nTTTT\n")
    pysam.faidx(str(p))
    return pysam.FastaFile(str(p))


# ---------------------------------------------------------------------------
# Required-contig collection
# ---------------------------------------------------------------------------

def test_collect_required_contigs_tracks_max_extent(tmp_path):
    d = tmp_path / "single_copy"
    d.mkdir()
    (d / "g1.fna").write_text(">g1|chr1:0-9|+\nACGT\n")
    (d / "g2.fna").write_text(">g2|chr1:5-19|-\nACGT\n")
    (d / "g3.fna").write_text(">g3|chr2:0-4|+\nAC\n")

    # end1 = end0 + 1; chr1's max extent comes from g2 (19 -> 20).
    assert acb.collect_required_contigs(d) == {"chr1": 20, "chr2": 5}


# ---------------------------------------------------------------------------
# Coordinate-identity guard
# ---------------------------------------------------------------------------

def test_validate_consensus_contigs_ok(tmp_path):
    fasta = _make_fasta(tmp_path)
    # chr1 has length 10; requiring up to 5 is fine.
    acb.validate_consensus_contigs({"chr1": 5}, fasta, "consensus.fasta")  # no exit


def test_validate_consensus_contigs_too_short_exits(tmp_path):
    fasta = _make_fasta(tmp_path)
    with pytest.raises(SystemExit):
        acb.validate_consensus_contigs({"chr1": 20}, fasta, "consensus.fasta")


def test_validate_consensus_contigs_missing_exits(tmp_path):
    fasta = _make_fasta(tmp_path)
    with pytest.raises(SystemExit):
        acb.validate_consensus_contigs({"chrX": 1}, fasta, "consensus.fasta")


# ---------------------------------------------------------------------------
# CDS slicing from the consensus genome
# ---------------------------------------------------------------------------

def test_build_cds_from_consensus_plus_strand_single_exon(tmp_path):
    fasta = _make_fasta(tmp_path)
    # chr1 = ACGTACGTAC; exon 1..3 -> "ACG".
    assert acb.build_cds_from_consensus(fasta, "chr1", "+", [(1, 3)]) == "ACG"


def test_build_cds_from_consensus_plus_strand_multi_exon(tmp_path):
    fasta = _make_fasta(tmp_path)
    # exons 1..2 ("AC") + 5..6 ("AC") -> "ACAC".
    assert acb.build_cds_from_consensus(fasta, "chr1", "+", [(1, 2), (5, 6)]) == "ACAC"


def test_build_cds_from_consensus_minus_strand_reverse_complements(tmp_path):
    fasta = _make_fasta(tmp_path)
    # exon 1..3 "ACG" reverse-complemented -> "CGT".
    assert acb.build_cds_from_consensus(fasta, "chr1", "-", [(1, 3)]) == "CGT"


def test_build_cds_from_consensus_missing_chrom_raises(tmp_path):
    fasta = _make_fasta(tmp_path)
    with pytest.raises(ValueError):
        acb.build_cds_from_consensus(fasta, "chrX", "+", [(1, 3)])
