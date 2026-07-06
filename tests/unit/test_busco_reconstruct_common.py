"""Unit tests for workflow/scripts/busco_reconstruct_common.py.

These lock down the subtle, easy-to-break parts shared by apply_vcf_to_busco.py
and apply_consensus_to_busco.py: exon-aware coordinate mapping (especially the
minus-strand reverse-complement offset), IUPAC/translation, MetaEuk header
parsing, and CDS-exon resolution from GFF.
"""
import textwrap

import pytest

import busco_reconstruct_common as brc


# ---------------------------------------------------------------------------
# Small primitives
# ---------------------------------------------------------------------------

def test_reverse_complement_basic():
    assert brc.reverse_complement("AAGG") == "CCTT"
    assert brc.reverse_complement("ACGT") == "ACGT"  # palindrome


def test_reverse_complement_iupac_and_case():
    # R (A/G) complements to Y (C/T); case is preserved.
    assert brc.reverse_complement("aR") == "Yt"


@pytest.mark.parametrize("bases,expected", [
    ({"A"}, "A"),
    ({"A", "G"}, "R"),
    ({"C", "T"}, "Y"),
    ({"A", "C", "G", "T"}, "N"),
])
def test_bases_to_iupac(bases, expected):
    assert brc.bases_to_iupac(bases) == expected


def test_translate_sequence_stops_at_stop_codon():
    assert brc.translate_sequence("ATGGCC") == "MA"
    assert brc.translate_sequence("ATGTAA") == "M"   # TAA stop -> break


def test_translate_ambiguous_codon():
    assert brc.translate_ambiguous_codon("GCN") == "A"   # all 4 -> Ala
    assert brc.translate_ambiguous_codon("ATR") == "X"   # Ile / Met disagree -> X


# ---------------------------------------------------------------------------
# MetaEuk header parsing
# ---------------------------------------------------------------------------

def test_parse_metaeuk_fna_header():
    header = ">8689at40674_1868482_0:0013eb|NC_037328.1:120896888-120904250|-"
    info = brc.parse_metaeuk_fna_header(header)
    assert info["chrom"] == "NC_037328.1"
    assert info["start0"] == 120896888
    assert info["end0"] == 120904250
    assert info["strand"] == "-"


def test_parse_metaeuk_fna_header_rejects_bad_format():
    with pytest.raises(ValueError):
        brc.parse_metaeuk_fna_header(">no_pipes_here")


# ---------------------------------------------------------------------------
# CDS exon resolution from GFF
# ---------------------------------------------------------------------------

def _write_gff(tmp_path):
    gff = tmp_path / "genes.gff"
    gff.write_text(textwrap.dedent("""\
        # comment line
        chr1\tMetaEuk\tCDS\t10\t12\t.\t+\t.\tTarget_ID=gene1;extra=x
        chr1\tMetaEuk\tgene\t10\t22\t.\t+\t.\tTarget_ID=gene1
        chr1\tMetaEuk\tCDS\t20\t22\t.\t+\t.\tTarget_ID=gene1;extra=y
    """))
    return gff


def test_load_cds_index_only_keeps_cds_rows(tmp_path):
    index = brc.load_cds_index(_write_gff(tmp_path))
    assert index == {("gene1", "chr1"): [(10, 12), (20, 22)]}


def test_get_cds_exons_filters_to_gene_range(tmp_path):
    index = brc.load_cds_index(_write_gff(tmp_path))
    exons, source = brc.get_cds_exons("gene1", "chr1", 10, 22, index, {})
    assert exons == [(10, 12), (20, 22)]
    assert source == "rerun"


def test_get_cds_exons_falls_back_to_initial(tmp_path):
    index = brc.load_cds_index(_write_gff(tmp_path))
    exons, source = brc.get_cds_exons("gene1", "chr1", 10, 22, {}, index)
    assert source == "initial"


def test_get_cds_exons_raises_when_absent():
    with pytest.raises(ValueError):
        brc.get_cds_exons("ghost", "chr1", 1, 100, {}, {})


# ---------------------------------------------------------------------------
# Exon-aware SNP application — the trickiest logic
# ---------------------------------------------------------------------------

def test_apply_snps_plus_strand_single_exon():
    # exon 10..13 (len 4); SNP at genomic 11 -> FNA index 1, base as-is.
    out = brc.apply_snps_exon_aware("AAAA", "+", [(10, 13)], {0: {11: "G"}})
    assert out == "AGAA"


def test_apply_snps_minus_strand_single_exon():
    # Minus strand: FNA index 0 <-> genomic 13, index 2 <-> genomic 11.
    # A SNP at genomic 11 (forward base G) lands at FNA index 2 and is complemented.
    out = brc.apply_snps_exon_aware("TTTT", "-", [(10, 13)], {0: {11: "G"}})
    assert out == "TTCT"


def test_apply_snps_plus_strand_multi_exon_offsets():
    # exon0 10..12 (len3), exon1 20..22 (len3); SNPs 11->C (idx1), 22->T (idx5).
    out = brc.apply_snps_exon_aware(
        "AAAAAA", "+", [(10, 12), (20, 22)], {0: {11: "C"}, 1: {22: "T"}}
    )
    assert out == "ACAAAT"


def test_apply_snps_minus_strand_multi_exon_offsets():
    # Minus-strand exons are listed highest-coord first (MetaEuk order).
    # exon0 20..22: SNP 21 (fwd G) -> idx1, complemented to C.
    # exon1 10..12: SNP 11 (fwd C) -> idx4, complemented to G.
    out = brc.apply_snps_exon_aware(
        "AAAAAA", "-", [(20, 22), (10, 12)], {0: {21: "G"}, 1: {11: "C"}}
    )
    assert out == "ACAAGA"


def test_apply_snps_no_changes_returns_uppercased_original():
    out = brc.apply_snps_exon_aware("acgt", "+", [(1, 4)], {})
    assert out == "ACGT"
