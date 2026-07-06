"""Unit tests for workflow/scripts/add_altref_to_alignment.py.

Covers the gap-aware insertion logic that places AltRef sequences into a raw
alignment by copying the reference's gap pattern (valid only because
reconstruction is substitution-only, so ungapped lengths match).
"""
from pathlib import Path

import add_altref_to_alignment as aa


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def test_read_fasta_uses_bare_id_and_joins_wrapped_lines(tmp_path):
    fa = tmp_path / "in.fna"
    fa.write_text(">sp1 some annotation\nAC\nGT\n>sp2\nTTTT\n")
    seqs = aa.read_fasta(fa)
    assert seqs == {"sp1": "ACGT", "sp2": "TTTT"}


def test_write_then_read_fasta_roundtrip(tmp_path):
    seqs = {"a": "A" * 130, "b": "CGT"}
    out = tmp_path / "out.fna"
    aa.write_fasta(seqs, out)
    # 60 chars per line for the long sequence.
    lines = out.read_text().splitlines()
    assert lines[0] == ">a"
    assert len(lines[1]) == 60
    assert aa.read_fasta(out) == seqs


# ---------------------------------------------------------------------------
# Gap insertion
# ---------------------------------------------------------------------------

def test_insert_gaps_places_dashes_at_positions():
    assert aa.insert_gaps("ACGT", [1, 3], 6) == "A-C-GT"


def test_insert_gaps_no_gaps_returns_sequence():
    assert aa.insert_gaps("ACGT", [], 4) == "ACGT"


def test_get_gap_positions():
    assert aa.get_gap_positions("A-C-GT") == [1, 3]
    assert aa.get_gap_positions("ACGT") == []


# ---------------------------------------------------------------------------
# AltRef gene file discovery
# ---------------------------------------------------------------------------

def _make_altref_gene(busco_dir: Path, species: str, bare_id: str, seq: str) -> Path:
    d = busco_dir / species / "busco_sequences" / "single_copy_busco_sequences"
    d.mkdir(parents=True)
    f = d / f"{bare_id}.fna"
    f.write_text(f">{bare_id}\n{seq}\n")
    return f


def test_find_altref_gene_file_strips_merged_suffix(tmp_path):
    busco = tmp_path / "busco"
    expected = _make_altref_gene(busco, "S1.ref1.AltRef", "geneX", "AAGG")
    found = aa.find_altref_gene_file(busco, "S1.ref1.AltRef", "geneX.merged")
    assert found == expected


def test_find_altref_gene_file_missing_returns_none(tmp_path):
    busco = tmp_path / "busco"
    busco.mkdir()
    assert aa.find_altref_gene_file(busco, "S1.ref1.AltRef", "absent.merged") is None


# ---------------------------------------------------------------------------
# add_altref_sequences — the integration point
# ---------------------------------------------------------------------------

def test_add_altref_sequences_inserts_with_ref_gap_pattern(tmp_path):
    busco = tmp_path / "busco"
    # ref1 aligned as "A-CGT-": gaps at positions 1 and 5, ungapped "ACGT" (len 4).
    _make_altref_gene(busco, "S1.ref1.AltRef", "geneX", "AAGG")  # len 4 -> matches
    raw_aln = {"ref1": "A-CGT-", "sp2": "AACGTA"}

    extended = aa.add_altref_sequences(
        raw_aln=raw_aln,
        busco_dir=busco,
        ref_to_altrefs={"ref1": ["S1.ref1.AltRef"]},
        gene_id="geneX.merged",
    )

    # AAGG with gaps copied from ref1's pattern (positions 1 and 5), uppercased.
    assert extended["S1.ref1.AltRef"] == "A-AGG-"
    # original sequences are preserved.
    assert extended["ref1"] == "A-CGT-"
    assert extended["sp2"] == "AACGTA"


def test_add_altref_sequences_skips_length_mismatch(tmp_path):
    busco = tmp_path / "busco"
    # ungapped ref length is 4, but the AltRef gene is length 5 -> skipped.
    _make_altref_gene(busco, "S1.ref1.AltRef", "geneX", "AAGGT")
    raw_aln = {"ref1": "A-CGT-"}

    extended = aa.add_altref_sequences(
        raw_aln=raw_aln,
        busco_dir=busco,
        ref_to_altrefs={"ref1": ["S1.ref1.AltRef"]},
        gene_id="geneX.merged",
    )
    assert "S1.ref1.AltRef" not in extended


def test_add_altref_sequences_skips_missing_reference(tmp_path):
    busco = tmp_path / "busco"
    busco.mkdir()
    raw_aln = {"sp2": "ACGT"}
    extended = aa.add_altref_sequences(
        raw_aln=raw_aln,
        busco_dir=busco,
        ref_to_altrefs={"refMissing": ["S1.refMissing.AltRef"]},
        gene_id="geneX.merged",
    )
    assert extended == {"sp2": "ACGT"}


# ---------------------------------------------------------------------------
# Reference removal
# ---------------------------------------------------------------------------

def test_remove_ref_sequences_drops_only_refs():
    aln = {"ref1": "AAAA", "sp2": "CCCC", "S1.ref1.AltRef": "GGGG"}
    filtered = aa.remove_ref_sequences(aln, ["ref1"])
    assert filtered == {"sp2": "CCCC", "S1.ref1.AltRef": "GGGG"}
