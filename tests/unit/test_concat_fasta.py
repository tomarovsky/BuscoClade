"""Unit tests for workflow/scripts/concat_fasta.py.

Verifies that per-gene FASTA files are concatenated by shared header, in file
order, and that empty inputs are tolerated (dropped genes on prank routes).
"""
from types import SimpleNamespace

from Bio import SeqIO

import concat_fasta as cf


def _fasta_dict(path):
    return {r.id: str(r.seq) for r in SeqIO.parse(str(path), "fasta")}


def test_main_concatenates_by_header_in_input_order(tmp_path):
    f1 = tmp_path / "gene1.fna"
    f2 = tmp_path / "gene2.fna"
    f1.write_text(">sp1\nAAA\n>sp2\nCCC\n")
    f2.write_text(">sp1\nGGG\n>sp2\nTTT\n")
    out = tmp_path / "concat.fna"

    cf.args = SimpleNamespace(input=[str(f1), str(f2)], output=str(out))
    cf.main()

    assert _fasta_dict(out) == {"sp1": "AAAGGG", "sp2": "CCCTTT"}


def test_main_tolerates_empty_input_files(tmp_path):
    f1 = tmp_path / "gene1.fna"
    empty = tmp_path / "gene_empty.fna"
    f1.write_text(">sp1\nAAA\n>sp2\nCCC\n")
    empty.write_text("")  # e.g. a prank gene that timed out
    out = tmp_path / "concat.fna"

    cf.args = SimpleNamespace(input=[str(f1), str(empty)], output=str(out))
    cf.main()

    assert _fasta_dict(out) == {"sp1": "AAA", "sp2": "CCC"}
