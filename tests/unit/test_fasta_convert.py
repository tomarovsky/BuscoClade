"""Unit tests for workflow/scripts/fasta_convert.py.

Covers the FASTA -> {phylip, stockholm, nexus} conversions, focusing on the
PHYLIP numeric-ID remapping (and its .map sidecar) and the NEXUS MrBayes-block
requirement.
"""
from types import SimpleNamespace

import pytest

import fasta_convert as fc


def _write_alignment(tmp_path):
    fa = tmp_path / "concat.fna"
    # Equal-length sequences so phylip/nexus alignment writers accept them.
    fa.write_text(">Homo_sapiens\nACGTACGT\n>Pan\nACGTAAGT\n")
    return fa


# ---------------------------------------------------------------------------
# PHYLIP
# ---------------------------------------------------------------------------

def test_convert_phylip_writes_zero_padded_ids_and_map(tmp_path):
    fa = _write_alignment(tmp_path)
    out = tmp_path / "concat.phy"
    fc.convert_phylip(SimpleNamespace(input=str(fa), output=str(out)))

    map_lines = (tmp_path / "concat.phy.map").read_text().splitlines()
    assert map_lines[0] == "Short_ID\tOriginal_ID"
    assert map_lines[1] == "000001\tHomo_sapiens"
    assert map_lines[2] == "000002\tPan"

    # The PHYLIP body carries the remapped numeric IDs, not the original names.
    body = out.read_text()
    assert "000001" in body
    assert "Homo_sapiens" not in body


def test_convert_phylip_empty_input_exits(tmp_path):
    fa = tmp_path / "empty.fna"
    fa.write_text("")
    with pytest.raises(SystemExit):
        fc.convert_phylip(SimpleNamespace(input=str(fa), output=str(tmp_path / "x.phy")))


# ---------------------------------------------------------------------------
# Stockholm
# ---------------------------------------------------------------------------

def test_convert_stockholm(tmp_path):
    fa = _write_alignment(tmp_path)
    out = tmp_path / "concat.sth"
    fc.convert_stockholm(SimpleNamespace(input=str(fa), output=str(out)))
    text = out.read_text()
    assert text.startswith("# STOCKHOLM")
    assert "Homo_sapiens" in text


# ---------------------------------------------------------------------------
# NEXUS
# ---------------------------------------------------------------------------

def test_convert_nexus_appends_block(tmp_path):
    fa = _write_alignment(tmp_path)
    block = tmp_path / "mrbayes.block"
    block.write_text("begin mrbayes;\n  lset nst=6;\nend;\n")
    out = tmp_path / "concat.nex"

    fc.convert_nexus(SimpleNamespace(input=str(fa), output=str(out), block=str(block)))
    text = out.read_text()
    assert "#NEXUS" in text
    assert "begin mrbayes;" in text


def test_convert_nexus_requires_block(tmp_path):
    fa = _write_alignment(tmp_path)
    with pytest.raises(SystemExit):
        fc.convert_nexus(SimpleNamespace(input=str(fa), output=str(tmp_path / "x.nex"), block=None))
