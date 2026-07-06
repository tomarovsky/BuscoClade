"""Unit tests for workflow/scripts/phylip_tree_namefix.py.

Covers loading the Short_ID/Original_ID map and restoring original taxon names
in a NEWICK tree produced from a strict-PHYLIP run.
"""
import pytest

import phylip_tree_namefix as ptn


# ---------------------------------------------------------------------------
# Map loading
# ---------------------------------------------------------------------------

def test_load_mapping_parses_valid_and_skips_malformed(tmp_path):
    m = tmp_path / "concat.phy.map"
    m.write_text(
        "Short_ID\tOriginal_ID\n"
        "000001\tHomo_sapiens\n"
        "malformed_line_without_tab\n"
        "000002\tPan\n"
    )
    mapping = ptn.load_mapping(str(m))
    assert mapping == {"000001": "Homo_sapiens", "000002": "Pan"}


def test_load_mapping_missing_file_exits(tmp_path):
    with pytest.raises(SystemExit):
        ptn.load_mapping(str(tmp_path / "nope.map"))


def test_load_mapping_header_only_exits(tmp_path):
    m = tmp_path / "empty.map"
    m.write_text("Short_ID\tOriginal_ID\n")
    with pytest.raises(SystemExit):
        ptn.load_mapping(str(m))


# ---------------------------------------------------------------------------
# End-to-end rename
# ---------------------------------------------------------------------------

def test_main_replaces_ids_in_tree(tmp_path, monkeypatch):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(000001:0.1,000002:0.2);\n")
    m = tmp_path / "tree.map"
    m.write_text("Short_ID\tOriginal_ID\n000001\tHomo_sapiens\n000002\tPan\n")
    out = tmp_path / "renamed.nwk"

    monkeypatch.setattr(
        "sys.argv",
        ["phylip_tree_namefix.py", "-i", str(tree), "-m", str(m), "-o", str(out)],
    )
    ptn.main()

    result = out.read_text()
    assert "Homo_sapiens" in result
    assert "Pan" in result
    assert "000001" not in result
    assert "000002" not in result
