"""Unit tests for workflow/scripts/merge_common_ids.py.

Verifies that per-species single-copy sequences for each common BUSCO id are
merged into one FASTA per gene, headed by the species name (derived from the
grandparent directory of each single_copy_busco_sequences dir).
"""
from types import SimpleNamespace

import merge_common_ids as mci


def _species_dir(root, species):
    d = root / species / "busco_sequences" / "single_copy_busco_sequences"
    d.mkdir(parents=True)
    return d


def test_main_merges_species_by_gene(tmp_path):
    sp1 = _species_dir(tmp_path, "sp1")
    sp2 = _species_dir(tmp_path, "sp2")

    # Each BUSCO file has a header line (dropped) followed by the sequence.
    (sp1 / "gene1.fna").write_text(">h\nAAAA\n")
    (sp1 / "gene1.faa").write_text(">h\nKK\n")
    (sp2 / "gene1.fna").write_text(">h\nCCCC\n")
    (sp2 / "gene1.faa").write_text(">h\nLL\n")

    common = tmp_path / "common.ids"
    common.write_text("gene1\n")

    outdir = tmp_path / "merged"  # must not pre-exist (script calls mkdir())

    mci.args = SimpleNamespace(
        common_ids=str(common),
        single_copy_files=[str(sp1), str(sp2)],
        outdir=str(outdir),
    )
    mci.main()

    assert (outdir / "gene1.merged.fna").read_text() == ">sp1\nAAAA\n>sp2\nCCCC\n"
    assert (outdir / "gene1.merged.faa").read_text() == ">sp1\nKK\n>sp2\nLL\n"
