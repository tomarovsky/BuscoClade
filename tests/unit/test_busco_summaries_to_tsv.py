"""Unit tests for workflow/scripts/busco_summaries_to_tsv.py.

Covers parsing the BUSCO short_summary results line into the S/D/F/M/N columns
and the species-name derivation from the filename.
"""
from types import SimpleNamespace

import pandas as pd

import busco_summaries_to_tsv as bs


# A realistic BUSCO short_summary layout: the C:...[S,D],F,M,n results string
# sits on line index 8 (the 9th line), which is what the script reads.
SHORT_SUMMARY = (
    "# BUSCO version is: 6.0.0\n"
    "# The lineage dataset is: eukaryota_odb10\n"
    "# Summarized benchmarking in BUSCO notation for file genome.fna\n"
    "# BUSCO was run in mode: genome\n"
    "# Gene predictor used: metaeuk\n"
    "\n"
    "\t***** Results: *****\n"
    "\n"
    "\tC:98.6%[S:97.2%,D:1.4%],F:0.5%,M:0.9%,n:255\n"
    "\t251\tComplete BUSCOs (C)\n"
)


def test_main_extracts_sdfmn_and_species(tmp_path):
    # Filename must be short_summary_<Genus_species>.txt for the [14:-4] slice.
    f = tmp_path / "short_summary_Homo_sapiens.txt"
    f.write_text(SHORT_SUMMARY)
    out = tmp_path / "busco_summaries.tsv"

    bs.args = SimpleNamespace(input=[str(f)], output=str(out))
    bs.main()

    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["Species", "S", "D", "F", "M", "N"]
    row = df.iloc[0]
    assert row["Species"] == "Homo sapiens"
    assert row["S"] == 97.2
    assert row["D"] == 1.4
    assert row["F"] == 0.5
    assert row["M"] == 0.9
    assert row["N"] == 255


def test_main_skips_incomplete_files(tmp_path):
    good = tmp_path / "short_summary_Homo_sapiens.txt"
    good.write_text(SHORT_SUMMARY)
    truncated = tmp_path / "short_summary_Pan_troglodytes.txt"
    truncated.write_text("# only a couple of lines\n# nothing useful\n")
    out = tmp_path / "busco_summaries.tsv"

    bs.args = SimpleNamespace(input=[str(good), str(truncated)], output=str(out))
    bs.main()

    df = pd.read_csv(out, sep="\t")
    assert list(df["Species"]) == ["Homo sapiens"]
