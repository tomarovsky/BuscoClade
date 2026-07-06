#!/usr/bin/env python3
__author__ = 'tomarovsky'

"""
gene_count_report.py

Summarizes, per species, how many single-copy BUSCO genes it contributed and how
many survived into the common (intersection) set used to build the supermatrix.

This makes gene loss visible — a reconstructed sample (VCF- or consensus-derived)
that reconstructs poorly, or any species with low BUSCO completeness, shrinks
`common_ids` for everyone. The species with `single_copy` closest to the common
total are the bottleneck limiting the tree's data.

Usage:
    gene_count_report.py \\
        --single_copy_dir results/ids/species_ids/single_copy \\
        --multi_copy_dir  results/ids/species_ids/multi_copy \\
        --common          results/ids/common_ids/common.ids \\
        --output          results/ids/gene_counts.tsv
"""

import argparse
import sys
from pathlib import Path


def read_ids(path: Path) -> set:
    if not path.is_file():
        return set()
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def classify_source(species: str) -> str:
    if species.endswith(".AltRef"):
        return "reconstruct:vcf"
    if species.endswith(".Consensus"):
        return "reconstruct:consensus"
    return "genome"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single_copy_dir", required=True, type=Path)
    parser.add_argument("--multi_copy_dir", required=True, type=Path)
    parser.add_argument("--common", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    common = read_ids(args.common)
    n_common = len(common)

    rows = []
    for sc_file in sorted(args.single_copy_dir.glob("*.ids")):
        species = sc_file.stem
        sc = read_ids(sc_file)
        mc = read_ids(args.multi_copy_dir / f"{species}.ids")
        in_common = len(sc & common)
        rows.append({
            "species": species,
            "source": classify_source(species),
            "single_copy": len(sc),
            "multi_copy": len(mc),
            "in_common": in_common,
            "single_copy_not_in_common": len(sc) - in_common,
        })

    # Smallest single_copy first — these species constrain the common set.
    rows.sort(key=lambda r: r["single_copy"])

    columns = ["species", "source", "single_copy", "multi_copy",
               "in_common", "single_copy_not_in_common"]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write(f"# common_genes\t{n_common}\n")
        out.write("\t".join(columns) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in columns) + "\n")

    print(f"Common single-copy BUSCO genes (supermatrix size): {n_common}")
    if rows:
        bottleneck = rows[0]
        print(
            f"Fewest single-copy genes: {bottleneck['species']} "
            f"({bottleneck['single_copy']}, source={bottleneck['source']})"
        )
    print(f"Wrote per-species gene counts to {args.output}")


if __name__ == "__main__":
    main()
