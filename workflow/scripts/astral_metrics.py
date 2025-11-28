#!/usr/bin/env python3
__author__ = 'tomarovsky'

"""
Script for extracting node statistics from phylogenetic trees
"""

import argparse
import csv
import os
import sys

from ete3 import Tree


def newick_to_nhx(newick_file) -> str:
    """Convert Newick format to NHX format"""
    with open(newick_file, 'r') as file:
        tree_string = ''
        newick = file.readline().replace("_", " ").strip().split("'")
        tree_string += newick[0]
        for i in range(1, len(newick), 2):
            line = ''
            flag = True
            for s in newick[i+1]:
                if s == ")" or s == ",":
                    if flag:
                        nhx = newick[i].replace(',', '.').replace(';', ':')[1:]
                        line += f"[&&NHX:{nhx}{s}"
                        flag = False
                    else:
                        line += s
                else:
                    line += s
            tree_string += line
        return tree_string


def main():
    parser = argparse.ArgumentParser(description='Extract node statistics from phylogenetic tree')
    parser.add_argument('-i', '--input', required=True, help='Input tree file (Newick format)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)

    metrics = ["q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN"]

    # Load tree
    tree = Tree(newick_to_nhx(args.input))

    # Ladderize
    tree.ladderize(direction=True)

    # Extract node statistics
    stats_data = []

    for node in tree.traverse():
        if not node.is_leaf() and hasattr(node, "node_number"):
            node_stats = {"node": node.node_number}
            for metric in metrics:
                value = getattr(node, metric, "")
                try:
                    value = float(value)
                except:
                    value = ""
                node_stats[metric] = value
            stats_data.append(node_stats)

    # Sort by node number and write to TSV
    stats_data.sort(key=lambda x: x["node"])

    with open(args.output, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["node"] + metrics, delimiter="\t")
        writer.writeheader()
        for row in stats_data:
            writer.writerow(row)

    print(f"Node statistics saved to {args.output}")


if __name__ == "__main__":
    main()
