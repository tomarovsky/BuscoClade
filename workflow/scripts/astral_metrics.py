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
    with open(newick_file, "r") as file:
        tree_string = ""
        newick = file.readline().replace("_", " ").strip().split("'")
        tree_string += newick[0]
        for i in range(1, len(newick), 2):
            line = ""
            flag = True
            for s in newick[i + 1]:
                if s == ")" or s == ",":
                    if flag:
                        nhx = newick[i].replace(",", ".").replace(";", ":")[1:]
                        line += f"[&&NHX:{nhx}{s}"
                        flag = False
                    else:
                        line += s
                else:
                    line += s
            tree_string += line
        return tree_string


def add_node_numbers(tree):
    node_counter = 1
    for node in tree.traverse():
        if not node.is_leaf():
            node.node_number = node_counter
            node_counter += 1


def main():
    parser = argparse.ArgumentParser(description='Extract node statistics from phylogenetic tree')
    parser.add_argument('-i', '--input', required=True, help='Input tree file (Newick format)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-g', '--outgroup', default=False, help="outgroup species name (default = unrooted)")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)

    try:
        # Load tree
        tree = Tree(newick_to_nhx(args.input))

        # Set outgroup
        if args.outgroup:
            outgroup_names = [name.strip() for name in args.outgroup.split(',')]
            target_nodes = []
            for name in outgroup_names:
                node = tree.search_nodes(name=name)
                if node:
                    target_nodes.append(node[0])
                else:
                    print(f"Warning: Outgroup species '{name}' not found. Skipping.")
            if target_nodes:
                if len(target_nodes) == 1:
                    tree.set_outgroup(target_nodes[0])
                    print(f"Rooted using: {target_nodes[0].name}")
                else:
                    mrca = tree.get_common_ancestor(target_nodes)
                    try:
                        tree.set_outgroup(target_nodes)
                        outgroup_list = [node.name for node in target_nodes]
                        print(f"Rooted using group (MRCA of {len(outgroup_list)} species): {', '.join(outgroup_list)}")
                    except Exception as e:
                        print(f"Error during set_outgroup: {e}. Attempting manual MRCA rooting.")
                        if mrca != tree:
                            tree.set_outgroup(mrca)
                            print(f"Rooted manually using MRCA of {len(target_nodes)} species.")
                        else:
                            print("Warning: MRCA is the current root. Cannot root, unrooting.")
                            tree.unroot()

            else:
                print("Warning: No outgroup species found. Unrooting.")
                tree.unroot()
        else:
            tree.unroot()

        # Ladderize
        tree.ladderize(direction=True)

        # Add node numbers
        add_node_numbers(tree)

        # Extract node statistics
        metrics = ["q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3", "QC", "EN"]
        stats_data = []

        for node in tree.traverse():
            if not node.is_leaf() and hasattr(node, "node_number"):
                node_stats = {"node": node.node_number}
                for metric in metrics:
                    value = getattr(node, metric, "")
                    try:
                        value = float(value)
                    except (ValueError, TypeError):
                        value = ""
                    node_stats[metric] = value
                stats_data.append(node_stats)

        # Sort by node number
        stats_data.sort(key=lambda x: x["node"])

        with open(args.output, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["node"] + metrics, delimiter="\t")
            writer.writeheader()
            for row in stats_data:
                writer.writerow(row)

        print(f"Successfully extracted statistics for {len(stats_data)} nodes")
        print(f"Output saved to: {args.output}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
