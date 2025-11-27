#!/usr/bin/env python3
__author__ = "tomarovsky"
import argparse
import sys


def load_mapping(map_file_path):
    mapping = {}

    try:
        with open(map_file_path, "r") as f:
            next(f)

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) == 2:
                    short_id = parts[0].strip()
                    original_id = parts[1].strip()
                    mapping[short_id] = original_id
                else:
                    print(f"Warning: Skipping malformed line in map file: {line}", file=sys.stderr)

    except FileNotFoundError:
        sys.exit(f"Error: Map file not found at {map_file_path}")

    if not mapping:
        sys.exit(f"Error: Map file {map_file_path} is empty or has no valid entries.")

    return mapping

def main():
    parser = argparse.ArgumentParser(description="Rename taxa in NEWICK tree using a .map file.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input NEWICK tree file.")
    parser.add_argument("-m", "--mapfile", type=str, required=True, help="Input .map file (Short_ID <tab> Original_ID).")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output NEWICK tree file.")
    args = parser.parse_args()

    print(f"Loading mapping data from {args.mapfile}...")
    name_map = load_mapping(args.mapfile)

    try:
        with open(args.input, "r") as file:
            nwk_tree = file.read().replace("\n", "").strip()
    except FileNotFoundError:
        sys.exit(f"Error: NEWICK input file not found at {args.input}")

    print("Performing replacements in NEWICK tree...")
    replacements_count = 0

    for short_id, full_id in sorted(name_map.items(), key=lambda item: len(item[0]), reverse=True):
        if short_id in nwk_tree:
            nwk_tree = nwk_tree.replace(short_id, full_id)
            replacements_count += 1

    if replacements_count < len(name_map):
        print(f"Warning: Only {replacements_count} taxa were renamed, but {len(name_map)} entries were in the map file. Check if all taxa from map are present in the tree.", file=sys.stderr)

    try:
        with open(args.output, "w") as nwk_namefix_tree:
            nwk_namefix_tree.write(nwk_tree)
        print(f"Success! Renamed {replacements_count} taxa and saved the new tree to {args.output}")
    except IOError as e:
        sys.exit(f"Error writing output file {args.output}: {e}")


if __name__ == "__main__":
    main()
