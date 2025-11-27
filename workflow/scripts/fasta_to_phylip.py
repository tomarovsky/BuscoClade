#!/usr/bin/env python3
__author__ = "tomarovsky"

import argparse
import sys

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to strict PHYLIP with zero-padded numeric IDs")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")
    args = parser.parse_args()

    records = []
    mapping = []

    try:
        for i, record in enumerate(SeqIO.parse(args.input, "fasta"), 1):
            original_id = record.id

            new_id = f"{i:06}"

            mapping.append(f"{new_id}\t{original_id}")

            record.id = new_id
            record.description = ""
            records.append(record)

        if not records:
            sys.exit("Input file is empty or has an invalid format.")

        count = SeqIO.write(records, args.output, "phylip")
        print(f"Done! Converted {count} sequences to {args.output}")

        map_file = args.output + ".map"
        with open(map_file, "w") as f:
            f.write("Short_ID\tOriginal_ID\n")
            f.write("\n".join(mapping))
            f.write("\n")
        print(f"Mapping saved to {map_file}")

    except FileNotFoundError:
        sys.exit(f"Error: File {args.input} not found.")
    except Exception as e:
        sys.exit(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
