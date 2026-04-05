#!/usr/bin/env python3
"""
fasta_convert.py — convert a concatenated FASTA alignment to Stockholm, PHYLIP, or NEXUS.

Replaces fasta_to_stockholm.py, fasta_to_phylip.py, and fasta_to_nexus.py.

Usage
-----
# Stockholm
fasta_convert.py -i concat.fna -o concat.sth -f stockholm

# PHYLIP (strict, zero-padded numeric IDs; writes <output>.map alongside)
fasta_convert.py -i concat.fna -o concat.phy -f phylip

# NEXUS (appends a MrBayes block from --block)
fasta_convert.py -i concat.fna -o concat.nex -f nexus --block mrbayes.block
"""
__author__ = "tomarovsky"

import argparse
import sys
from Bio import AlignIO, SeqIO


def convert_stockholm(args):
    records = SeqIO.parse(args.input, "fasta")
    count = SeqIO.write(records, args.output, "stockholm")
    print(f"Written {count} sequences to {args.output}")


def convert_phylip(args):
    records = []
    mapping = []

    for i, record in enumerate(SeqIO.parse(args.input, "fasta"), 1):
        original_id = record.id
        new_id = f"{i:06}"
        mapping.append(f"{new_id}\t{original_id}")
        record.id = new_id
        record.description = ""
        records.append(record)

    if not records:
        sys.exit("Error: input file is empty or has an invalid format.")

    count = SeqIO.write(records, args.output, "phylip")
    print(f"Written {count} sequences to {args.output}")

    map_file = args.output + ".map"
    with open(map_file, "w") as f:
        f.write("Short_ID\tOriginal_ID\n")
        f.write("\n".join(mapping) + "\n")
    print(f"ID mapping saved to {map_file}")


def convert_nexus(args):
    if not args.block:
        sys.exit("Error: --block (MrBayes block file) is required for NEXUS output.")

    with open(args.output, "w") as outfile, open(args.block, "r") as blockfile:
        AlignIO.convert(args.input, "fasta", outfile, "nexus", "DNA")
        outfile.write("\n")
        for line in blockfile:
            outfile.write(line)
        outfile.write("\n")
    print(f"Written NEXUS alignment with MrBayes block to {args.output}")


FORMATS = {
    "stockholm": convert_stockholm,
    "phylip":    convert_phylip,
    "nexus":     convert_nexus,
}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    req = parser.add_argument_group("Required")
    req.add_argument("-i", "--input",  required=True, help="Input FASTA file")
    req.add_argument("-o", "--output", required=True, help="Output file")
    req.add_argument("-f", "--format", required=True, choices=FORMATS,
                     help="Output format: stockholm, phylip, or nexus")

    opt = parser.add_argument_group("NEXUS-specific")
    opt.add_argument("--block", help="Path to MrBayes block file (nexus only)")

    args = parser.parse_args()

    try:
        FORMATS[args.format](args)
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}")
    except Exception as e:
        sys.exit(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
