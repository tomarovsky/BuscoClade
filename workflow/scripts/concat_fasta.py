#!/usr/bin/env python3
__author__ = "tomarovsky"
from collections import defaultdict
from Bio import SeqIO
import argparse


def main():
    sequence_map = defaultdict(str)
    for i in args.input:
        for sequence in SeqIO.parse(i, "fasta"):
            sequence_map[sequence.name] += str(sequence.seq)

    outfile = open(args.output, "w")
    for key, value in sequence_map.items():
        outfile.write(f">{key}\n{value}\n")
    outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="script for concatenate FASTA files into one big FASTA file " "by concatenating sequences with the same identifier"
    )
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-i", "--input", type=str, nargs="+", help="input list of concat FASTA files with the same headers")
    group_required.add_argument("-o", "--output", type=str, help="output FASTA file name")
    args = parser.parse_args()
    main()
