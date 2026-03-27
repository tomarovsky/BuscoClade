#!/usr/bin/env python3
__author__ = "tomarovsky"
import argparse
from collections import defaultdict
from Bio import SeqIO

def main():
    sequence_map = defaultdict(str)
    empty_files = 0
    non_empty_files = 0
    for i in args.input:
        records = list(SeqIO.parse(i, "fasta"))
        if records:
            non_empty_files += 1
            for sequence in records:
                sequence_map[sequence.name] += str(sequence.seq)
        else:
            empty_files += 1

    print(f"Total input files: {len(args.input)}")
    print(f"Non-empty files: {non_empty_files}")
    print(f"Empty files: {empty_files}")

    outfile = open(args.output, "w")
    for key, value in sequence_map.items():
        outfile.write(f">{key}\n{value}\n")
    outfile.close()

    lengths = set(len(v) for v in sequence_map.values())
    if len(lengths) == 1:
        print(f"Sequence length in output FASTA: {lengths.pop()}")
    else:
        print(f"WARNING: sequences have different lengths: {sorted(lengths)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="script for concatenate FASTA files into one big FASTA file by concatenating sequences with the same identifier"
    )
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-i", "--input", type=str, nargs="+", help="input list of concat FASTA files with the same headers")
    group_required.add_argument("-o", "--output", type=str, help="output FASTA file name")
    args = parser.parse_args()
    main()
