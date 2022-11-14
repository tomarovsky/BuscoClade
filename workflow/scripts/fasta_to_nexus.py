#!/usr/bin/env python3
from Bio import AlignIO
from sys import stdin
import argparse


def main():
    with open(args.output, "a") as outfile, open(args.block, "r") as blockfile:
        AlignIO.convert(args.input, "fasta", outfile, "nexus", args.type)
        outfile.write("\n")
        for line in blockfile:
            outfile.write(line)
        outfile.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for converting FASTA format to NEXUS format")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin, help="input concat FASTA file or stdin")
    group_required.add_argument('-t', '--type', type=str, help="molecular type (DNA, RNA or protein)")
    group_required.add_argument('-b', '--block', type=str, help="MrBayes block text file")
    group_required.add_argument('-o', '--output', type=str, help="output NEXUS file name")
    args = parser.parse_args()
    main()
