#!/usr/bin/env python3
__author__ = 'tomarovsky'
from Bio import SeqIO
import argparse


def main():
    records = SeqIO.parse(args.input, "fasta")
    count = SeqIO.write(records, args.output, "stockholm")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for converting FASTA format to Stockholm format")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="input concat FASTA file or stdin")
    group_required.add_argument('-o', '--output', type=str, help="output NEXUS file name")
    args = parser.parse_args()
    main()
