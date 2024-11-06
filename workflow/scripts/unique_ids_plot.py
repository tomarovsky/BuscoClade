#!/usr/bin/env python
__author__ = "tomarovsky"
import os
import matplotlib.pyplot as plt
from supervenn import supervenn

plt.ioff()
import argparse


def read_species_ids(file_path):
    with open(file_path, 'r') as file:
        return {line.strip() for line in file}


def main():
    sets = [read_species_ids(file) for file in args.species_ids_files]
    labels = [os.path.splitext(os.path.basename(file))[0] for file in args.species_ids_files]
    fig, ax = plt.subplots(1, 1, figsize=(20, 10), dpi=300)

    supervenn(sets,
              labels,
              ax=ax,
              sets_ordering='minimize gaps',
              chunks_ordering='minimize gaps',
              )

    plt.savefig(args.outplot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to get plot and csv file of unique species ids")
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-s", "--species_ids_files", type=str, nargs="+", help="space separated list of species_ids_files")
    group_required.add_argument("--outplot", type=str, help="output plot file name")
    args = parser.parse_args()
    main()
