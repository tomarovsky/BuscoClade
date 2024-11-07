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
    labels = [os.path.splitext(os.path.basename(file))[0] for file in args.single_copy_ids_files]
    single_copy_sets = [read_species_ids(file) for file in args.single_copy_ids_files]
    multi_copy_sets = [read_species_ids(file) for file in args.multi_copy_ids_files]
    fig, ax = plt.subplots(2, 1, figsize=(30, len(labels) * 2), dpi=300)

    supervenn(single_copy_sets,
              labels,
              ax=ax[0],
              sets_ordering='size',
              chunks_ordering='size',
              min_width_for_annotation=50,
              rotate_col_annotations=True,
              col_annotations_area_height=1.2,
              )

    supervenn(multi_copy_sets,
              labels,
              ax=ax[1],
              sets_ordering='size',
              chunks_ordering='size',
              min_width_for_annotation=10,
              rotate_col_annotations=True,
              col_annotations_area_height=1.2,
              )

    ax[0].set_title("Single copy BUSCOs", fontsize=16, fontweight='bold')
    ax[1].set_title("Multi copy BUSCOs", fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.savefig(args.outplot)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to get plot and csv file of unique species ids")
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-s", "--single_copy_ids_files", type=str, nargs="+", help="space separated list of single copy_ids_files")
    group_required.add_argument("-m", "--multi_copy_ids_files", type=str, nargs="+", help="space separated list of multi copy ids_files")
    group_required.add_argument("--outplot", type=str, help="output plot file name")
    args = parser.parse_args()
    main()
