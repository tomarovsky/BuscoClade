#!/usr/bin/env python
__author__ = "tomarovsky"
import os

import matplotlib.pyplot as plt
from supervenn import supervenn

plt.ioff()
import argparse


def read_species_ids(file_path):
    with open(file_path, "r") as file:
        return {line.strip() for line in file if line.strip()}


def plot_supervenn(sets, labels, ax, title, min_width_for_annotation, widths_minmax_ratio=None):
    # Filter out empty sets
    filtered = [(s, l) for s, l in zip(sets, labels) if s]
    if not filtered:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=14, transform=ax.transAxes)
        ax.set_title(title, fontsize=16, fontweight="bold")
        return
    filtered_sets, filtered_labels = zip(*filtered)
    kwargs = dict(
        sets_ordering=None,
        chunks_ordering="size",
        min_width_for_annotation=min_width_for_annotation,
        rotate_col_annotations=True,
        col_annotations_area_height=1.4,
    )
    if widths_minmax_ratio is not None:
        kwargs["widths_minmax_ratio"] = widths_minmax_ratio
    supervenn(list(filtered_sets), list(filtered_labels), ax=ax, **kwargs)
    ax.set_title(title, fontsize=16, fontweight="bold")


def main():
    labels = [os.path.splitext(os.path.basename(file))[0] for file in args.single_copy_ids_files]
    single_copy_sets = [read_species_ids(file) for file in args.single_copy_ids_files]
    multi_copy_sets = [read_species_ids(file) for file in args.multi_copy_ids_files]
    fig, ax = plt.subplots(2, 1, figsize=(30, len(labels)), dpi=300)
    plot_supervenn(single_copy_sets, labels, ax[0], "Single copy BUSCOs", min_width_for_annotation=50)
    plot_supervenn(multi_copy_sets, labels, ax[1], "Multi copy BUSCOs", min_width_for_annotation=10, widths_minmax_ratio=0.005)
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
