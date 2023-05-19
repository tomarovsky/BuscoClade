#!/usr/bin/env python
__author__ = 'tomarovsky'
import matplotlib
import matplotlib.pyplot as plt
from collections import Counter
plt.ioff()
import pandas as pd
import argparse


def unique_ids_plot(species_ids_countdict):
    df = pd.DataFrame({"Species": [len(v) for v in species_ids_countdict.values()]}, index=species_ids_countdict.keys())
    df.plot(kind="bar")
    plt.ylabel("Unique BUSCOs")
    plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
    plt.tight_layout()
    plt.savefig(args.outplot)


def main():
    species_dict = {}
    for file in args.species_ids_files:
        species_dict[file.split('/')[-1]] = []
        with open(file, 'r') as f:
            for line in f:
                species_dict[file.split('/')[-1]].append(line.rstrip())

    species_dict_values = []
    for v in species_dict.values():
        for i in v:
            species_dict_values.append(i)

    species_ids_countdict = Counter(species_dict_values)
    unique_ids_list = []
    for k, v in species_ids_countdict.items():
        print(k, v)
        if v == 1:
            unique_ids_list.append(k)

    species_unique_ids_dict = {}
    for k, v in species_dict.items():
        if k not in species_unique_ids_dict:
            species_unique_ids_dict[k] = []
        for u in unique_ids_list:
            if u in v:
                species_unique_ids_dict[k].append(u)

    df = pd.DataFrame.from_dict(species_unique_ids_dict, orient='index')
    df = df.transpose()
    df.to_csv(args.outcsv, sep='\t', index=False)
    unique_ids_plot(species_unique_ids_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to get plot and csv file of unique species ids")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-s', '--species_ids_files', type=str, nargs='+',
                                help="space separated list of species_ids_files")
    group_required.add_argument('--outplot', type=str, help="output plot file name")
    group_required.add_argument('--outcsv', type=str, help="output csv file name")
    args = parser.parse_args()
    main()