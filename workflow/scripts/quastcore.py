#!/usr/bin/env python3
__author__ = "tomarovsky"

import bz2
import gzip
from argparse import ArgumentParser
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd


def metaopen(filename, mode="rt", buffering=None):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    elif filename.endswith(".bz2"):
        return bz2.open(filename, mode)
    return open(filename, mode, buffering=buffering) if buffering else open(filename, mode)


def parse_fasta(path, buffering=None):
    print("parse_sequences started")
    seqs = OrderedDict()
    lengths_dict = {}
    header = None

    with metaopen(path, "rt", buffering) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].split(" ")[0]
                seqs[header] = []
            else:
                seqs[header].append(line.rstrip())

    for name, chunks in seqs.items():
        seqs[name] = "".join(chunks)
        lengths_dict[name] = len(seqs[name])

    lengths_df = pd.DataFrame.from_dict({"lengths": lengths_dict})
    return seqs, lengths_df


def contig_count(df, min_contig) -> int:
    print("contig_count started")
    return len(df[df["lengths"] >= min_contig].index)


def largest_contig_length(df) -> int:
    print("largest_contig_length started")
    return df["lengths"].max()


def total_length(df, min_contig) -> int:
    print("total_length started")
    return df[df["lengths"] >= min_contig]["lengths"].sum()


def n_amount(df, sequences_dict, min_contig):
    print("n_amount started")
    count = 0
    index_list = df[df["lengths"] >= min_contig].index
    for i in index_list:
        count += sequences_dict[i].upper().count("N")
    return count


def gc_content(df, sequences_dict, min_contig) -> float:
    print("gc_content started")
    count = 0
    index_list = df[df["lengths"] >= min_contig].index
    for i in index_list:
        count += sequences_dict[i].upper().count("G") + sequences_dict[i].upper().count("C")

    t_len = total_length(df, min_contig)
    result = round(count / t_len * 100, 2) if t_len > 0 else 0
    return result


def n_l_statistics(df, percent: int, min_contig: int) -> list:
    print("n_l_statistics started")
    lengths = df[df["lengths"] >= min_contig]["lengths"].to_numpy(dtype=np.int64)
    if lengths.size == 0:
        return [None, None]

    lengths = np.sort(lengths)[::-1]
    total = lengths.sum()
    target = total * (percent / 100.0)

    csum = np.cumsum(lengths)
    idx = int(np.searchsorted(csum, target, side="left"))

    n_xx = int(lengths[idx])
    l_xx = idx + 1
    return [n_xx, l_xx]


def main():
    parser = ArgumentParser(description="This program annotates sequence.fasta. Provides various features.")

    group_required = parser.add_argument_group("Required options")
    group_required.add_argument(
        "-i", "--input", type=str, nargs="+", required=True, help="Input files. One or more files to run the program. Supports .gz, .bz2"
    )

    group_additional = parser.add_argument_group("Additional options")
    group_additional.add_argument("-b", "--buffering", metavar="INT", type=int, default=10000000, help="Text buffering. Default = 10000000")
    group_additional.add_argument(
        "-m",
        "--min-contig",
        metavar="INT",
        type=int,
        nargs="+",
        default=[0, 150, 500, 1000],
        help="Minimum contig length thresholds. Default: 0 150 500 1000",
    )
    group_additional.add_argument(
        "-n", "--nl-statistics", metavar="INT", type=int, nargs="+", default=[50, 75], help="Values for N and L statistics. Default: 50 75"
    )
    group_additional.add_argument("-o", "--output", type=str, default="quast_stats.csv", help="Output file name. Default: quast_stats.csv")
    group_additional.add_argument(
        "-p", "--print", type=lambda x: str(x).lower() in ["true", "yes", "1"], default=True, help="Print to terminal? Default: TRUE"
    )

    args = parser.parse_args()
    data_dict = defaultdict(list)

    for file in args.input:
        seq_dict, lengths_df = parse_fasta(file, args.buffering)

        for mc in args.min_contig:
            data_dict[f"Contigs(>={mc})"].append(contig_count(lengths_df, mc))

        for mc in args.min_contig:
            data_dict[f"Totallength(>={mc})"].append(total_length(lengths_df, mc))

        for mc in args.min_contig:
            data_dict[f"GC(%)(>={mc})"].append(gc_content(lengths_df, seq_dict, mc))

        for mc in args.min_contig:
            data_dict[f"N_amount(>={mc})"].append(n_amount(lengths_df, seq_dict, mc))

        data_dict["Largestcontig"].append(largest_contig_length(lengths_df))

        for p in args.nl_statistics:
            for mc in args.min_contig:
                n_l_res = n_l_statistics(lengths_df, p, mc)
                data_dict[f"N{p}stats(>={mc})"].append(n_l_res[0])
                data_dict[f"L{p}stats(>={mc})"].append(n_l_res[1])

    df_final = pd.DataFrame(data_dict, index=args.input)

    if args.print:
        print(df_final.T.astype(str))

    if args.output:
        df_final.to_csv(args.output, encoding="utf-8")


if __name__ == "__main__":
    main()
