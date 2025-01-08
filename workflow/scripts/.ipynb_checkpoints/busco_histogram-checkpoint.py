#!/usr/bin/env python3
__author__ = "tomarovsky"
import argparse

import matplotlib.pyplot as plt
import pandas as pd


def main():
    data = pd.read_csv(args.input, sep="\t").sort_values(by=['M', 'F', 'D'], ascending=False).reset_index(drop=True)

    fig_height = len(data) * 0.3
    plt.figure(figsize=(13, fig_height), dpi=300)
    plt.subplots_adjust(top=0.8)

    position = range(len(data))
    bar_width = 0.8
    plt.barh(position, data["S"], height=bar_width, label="Complete and single-copy BUSCOs (S)", color=args.colors[0])
    plt.barh(position, data["D"], height=bar_width, left=data["S"], label="Complete and duplicated BUSCOs (D)", color=args.colors[1])
    plt.barh(position, data["F"], height=bar_width, left=data["S"] + data["D"], label="Fragmented BUSCOs (F)", color=args.colors[2])
    plt.barh(position, data["M"], height=bar_width, left=data["S"] + data["D"] + data["F"], label="Missing BUSCOs (M)", color=args.colors[3])

    plt.legend(ncol=4, loc=(-0.005, 0.97), handlelength=0.8, frameon=False)

    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)

    plt.yticks(range(len(data["Species"].to_list())), ["" for i in range(len(data["Species"].to_list()))])
    plt.xticks([0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%"])
    plt.xlim(0, 100)

    for i, row in data.iterrows():
        plt.text(-1, i, f'{row["Species"]}', va="center", ha="right", fontweight="medium", style="italic")
        plt.text(
            1,
            i,
            f"S: {row['S']}%   |   D: {row['D']}%   |   F: {row['F']}%   |   M: {row['M']}%",
            va="center",
            ha="left",
            fontweight="bold",
            color="white",
            fontsize=8,
        )

    plt.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BUSCO summaries visualization")
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-i", "--input", type=str, help="input TSV file from busco_summaries_to_tsv.py")
    group_required.add_argument("-o", "--output", type=str, help="output SVG file name")
    group_additional = parser.add_argument_group("Additional options")
    group_additional.add_argument(
        "-c",
        "--colors",
        type=lambda s: list(map(str, s.split(","))),
        default=["#23b4e8", "#008dbf", "#fbbc04", "#ea4335"],
        help="comma-separated list of colors per BUSCO metrics",
    )
    args = parser.parse_args()
    main()
