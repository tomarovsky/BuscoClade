#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys
from argparse import ArgumentParser
from pathlib import Path

from ete3 import AttrFace, CircleFace, NodeStyle, TextFace, Tree, TreeStyle, faces


def newick_to_nhx(newick_file) -> str:
    with open(newick_file, 'r') as file:
        tree_string = ''
        newick = file.readline().replace("_", " ").strip().split("'")
        tree_string += newick[0]
        for i in range(1, len(newick), 2):
            line = ''
            flag = True
            for s in newick[i+1]:
                if s == ")" or s == ",":
                    if flag is True:
                        nhx = newick[i].replace(',', '.').replace(';', ':')[1:]
                        line += f"[&&NHX:{nhx}{s}"
                        flag = False
                    else:
                        line += s
                else:
                    line += s
            tree_string += line
        # print(tree_string)
        return tree_string


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"
        # node.dist = 0 # ASTRAL does not generate terminal branch lengths


def add_legend(ts):
    legend_items = [
        (" >90", "LimeGreen"),
        (" 71-90", "#008cf0"),
        (" 51-70", "#883ac2"),
        (" ≤50", "#ff0000"),
    ]

    for text, color in legend_items:
        circle = CircleFace(5, color=color, style="circle")
        label = TextFace(text, fsize=10, fgcolor="black")

        ts.legend.add_face(circle, column=0)
        ts.legend.add_face(label, column=1)

    ts.legend_position = (0, 0)

def process_tree(args):
    """Main logic for processing the tree."""
    input_path = Path(args.input)
    output_prefix = args.output if args.output else input_path.stem

    try:
        t = Tree(newick_to_nhx(args.input))
    except Exception as e:
        sys.exit(f"Error reading tree file: {e}")

    # 1. Rooting
    if args.outgroup:
        outgroup_names = [name.strip() for name in args.outgroup.split(',')]
        target_nodes = []
        for name in outgroup_names:
            node = t.search_nodes(name=name)
            if node:
                target_nodes.append(node[0])
            else:
                print(f"Warning: Outgroup species '{name}' not found. Skipping.")
        if target_nodes:
            if len(target_nodes) == 1:
                t.set_outgroup(target_nodes[0])
                print(f"Rooted using: {target_nodes[0].name}")
            else:
                mrca = t.get_common_ancestor(target_nodes)
                try:
                    t.set_outgroup(target_nodes)
                    outgroup_list = [node.name for node in target_nodes]
                    print(f"Rooted using group (MRCA of {len(outgroup_list)} species): {', '.join(outgroup_list)}")
                except Exception as e:
                    print(f"Error during set_outgroup: {e}. Attempting manual MRCA rooting.")
                    if mrca != t:
                        t.set_outgroup(mrca)
                        print(f"Rooted manually using MRCA of {len(target_nodes)} species.")
                    else:
                        print("Warning: MRCA is the current root. Cannot root, unrooting.")
                        t.unroot()

        else:
            print("Warning: No outgroup species found. Unrooting.")
            t.unroot()
    else:
        t.unroot()

    # 2. Ladderize (sort branches)
    t.ladderize(direction=True)

    # 3. Normalize leaf names
    for leaf in t.iter_leaves():
        leaf.name = (leaf.name
                    .replace("_", " ")
                    .replace("GCA ", "GCA_")
                    .replace("GCF ", "GCF_"))


    ts = TreeStyle()
    ts.mode = "r"
    ts.layout_fn = mylayout
    ts.show_leaf_name = False

    # Line style and metrics
    for n in t.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "Blue"
        nstyle["size"] = 0
        n.set_style(nstyle)
        if hasattr(n,"q1"):
            for metric in args.metrics:
                value = float(getattr(n, metric))
                normalized_value = value / args.number_of_genes * 100 if metric == 'EN' else value * 100
                for threshold in sorted(args.thresholds_and_colors.keys()):
                    if normalized_value >= threshold:
                        color = args.thresholds_and_colors[threshold] if metric in args.colored_metrics_whitelist else "Black"

                if len(metric) <= 2:
                    n.add_face(TextFace(f" {metric}  = "), column=1, position="branch-top")
                else:
                    n.add_face(TextFace(f" {metric} = "), column=1, position="branch-top")

                if metric == 'EN':
                    n.add_face(TextFace(f"{value:.2f} ({normalized_value:.2f}%) ", fgcolor = color), column=2, position="branch-top")
                else:
                    n.add_face(TextFace(f"{value:.0f} ", fgcolor=color), column=2, position="branch-top")

    add_legend(ts)

    # Render figure
    render_params = {
        "w": args.width,
        "units": "px",
        "tree_style": ts
    }
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = -4
    for f in args.output_formats:
        t.render(f"{args.output}.{f}", **render_params)

    if args.show:
        t.show(tree_style=ts)


def main():
    parser = ArgumentParser(description="script to visualize ASTRAL lll phylogenetic trees using ete3 (required python3 < 3.10)")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="NEWICK treefile from Astral lll with full annotation option (-t 2)")
    group_required.add_argument('-o', '--output', type=str, help="outfile prefix")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-g', '--outgroup', type=str, default=False, help="outgroup species name (default = unrooted)")

    # colorification:
    group_additional.add_argument('-m', '--metrics', type=lambda s: list(map(str, s.split(","))),
                    default=['q1', 'q2', 'pp1', 'pp2', 'EN'], help="comma-separated list of necessary metrics")
    group_additional.add_argument('--thresholds_and_colors', type=lambda s: dict(zip([int(s) for i in s[::2]], s[1::2])),
                    default={90: 'LimeGreen', 70: '#008cf0', 50: '#883ac2', 0: '#ff0000'}, help="colors per metrics"
                    "Example input: '90,LimeGreen,70,Gold,50,OrangeRed,0,Red'"
                    "This means that normalized values above 90 will be colored LimeGreen, values above 70 will be colored Gold, etc.")
    group_additional.add_argument('-n', '--number_of_genes', type=int,
                    help="total number of gene trees in ASTRAL input treefile (necessary to normalize 'EN' option value)")
    group_additional.add_argument('-w', '--colored_metrics_whitelist', type=lambda s: list(map(str, s.split(","))),
                    default=['EN'], help="comma-separated list of metrics for colorification (default metric color is 'Black')")

    # figure options:
    group_additional.add_argument('--width', type=int, default=800, help="width for result rendering")
    group_additional.add_argument('--show', action="store_true", help="option to show tree using GUI")
    group_additional.add_argument("-e", "--output_formats", dest="output_formats", type=lambda s: s.split(","),
                    default=("svg"), help="Comma-separated list of formats (supported by ete3) of output figure. Default: svg")
    args = parser.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"Error: File '{args.input}' not found.")

    process_tree(args)


if __name__ == "__main__":
    main()
