#!/usr/bin/env python3
__author__ = "tomarovsky"

import sys
from argparse import ArgumentParser
from pathlib import Path

from ete3 import AttrFace, CircleFace, NodeStyle, TextFace, Tree, TreeStyle, faces


def newick_to_nhx(newick_file) -> str:
    with open(newick_file, "r") as file:
        tree_string = ""
        newick = file.readline().strip().split("'")
        tree_string += newick[0]
        for i in range(1, len(newick), 2):
            line = ""
            flag = True
            for s in newick[i + 1]:
                if s == ")" or s == ",":
                    if flag:
                        nhx = newick[i].replace(",", ".").replace(";", ":")[1:]
                        line += f"[&&NHX:{nhx}{s}"
                        flag = False
                    else:
                        line += s
                else:
                    line += s
            tree_string += line
        return tree_string


def filter_tree(node, prefixes):
    """Filter tree based on prefix list"""
    if node.is_leaf():
        return any(node.name.startswith(prefix) for prefix in prefixes)
    keep_children = []
    for child in node.children:
        if filter_tree(child, prefixes):
            keep_children.append(child)
    node.children = keep_children
    return len(node.children) > 0


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"


def add_legend(ts):
    legend_items = [
        (" >0.9", "LimeGreen"),
        (" 0.71-0.9", "#008cf0"),
        (" 0.51-0.7", "#883ac2"),
        (" ≤0.5", "#ff0000"),
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

    # Rooting
    if args.outgroup:
        outgroup_names = [name.strip() for name in args.outgroup.split(",")]
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

    # Filter tree if prefixes provided
    if args.filter_prefixes:
        filter_tree(t, args.filter_prefixes)
        print(f"Tree filtered using prefixes: {args.filter_prefixes}")

    # Ladderize (sort branches)
    t.ladderize(direction=True)

    # Normalize leaf names
    for leaf in t.iter_leaves():
        leaf.name = leaf.name.replace("_", " ").replace("GCA ", "GCA_").replace("GCF ", "GCF_")

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
        if hasattr(n, "CULength"):
            for metric in args.metrics:
                value = float(getattr(n, metric))
                if metric in args.colored_metrics_whitelist:
                    color = "Black"
                    for threshold in sorted(args.thresholds_and_colors.keys()):
                        if value >= threshold:
                            color = args.thresholds_and_colors[threshold]
                else:
                    color = "Black"

                n.add_face(TextFace(f" {metric} ="), column=1, position="branch-top")
                n.add_face(TextFace(f"{value:.2f} ", fgcolor=color), column=2, position="branch-top")

    add_legend(ts)

    # Render figure
    render_params = {"dpi": args.dpi, "units": "px", "tree_style": ts}
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = -4
    for f in args.output_formats:
        t.render(f"{output_prefix}.{f}", **render_params)

    if args.show:
        t.show(tree_style=ts)


def main():
    parser = ArgumentParser(description="script to visualize ASTRAL-IV phylogenetic trees using ete3 (required python3 < 3.10)")
    group_required = parser.add_argument_group("Required options")
    group_required.add_argument("-i", "--input", type=str, help="extended NEWICK treefile (from Astral-IV with -u 2)")
    group_required.add_argument("-o", "--output", type=str, help="outfile prefix")
    group_additional = parser.add_argument_group("Additional options")
    group_additional.add_argument("-g", "--outgroup", type=str, default=False, help="outgroup species name (default = unrooted)")
    group_additional.add_argument("-f", "--filter_prefixes", type=str, default=False, help="comma-separated list of leaf name prefixes to filter out")

    # colorification:
    group_additional.add_argument(
        "-m",
        "--metrics",
        type=lambda s: list(map(str, s.split(","))),
        default=["CULength", "SULength", "f1", "f2", "f3", "localPP", "pp1", "pp2", "pp3", "q1", "q2", "q3"],
        help="comma-separated list of necessary metrics",
    )
    group_additional.add_argument(
        "--thresholds_and_colors",
        type=lambda s: dict(zip([float(s) for i in s[::2]], s[1::2])),
        default={0.9: "LimeGreen", 0.7: "#008cf0", 0.5: "#883ac2", 0: "#ff0000"},
        help="colors per metrics"
        "Example input: '0.9,LimeGreen,0.7,Gold,0.5,OrangeRed,0,Red'"
        "This means that values above 0.9 will be colored LimeGreen, values above 0.7 will be colored Gold, etc.",
    )
    group_additional.add_argument(
        "-w",
        "--colored_metrics_whitelist",
        type=lambda s: list(map(str, s.split(","))),
        default=["q1", "pp1"],
        help="comma-separated list of metrics for colorification (default metric color is 'Black')",
    )

    # figure options:
    group_additional.add_argument("--dpi", type=int, default=300, help="dpi for result rendering")
    group_additional.add_argument("--show", action="store_true", help="option to show tree using GUI")
    group_additional.add_argument(
        "-e",
        "--output_formats",
        dest="output_formats",
        type=lambda s: s.split(","),
        default=("svg"),
        help="Comma-separated list of formats (supported by ete3) of output figure. Default: svg",
    )
    args = parser.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"Error: File '{args.input}' not found.")

    process_tree(args)


if __name__ == "__main__":
    main()
