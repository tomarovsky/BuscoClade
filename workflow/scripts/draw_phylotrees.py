#!/usr/bin/env python3
__author__ = 'tomarovsky'

import sys
from argparse import ArgumentParser
from pathlib import Path

from ete3 import AttrFace, CircleFace, NodeStyle, TextFace, Tree, TreeStyle, faces


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, column=0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"
    else:
        if hasattr(node, "support"):
            if node.support > 90:
                color = "LimeGreen"
            elif node.support > 70:
                color = "#008cf0"
            elif node.support > 50:
                color = "#883ac2"
            else:
                color = "#ff0000"

            support_circle = CircleFace(4, color=color, style="circle")
            faces.add_face_to_node(support_circle, node, column=0, position="float")

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
        t = Tree(str(input_path))
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

    # 4. Set style
    ts = TreeStyle()
    ts.mode = "r" # Rectangular
    ts.layout_fn = mylayout
    ts.show_leaf_name = False

    # Общий стиль линий
    for n in t.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "Blue" # Line color
        nstyle["size"] = 0
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        n.set_style(nstyle)

    add_legend(ts)

    # 5. Render files
    # Define rendering parameters
    render_params = {
        "w": args.width,
        "units": "px",
        "tree_style": ts
    }

    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.branch_vertical_margin = -12
    t.render(f"{output_prefix}.length_and_support_tree.svg", **render_params)

    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -4
    t.render(f"{output_prefix}.only_support_tree.svg", **render_params)

    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = -2
    t.render(f"{output_prefix}.only_tree.svg", **render_params)

    print(f"Successfully generated 3 SVG files with prefix: {output_prefix}")

    if args.show:
        t.show(tree_style=ts)

def main():
    parser = ArgumentParser(description="Visualize phylogenetic trees using ete3 with custom styling.")

    req = parser.add_argument_group("Required options")
    req.add_argument("-i", "--input", type=str, required=True, help="Path to input NEWICK file")

    opt = parser.add_argument_group("Optional options")
    opt.add_argument("-o", "--output", type=str, help="Output file prefix (default: same as input filename)")
    opt.add_argument("-g", "--outgroup", type=str, help="Name of the outgroup species (if not found, uses midpoint rooting)")
    opt.add_argument("--width", type=int, default=800, help="Width of the image in pixels (default: 800)")
    opt.add_argument("--show", action="store_true", help="Open an interactive GUI window")

    args = parser.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"Error: File '{args.input}' not found.")

    process_tree(args)

if __name__ == "__main__":
    main()
