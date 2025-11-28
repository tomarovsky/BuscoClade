#!/usr/bin/env python3
__author__ = 'tomarovsky'

import argparse
import os
import sys
import tempfile

import matplotlib.pyplot as plt
from ete3 import CircleFace, ImgFace, NodeStyle, TextFace, Tree, TreeStyle, faces


def newick_to_nhx(newick_file) -> str:
    """Convert Newick format to NHX format"""
    with open(newick_file, "r") as file:
        tree_string = ""
        newick = file.readline().replace("_", " ").strip().split("'")
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


def add_legend(ts):
    """Add legend to tree style"""
    legend_items = [
        (" 100", "blue"),
        (" 90-99", "#009022"),
        (" 70-89", "#e6c700"),
        (" 0-69", "#e60000"),
    ]

    # empty_circle = CircleFace(5, color="#000000", style="circle")
    # ts.legend.add_face(empty_circle, column=0)
    ts.legend.add_face(TextFace(" ", fsize=14), column=0)
    ts.legend.add_face(TextFace(" Main local posterior probability (pp1):", fsize=10), column=1)

    for text, color in legend_items:
        circle = CircleFace(5, color=color, style="circle")
        label = TextFace(text, fsize=10, fgcolor="black")

        ts.legend.add_face(circle, column=0)
        ts.legend.add_face(label, column=1)

    ts.legend_position = (0, 0)



def create_tree_style():
    """Create and configure tree style"""
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.scale = 40
    ts.legend = TreeStyle().legend
    add_legend(ts)
    return ts


def main():
    parser = argparse.ArgumentParser(description='Visualize phylogenetic trees with ETE3')
    parser.add_argument('-i', '--input', required=True, help='Input tree file (Newick format)')
    parser.add_argument('-o', '--output', required=True, help='Output file prefix')
    parser.add_argument('-g', '--outgroup', default=False, help="outgroup species name (default = unrooted)")
    parser.add_argument('--ladderize', action='store_true', default=True, help='Ladderize tree (default: True)')
    parser.add_argument('--format', default='svg', help='Output format (svg, png, pdf, etc.)')
    parser.add_argument('--dpi', type=int, default=300, help='Output DPI for raster formats')
    parser.add_argument('--filter-prefixes', nargs='+', help='Filter leaves by prefixes')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)

    try:
        # Load and process tree
        tree = Tree(newick_to_nhx(args.input))

        # Set outgroup
        if args.outgroup:
            outgroup_names = [name.strip() for name in args.outgroup.split(',')]
            target_nodes = []
            for name in outgroup_names:
                node = tree.search_nodes(name=name)
                if node:
                    target_nodes.append(node[0])
                else:
                    print(f"Warning: Outgroup species '{name}' not found. Skipping.")
            if target_nodes:
                if len(target_nodes) == 1:
                    tree.set_outgroup(target_nodes[0])
                    print(f"Rooted using: {target_nodes[0].name}")
                else:
                    mrca = tree.get_common_ancestor(target_nodes)
                    try:
                        tree.set_outgroup(target_nodes)
                        outgroup_list = [node.name for node in target_nodes]
                        print(f"Rooted using group (MRCA of {len(outgroup_list)} species): {', '.join(outgroup_list)}")
                    except Exception as e:
                        print(f"Error during set_outgroup: {e}. Attempting manual MRCA rooting.")
                        if mrca != tree:
                            tree.set_outgroup(mrca)
                            print(f"Rooted manually using MRCA of {len(target_nodes)} species.")
                        else:
                            print("Warning: MRCA is the current root. Cannot root, unrooting.")
                            tree.unroot()

            else:
                print("Warning: No outgroup species found. Unrooting.")
                tree.unroot()
        else:
            tree.unroot()

        # Filter tree if prefixes provided
        if args.filter_prefixes:
            filter_tree(tree, args.filter_prefixes)

        # Ladderize if requested
        if args.ladderize:
            tree.ladderize(direction=True)

        # Node counter for numbering
        node_counter = {"count": 1}

        def my_layout(node):
            """Custom layout function"""
            node.img_style["size"] = 0
            node.img_style["fgcolor"] = "black"

            if node.is_leaf():
                face_text = " " + node.name
                face = TextFace(face_text, fsize=17, fstyle="italic", fgcolor="black")
                faces.add_face_to_node(face, node, column=0)
            else:
                # Support circles for pp1
                circle_size = 5
                if hasattr(node, "pp1"):
                    try:
                        pp1 = float(node.pp1)
                        if pp1 == 1.0:
                            support_color = "blue"
                        elif pp1 >= 0.9:
                            support_color = "#009022"
                        elif pp1 >= 0.7:
                            support_color = "#e6c700"
                        else:
                            support_color = "#e60000"
                        circle_face = CircleFace(radius=circle_size, color=support_color, style="circle")
                        faces.add_face_to_node(circle_face, node, column=1, position="float-behind")
                    except Exception as e:
                        print(f"Warning: Error processing pp1: {e}", file=sys.stderr)

                # Node numbering
                number = node_counter["count"]
                node.node_number = number
                number_face = TextFace(f" {number}", fsize=16, fgcolor="black")
                faces.add_face_to_node(number_face, node, column=0, position="float-behind")
                node_counter["count"] += 1

        # Create tree style and set layout
        ts = create_tree_style()
        ts.layout_fn = my_layout

        # Apply node styles
        nstyle = NodeStyle()
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_width"] = 1
        for n in tree.traverse():
            n.set_style(nstyle)

        # Render tree
        output_file = f"{args.output}.{args.format}"
        tree.render(output_file, tree_style=ts, dpi=args.dpi)
        print(f"Tree successfully rendered to: {output_file}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
