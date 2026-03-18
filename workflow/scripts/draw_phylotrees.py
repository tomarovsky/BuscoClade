#!/usr/bin/env python3
__author__ = "tomarovsky"

import sys
from argparse import ArgumentParser
from pathlib import Path

from ete3 import AttrFace, CircleFace, NodeStyle, TextFace, Tree, TreeStyle, faces

# Support color thresholds: (min_value, color)
SUPPORT_COLORS = [
    (90, "LimeGreen"),
    (70, "#008cf0"),
    (50, "#883ac2"),
    (0,  "#ff0000"),
]


def get_support_color(support):
    for threshold, color in SUPPORT_COLORS:
        if support > threshold:
            return color
    return SUPPORT_COLORS[-1][1]


def is_artificial_support_node(node):
    """Return True for the root and for nodes adjacent to root with a single-leaf outgroup."""
    if node.up is None:
        return True
    if node.up.up is None and any(child.is_leaf() for child in node.up.children):
        return True
    return False


def make_layout(show_support_text=False, show_branch_length=False):
    def layout(node):
        if node.is_leaf():
            face = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
            faces.add_face_to_node(face, node, column=0)
            node.img_style["size"] = 1
            node.img_style["shape"] = "circle"
            node.img_style["fgcolor"] = "black"
        else:
            if hasattr(node, "support") and not is_artificial_support_node(node):
                color = get_support_color(node.support)
                faces.add_face_to_node(CircleFace(3, color, style="circle"), node, column=0, position="float")
                if show_support_text:
                    faces.add_face_to_node(
                        TextFace(f" {node.support:.4g}", fsize=7, fgcolor="darkred"),
                        node, column=0, position="branch-top",
                    )

        if show_branch_length and node.up is not None:
            faces.add_face_to_node(
                TextFace(f" {node.dist:.4g}", fsize=7, fgcolor="grey"),
                node, column=0, position="branch-bottom",
            )

    return layout


def add_legend(ts):
    legend_items = [
        (" >90",   "LimeGreen"),
        (" 71-90", "#008cf0"),
        (" 51-70", "#883ac2"),
        (" ≤50",   "#ff0000"),
    ]
    for text, color in legend_items:
        ts.legend.add_face(CircleFace(5, color, style="circle"), column=0)
        ts.legend.add_face(TextFace(text, fsize=10, fgcolor="black"), column=1)
    ts.legend_position = (0, 0)


def normalize_phylip_support(t):
    """Rescale PHYLIP support values from 0-1 to 0-100."""
    for node in t.traverse():
        if not node.is_leaf() and hasattr(node, "support") and node.support is not None:
            node.support = round(node.support * 100, 2)


def root_tree(t, outgroup_str):
    """Root tree by outgroup. Supports single species or comma-separated list."""
    names = [name.strip() for name in outgroup_str.split(",")]
    target_nodes = []
    for name in names:
        found = t.search_nodes(name=name)
        if found:
            target_nodes.append(found[0])
        else:
            print(f"Warning: Outgroup '{name}' not found. Skipping.")

    if not target_nodes:
        print("Warning: No outgroup species found. Unrooting.")
        t.unroot()
        return

    if len(target_nodes) == 1:
        t.set_outgroup(target_nodes[0])
        print(f"Rooted using: {target_nodes[0].name}")
        return

    mrca = t.get_common_ancestor(target_nodes)
    try:
        t.set_outgroup(target_nodes)
        print(f"Rooted using group (MRCA of {len(target_nodes)} species): {', '.join(names)}")
    except Exception as e:
        print(f"Error during set_outgroup: {e}. Attempting manual MRCA rooting.")
        if mrca != t:
            t.set_outgroup(mrca)
            print(f"Rooted manually using MRCA of {len(target_nodes)} species.")
        else:
            print("Warning: MRCA is the current root. Unrooting.")
            t.unroot()


def normalize_leaf_names(t):
    for leaf in t.iter_leaves():
        leaf.name = (
            leaf.name
            .replace("'", "")
            .replace("_", " ")
            .replace("GCA ", "GCA_")
            .replace("GCF ", "GCF_")
        )


def build_tree_style():
    ts = TreeStyle()
    ts.mode = "r"
    ts.show_leaf_name = False
    ts.show_branch_support = False  # handled manually in layout
    ts.show_branch_length = False   # handled manually in layout
    return ts


def apply_node_style(t):
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "blue"
    nstyle["size"] = 0
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    for node in t.traverse():
        node.set_style(nstyle)


def render_trees(t, ts, output_prefix, dpi):
    render_params = {"dpi": dpi, "units": "px", "tree_style": ts}
    renders = [
        (f"{output_prefix}.length_and_support_tree.svg", -12, True,  True),
        (f"{output_prefix}.only_support_tree.svg",       -4,  True,  False),
        (f"{output_prefix}.only_tree.svg",               -2,  False, False),
    ]
    for filepath, margin, show_support, show_length in renders:
        ts.branch_vertical_margin = margin
        ts.layout_fn = make_layout(show_support_text=show_support, show_branch_length=show_length)
        t.render(filepath, **render_params)
    print(f"Successfully generated {len(renders)} SVG files with prefix: {output_prefix}")


def process_tree(args):
    input_path = Path(args.input)
    output_prefix = args.output if args.output else input_path.stem

    try:
        t = Tree(str(input_path))
    except Exception as e:
        sys.exit(f"Error reading tree file: {e}")

    if args.phylip_support:
        normalize_phylip_support(t)
        print("PHYLIP support values normalized: multiplied by 100 (0-1 -> 0-100 scale).")

    if args.outgroup:
        root_tree(t, args.outgroup)
    else:
        t.unroot()

    t.ladderize(direction=True)
    normalize_leaf_names(t)

    ts = build_tree_style()
    apply_node_style(t)
    add_legend(ts)

    render_trees(t, ts, output_prefix, args.dpi)

    if args.show:
        t.show(tree_style=ts)


def main():
    parser = ArgumentParser(description="Visualize phylogenetic trees using ete3 with custom styling.")

    req = parser.add_argument_group("Required options")
    req.add_argument("-i", "--input", required=True, help="Path to input Newick file")

    opt = parser.add_argument_group("Optional options")
    opt.add_argument("-o", "--output", help="Output file prefix (default: input filename stem)")
    opt.add_argument("-g", "--outgroup", help="Outgroup species name(s), comma-separated")
    opt.add_argument("--dpi", type=int, default=300, help="DPI for rendering (default: 300)")
    opt.add_argument("--show", action="store_true", help="Open an interactive GUI window")
    opt.add_argument("--phylip-support", action="store_true",
                     help="Rescale PHYLIP bootstrap support values from 0-1 to 0-100 scale")

    args = parser.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"Error: File '{args.input}' not found.")

    process_tree(args)


if __name__ == "__main__":
    main()
