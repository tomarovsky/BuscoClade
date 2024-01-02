#!/usr/bin/env python3
__author__ = 'tomarovsky'
from ete3 import TextFace, Tree, faces, AttrFace, TreeStyle, NodeStyle
from argparse import ArgumentParser


def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"

def main():
    tree = "".join(open(args.input).readlines()).strip()
    if "'" in tree:
        tree = tree.replace("'", "")
    t = Tree(tree)
    if args.outgroup:
        outgroup = args.outgroup.replace(" ", "_")
        if "," in args.outgroup:
            try:
                nodes_to_root = outgroup.split(",")
                common_ancestor = t.get_common_ancestor(*nodes_to_root)
                t.set_outgroup(common_ancestor)
            except:
                R = t.get_midpoint_outgroup()
                t.set_outgroup(R)
                nodes_to_root = outgroup.split(",")
                common_ancestor = t.get_common_ancestor(*nodes_to_root)
                t.set_outgroup(common_ancestor)
        else:
            t.set_outgroup(outgroup)
    else:
        t.unroot()
    for i in t.get_leaves(): # 'Homo_sapiens' -> 'Homo sapiens'
        i.name = i.name.replace("_", " ")
    ts = TreeStyle()
    ts.mode = "r"
    ts.layout_fn = mylayout
    ts.show_leaf_name = False
    for n in t.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "Blue"
        nstyle["size"] = 0
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        n.set_style(nstyle)

    ts.show_branch_length = True
    ts.show_branch_support = True
    ts.branch_vertical_margin = -12
    t.render(f"{args.output}.length_and_support_tree.svg", w=500, units="px", tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -4
    t.render(f"{args.output}.only_support_tree.svg", w=500, units="px", tree_style=ts)
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = -2
    t.render(f"{args.output}.only_tree.svg", w=500, units="px", tree_style=ts)
    if args.show:
        t.show(tree_style=ts)

if __name__ == "__main__":
    parser = ArgumentParser(description="script to visualize phylogenetic trees using ete3 (required python < 3.10)")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="NEWICK file")
    group_required.add_argument('-o', '--output', type=str, help="outfile prefix")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-g', '--outgroup', type=str, default=False, help="outgroup species name (default = unrooted)")
    group_additional.add_argument('--show', action="store_true", help="option to show tree using GUI")
    args = parser.parse_args()
    main()

