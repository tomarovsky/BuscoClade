#!/usr/bin/env python3
__author__ = 'tomarovsky'
from ete3 import TextFace, Tree, faces, AttrFace, TreeStyle, NodeStyle, CircleFace
from argparse import ArgumentParser
import re

def format_name(name):
    name = re.sub(r'\.(GCA|GCF)_\d+\.\d+$', '', name)
    name = re.sub(r'[_\.]', ' ', name)
    return name.strip()

    
    return name

def mylayout(node):
    if node.is_leaf():
        N = AttrFace("name", fgcolor="black", text_prefix="  ", fstyle="italic", fsize=12)
        faces.add_face_to_node(N, node, 0)
        # F = AttrFace("support", fsize=8, fgcolor="green", text_prefix="")
        # faces.add_face_to_node(F, node, 0)
        node.img_style["size"] = 1
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "Black"
    else:
        support_value = node.support
        if support_value > 90:
            color = "LimeGreen"
        elif support_value > 70:
            color = "#008cf0"
        elif support_value > 50:
            color = "#883ac2"
        else:
            color = "#ff0000"
        
        support_circle = CircleFace(radius=3, color=color, style="circle")
        faces.add_face_to_node(support_circle, node, column=0, position="float")
        node.img_style["size"] = 0 
        
        

def export_legend(palette):
    legend = TreeStyle()
    legend.show_leaf_name = False
    legend.title.add_face(TextFace("LegendColors of \nnormalized values", fsize=10, bold=True), column=0)

    for label, color in palette.items():
        legend.legend.add_face(faces.CircleFace(6, color), column=0)
        legend.legend.add_face(TextFace(f"  {label}", fsize=8, fgcolor="black"), column=1)

    return legend


def main():
    palette = {"91 - 100": "LimeGreen", "71-90": '#008cf0', "51-70": '#883ac2', "<50": '#ff0000'}
    
    t = Tree(args.input)
    for i in t.get_leaves():
        i.name = format_name(i.name)
        #i.name = i.name.replace("_", " ").replace("GCA ", "GCA_").replace("GCF ", "GCF_")
    if args.outgroup:
        outgroup_species = args.outgroup.split(',')
        if len(outgroup_species) == 1:
            t.set_outgroup(outgroup_species[0])
        else:
            common_ancestor = t.get_common_ancestor(outgroup_species)
            t.set_outgroup(common_ancestor)
    else:
        t.unroot()

    ts = TreeStyle()
    ts.mode = "r"
    ts.layout_fn = mylayout
    ts.show_scale = True
    ts.show_leaf_name = False
    ts.legend_position = 1
    ts.legend = export_legend(palette).legend

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
    t.ladderize(direction=1)
    t.render(f"{args.output}.length_and_support_tree.png", w=3000, units="px", tree_style=ts)

    ts.show_branch_length = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = -4
    t.ladderize(direction=1)
    t.render(f"{args.output}.only_support_tree.png", w=3000, units="px", tree_style=ts)

    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = 0
    t.ladderize(direction=1)
    t.render(f"{args.output}.only_tree.png", w=6000, units="px", tree_style=ts)

    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.branch_vertical_margin = 0
    t.ladderize(direction=1)
    t.render(f"{args.output}.dots.png", w=6000, units="px", tree_style=ts)


    if args.show:
        t.show(tree_style=ts)


if __name__ == "__main__":
    parser = ArgumentParser(description="script to visualize phylogenetic trees using ete3 (required python < 3.10)")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="NEWICK file")
    group_required.add_argument('-o', '--output', type=str, help="outfile name")
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-g', '--outgroup', type=str, default=False, help="outgroup species name (default = unrooted)")
    group_additional.add_argument('--show', action="store_true", help="option to show tree using GUI")
    args = parser.parse_args()
    main()
