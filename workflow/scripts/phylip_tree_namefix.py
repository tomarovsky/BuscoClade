#!/usr/bin/env python3
__author__ = 'tomarovsky'
import argparse

def main():
    with open(args.input, 'r') as file:
        nwk_tree = file.read().replace('\n', '')

    for species in args.species:
        species_prefix = species[:10]
        if species_prefix in nwk_tree:
            nwk_tree = nwk_tree.replace(species_prefix, species)
        else:
            print("Error! The first 10 characters of 2 or more species are identical.")
            print("Correct species names manually in <concat_alignment>.fna.phy and restart PHYLIP rules.")
            exit()

    with open(args.output, 'w') as nwk_namefix_tree:
        nwk_namefix_tree.write(nwk_tree)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix species names in NEWICK tree from PHYLIP")
    parser.add_argument('-i', '--input', type=str, required=True, help="input NEWICK tree")
    parser.add_argument('-s', '--species', type=str, nargs='+', required=True,
                        help="Comma-separated list of species names")
    parser.add_argument('-o', '--output', type=str, required=True, help="output NEWICK tree")
    args = parser.parse_args()
    main()
