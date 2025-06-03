#!/usr/bin/env python3
import sys
from Bio import Phylo

def extract_species(name):
    return name.split("_")[-2] + "_" + name.split("_")[-1]  # last two parts

def main(input_tree, output_tree):
    tree = Phylo.read(input_tree, "newick")
    for tip in tree.get_terminals():
        tip.name = extract_species(tip.name)
    Phylo.write(tree, output_tree, "newick")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: rename_tree_tips_by_species.py input.tree output.tree")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
