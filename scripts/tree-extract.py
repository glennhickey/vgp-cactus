#!/usr/bin/env python3
"""
Use a table to extract a subtree
"""

import os, sys
import subprocess
import argparse
import shutil
from collections import defaultdict
from Bio import Phylo
import pandas as pd

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Add some node information to a tree from a table')
    parser.add_argument('--tree', required=True,
                        help='tree to read in newick format')
    parser.add_argument('--table', required=True,
                        help='table with first column being node name in the tree')
    parser.add_argument('--category', required=True,
                        help='column label in table we want to add to the node names')
    parser.add_argument('--value', required=True,
                        help='select species with this value for category in the table')
    parser.add_argument('--suffix-category',
                        help='add suffix to name using this column from the table')
    parser.add_argument('--outgroups', type=int, default=0,
                        help='number of outgroups to add')
    
    args = parser.parse_args(command_line)


    tree = Phylo.read(args.tree, "newick")
    tree_string = ''
    with open(args.tree) as tree_file:
        for line in tree_file:
            tree_string += line.strip()

    table = pd.read_csv(args.table, sep='\t', index_col='# accession')

    selection = []
    unselection = []
    
    for leaf_clade in tree.get_terminals():
        value = table.at[leaf_clade.name, args.category]
        if value == args.value:
            selection.append(leaf_clade)
        else:
            unselection.append(leaf_clade.name)

    anc = tree.common_ancestor(selection)

    subtree = Phylo.BaseTree.Tree().from_clade(anc)
    # this will fail if selected group is not monophyletic
    assert len(subtree.get_terminals()) == len(selection)

    # make the distance matrix for outgroups
    dists = []
    for og in unselection:
        dists.append((tree.distance(anc, og), og))
    dists = sorted(dists, key=lambda a : a[0])

    # remove the top N outgroups from the unselection
    for dist, og in dists[:args.outgroups]:        
        unselection.remove(og)
        sys.stderr.write(f'Selection outgroup {og} with distance {dist}\n')

    # prune the final tree
    if len(selection):
        for u in unselection:
            tree.prune(u)
        for leaf_clade in tree.get_terminals():
            # add a suffix from the table to make the names human readable
            if args.suffix_category:
                suffix = table.at[leaf_clade.name, args.suffix_category]
                suffix = suffix.replace(' ', '_')
                leaf_clade.name += '-' + suffix
        Phylo.write(tree, sys.stdout, 'newick')
    
if __name__ == '__main__':
    main()
