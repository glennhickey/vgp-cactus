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

SEARCH_COLUMNS = [
    'Accession # for main haplotype',
    'Accession #s other high-quality haplotypes',
    'RefSeq annotation main haplotype',
    'UCSC Browser main haplotype',
]

def find_row_by_name(table, name):
    """Search specific columns in the table for a matching value.

    Returns the first row where any of the SEARCH_COLUMNS contains the given name.
    Raises KeyError if the name is not found.
    """
    for idx, row in table.iterrows():
        for col in SEARCH_COLUMNS:
            if col not in table.columns:
                continue
            cell_value = row[col]
            if pd.notna(cell_value) and str(cell_value) == name:
                return row
    raise KeyError(f"Name '{name}' not found in columns {SEARCH_COLUMNS}")

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Add some node information to a tree from a table')
    parser.add_argument('--tree', required=True,
                        help='tree to read in newick format')
    parser.add_argument('--table',
                        help='Big vgp table for sex chromosome information, ex tables/VGPPhase1-freeze-1.0.tsv')
    parser.add_argument('--category', required=True,
                        help='column label in table we want to add to the node names')
    parser.add_argument('--value', required=True, action='append',
                        help='select species with this value for category in the table (can be specified multiple times)')
    parser.add_argument('--ignore', type=str,
                        help='ignore species with this value for category when selecting outgroups')
    parser.add_argument('--suffix-category',
                        help='add suffix to name using this column from the table')
    parser.add_argument('--outgroups', type=int, default=0,
                        help='number of outgroups to add')
    parser.add_argument('--prune-excluded', action='store_true',
                        help='remove elements that share LCA with selection but don\'t match category/value')
    parser.add_argument('--rename-column', default='UCSC Browser main haplotype',
                        help='replace node names with value from this column (use empty string to disable)')

    args = parser.parse_args(command_line)

    # Ensure ignore value is not in the selection values
    if args.ignore:
        assert(args.ignore not in args.value), f"Ignore value '{args.ignore}' cannot be in selection values {args.value}"

    tree = Phylo.read(args.tree, "newick")
    tree_string = ''
    with open(args.tree) as tree_file:
        for line in tree_file:
            tree_string += line.strip()

    table = pd.read_csv(args.table, sep='\t')

    selection = []
    unselection = []
    pruned_ingroups = []
    ignore = []
    
    for leaf_clade in tree.get_terminals():
        row = find_row_by_name(table, leaf_clade.name)
        value = row[args.category]
        if value in args.value:
            selection.append(leaf_clade.name)
        else:
            unselection.append(leaf_clade.name)
        if args.ignore and args.ignore == value:
            ignore.append(leaf_clade.name)

    anc = tree.common_ancestor(selection)

    subtree = Phylo.BaseTree.Tree().from_clade(anc)
    assert len(subtree.get_terminals()) >= len(selection)
    if len(subtree.get_terminals()) > len(selection):
        # our subtree contains elements we didn't select
        # so presumably the lineage chosen in not mono-phyletic?
        # we report (but keep) the unspecified elements:
        selection_set = set(selection)
        unselection_set = set(unselection)
        for sub_leaf in subtree.get_terminals():
            if sub_leaf.name not in selection_set:
                assert sub_leaf.name in unselection_set
                row = find_row_by_name(table, sub_leaf.name)
                value = row[args.category]
                assert value not in args.value
                if not args.prune_excluded:
                    sys.stderr.write(f'Warning: genome {sub_leaf.name} with {args.category}={value} will be included in tree\n')
                    unselection_set.remove(sub_leaf.name)
                    selection_set.add(sub_leaf.name)
                else:
                    sys.stderr.write(f'Warning: genome {sub_leaf.name} with {args.category}={value} will be pruned from tree\n')
                    pruned_ingroups.append(sub_leaf.name)
        selection = list(selection_set)
        unselection = list(unselection_set)
                        

    if args.outgroups:
        # make the distance matrix for outgroup candidates, and clade root
        # this is slow, and only used for multiple outgroups
        og_candidates = set(unselection) - set(pruned_ingroups)
        og_candidates -= set(ignore)
        if args.outgroups > 1 and len(og_candidates) > 100:
            sys.stderr.write(f'Computing brute force distance matrix (N={len(og_candidates)}), this may take a while ....\n')
        matrix = defaultdict(dict)
        for og in og_candidates:
            matrix[og][anc] = tree.distance(og, anc)
            matrix[anc][og] = matrix[og][anc]
            if args.outgroups > 1:
                for og2 in og_candidates:
                    if og2 in matrix and og in matrix[og2]:
                        matrix[og][og2] = matrix[og2][og]
                    else:
                        matrix[og][og2] = tree.distance(og, og2)
                    
        # choose outgroups that are closest to anc, but not similar to each other
        # this is a very simple brute-force loop, but shold be okay for vgp tree
        og_selection = set()
        while len(og_selection) < args.outgroups and len(og_candidates) > 0:
            min_dist = sys.float_info.max
            min_og = None
            for og in og_candidates:
                dist_anc = matrix[og][anc]
                dist_og = min([matrix[og][og2] for og2 in og_selection]) if len(og_selection) else sys.float_info.max
                # simple heuristic: split weight over total distance and distance to outgroup set
                # todo: something more sophisticated?
                penalty = (1.0 - dist_og / dist_anc) * dist_anc if dist_og < dist_anc else 0.
                if dist_anc + penalty < min_dist:
                    min_dist = dist_anc + penalty
                    min_og = og
            og_selection.add(min_og)
            og_candidates.remove(min_og)
            sys.stderr.write(f'Selecting outgroup {min_og} with distance {min_dist}\n')

        for og in og_selection:
            unselection.remove(og)

    # prune the final tree
    if len(selection):
        for u in unselection:
            tree.prune(u)
        for leaf_clade in tree.get_terminals():
            row = find_row_by_name(table, leaf_clade.name)
            # optionally rename the node using a column from the table
            if args.rename_column:
                new_name = row[args.rename_column]
                if pd.notna(new_name) and str(new_name) != leaf_clade.name:
                    sys.stderr.write(f'Warning: renaming {leaf_clade.name} to {new_name}\n')
                    leaf_clade.name = str(new_name)
            # add a suffix from the table to make the names human readable
            if args.suffix_category:
                suffix = row[args.suffix_category]
                suffix = suffix.replace(' ', '_').rstrip('_')
                leaf_clade.name += '-' + suffix
        Phylo.write(tree, sys.stdout, 'newick')
    
if __name__ == '__main__':
    main()
