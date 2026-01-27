#!/usr/bin/env python3
"""
Label ancestral nodes in a Newick tree using taxonomic information from a table.

Takes a tree and VGP table as input, labels internal nodes based on shared
taxonomic categories, and outputs the labeled tree.

Example usage with combined labels for major clades:

    python3 label-tree.py --tree input.nwk --table VGPPhase1-freeze-1.0.tsv \
        --combine "Sauropsida:Lineage:Birds,Reptiles" \
        --combine "Amniota:Lineage:Mammals,Birds,Reptiles" \
        --combine "Tetrapoda:Lineage:Mammals,Birds,Reptiles,Amphibians"
"""

import sys
import argparse
from Bio import Phylo
import pandas as pd
import io
import re

SEARCH_COLUMNS = [
    'Accession # for main haplotype',
    'Accession #s other high-quality haplotypes',
    'RefSeq annotation main haplotype',
    'UCSC Browser main haplotype',
]

def build_accession_lookup(table):
    """Build a dictionary mapping accessions to table rows for fast lookup."""
    lookup = {}
    for idx, row in table.iterrows():
        for col in SEARCH_COLUMNS:
            if col not in table.columns:
                continue
            cell_value = row[col]
            if pd.notna(cell_value):
                lookup[str(cell_value)] = row
    return lookup

def sanitize_label(label):
    """Convert a taxonomic label to a valid camelCase identifier."""
    if pd.isna(label) or not label:
        return None
    label = str(label).strip()
    # Split on non-alphanumeric characters
    words = re.split(r'[^a-zA-Z0-9]+', label)
    # Filter empty strings and capitalize each word (camelCase with first letter cap)
    words = [w.capitalize() for w in words if w]
    return ''.join(words) if words else None

def parse_hierarchical_value(val):
    """Parse a value like 'Rodentia (Erethizontidae)' into ['Rodentia', 'Erethizontidae'].
    Returns list from most general to most specific."""
    if pd.isna(val):
        return []
    val = str(val).strip()
    # Check for parenthesized subcategory
    match = re.match(r'^([^(]+)\s*\(([^)]+)\)\s*$', val)
    if match:
        base = match.group(1).strip()
        sub = match.group(2).strip()
        return [base, sub]
    return [val]

def compute_global_value_counts(tree, columns, accession_map, accession_lookup):
    """Compute the total count of leaves in the tree for each (column, level, value).

    Returns:
        dict mapping (col, level, value) -> count of leaves with that value
    """
    global_counts = {}
    for leaf in tree.get_terminals():
        accession = accession_map.get(leaf.name)
        if not accession:
            continue
        row = accession_lookup.get(accession)
        if row is None:
            continue
        for col in columns:
            if col in row.index:
                parts = parse_hierarchical_value(row[col])
                if len(parts) >= 1:
                    key = (col, 0, parts[0])
                    global_counts[key] = global_counts.get(key, 0) + 1
                if len(parts) >= 2:
                    key = (col, 1, parts[1])
                    global_counts[key] = global_counts.get(key, 0) + 1
    return global_counts


def get_shared_label(clade, columns, accession_map, accession_lookup, global_value_counts, combined_labels=None):
    """Find the most specific shared taxonomic label for all leaves under this clade.

    A label is only returned if this clade contains ALL leaves in the tree with that
    taxonomic value (i.e., this is the true clade root for that category).

    Args:
        clade: the clade to find a label for
        columns: list of column names to check
        accession_map: dict mapping leaf names to accessions
        accession_lookup: dict mapping accessions to table rows
        global_value_counts: dict mapping (col, level, value) -> total count in tree
        combined_labels: list of (label, column, value_set) tuples for combined categories
    """
    leaves = clade.get_terminals()

    # Collect values for each column from all leaves, splitting parenthesized values
    # Structure: {(col, level): set of values}
    # level 0 = base (e.g., "Rodentia"), level 1 = sub (e.g., "Erethizontidae")
    column_values = {}
    column_counts = {}  # Track how many leaves contributed to each level
    # Also track per-value counts: {(col, level, value): count}
    value_counts = {}
    for col in columns:
        column_values[(col, 0)] = set()  # base level
        column_values[(col, 1)] = set()  # sub level (parenthesized)
        column_counts[(col, 0)] = 0
        column_counts[(col, 1)] = 0

    total_leaves = 0
    for leaf in leaves:
        accession = accession_map.get(leaf.name)
        if not accession:
            continue
        row = accession_lookup.get(accession)
        if row is None:
            continue
        total_leaves += 1
        for col in columns:
            if col in row.index:
                parts = parse_hierarchical_value(row[col])
                if len(parts) >= 1:
                    column_values[(col, 0)].add(parts[0])
                    column_counts[(col, 0)] += 1
                    key = (col, 0, parts[0])
                    value_counts[key] = value_counts.get(key, 0) + 1
                if len(parts) >= 2:
                    column_values[(col, 1)].add(parts[1])
                    column_counts[(col, 1)] += 1
                    key = (col, 1, parts[1])
                    value_counts[key] = value_counts.get(key, 0) + 1

    # Find the most specific level where all leaves share the same value
    # AND this clade contains ALL leaves in the tree with that value
    # Check columns in reverse order (most specific first), and within each column
    # check sub-level before base-level
    # Only consider a level if ALL leaves contributed to it
    for col in reversed(columns):
        # Check sub-level first (more specific)
        values = column_values[(col, 1)]
        if len(values) == 1 and column_counts[(col, 1)] == total_leaves:
            value = list(values)[0]
            key = (col, 1, value)
            # Only use this label if we have ALL leaves with this value in the tree
            if value_counts.get(key, 0) == global_value_counts.get(key, 0):
                return sanitize_label(value)
        # Check base-level
        values = column_values[(col, 0)]
        if len(values) == 1 and column_counts[(col, 0)] == total_leaves:
            value = list(values)[0]
            key = (col, 0, value)
            # Only use this label if we have ALL leaves with this value in the tree
            if value_counts.get(key, 0) == global_value_counts.get(key, 0):
                return sanitize_label(value)

    # If no single shared value, check combined labels
    # Combined labels are sorted by size (smallest/most specific first)
    if combined_labels:
        for label, col, value_set in combined_labels:
            if col not in columns:
                continue
            # Check if all leaves have values within this combined set
            leaf_values = column_values.get((col, 0), set())
            if leaf_values and column_counts[(col, 0)] == total_leaves:
                if leaf_values.issubset(value_set) and len(leaf_values) > 1:
                    # Check that we have ALL leaves for each value in our set
                    has_all = True
                    for value in leaf_values:
                        key = (col, 0, value)
                        if value_counts.get(key, 0) != global_value_counts.get(key, 0):
                            has_all = False
                            break
                    if has_all:
                        return label

    return None

def label_tree(tree, table, columns, accession_map, fallback_label=None, combined_labels=None):
    """Label all ancestral nodes in the tree.

    Args:
        tree: Bio.Phylo tree object
        table: pandas DataFrame with taxonomic information
        columns: list of column names to use for labeling (general to specific)
        accession_map: dict mapping leaf names to accessions
        fallback_label: prefix for unlabeled ancestors (e.g., 'Vertebrates')
        combined_labels: list of (label, column, value_set) tuples for combined categories
    """
    accession_lookup = build_accession_lookup(table)

    # Pre-compute global value counts (how many leaves in tree have each value)
    global_value_counts = compute_global_value_counts(tree, columns, accession_map, accession_lookup)

    # Clear existing ancestor names
    for anc_clade in tree.get_nonterminals():
        anc_clade.confidence = None
        anc_clade.name = None

    # First pass: compute taxonomic labels for each node (bottom-up)
    node_labels = {}  # clade -> taxonomic label (or None)
    for anc_clade in tree.get_nonterminals(order='postorder'):
        label = get_shared_label(anc_clade, columns, accession_map, accession_lookup, global_value_counts, combined_labels)
        node_labels[anc_clade] = label

    # Second pass: count how many nodes will use each label (top-down)
    label_counts = {}  # label -> total count

    def count_labels(clade, inherited_label):
        own_label = node_labels.get(clade)
        if own_label:
            label_counts[own_label] = label_counts.get(own_label, 0) + 1
            current_label = own_label
        elif inherited_label:
            label_counts[inherited_label] = label_counts.get(inherited_label, 0) + 1
            current_label = inherited_label
        else:
            current_label = None
        for child in clade.clades:
            if not child.is_terminal():
                count_labels(child, current_label)

    count_labels(tree.root, None)

    # Calculate padding width for each label
    label_widths = {label: len(str(count - 1)) for label, count in label_counts.items()}

    # Third pass: assign names with AncX numbering (top-down)
    # Each clade with a label gets Anc0, descendants without own label get Anc1, Anc2, etc.
    label_counters = {}  # label -> next counter value

    def assign_names(clade, inherited_label):
        own_label = node_labels.get(clade)
        if own_label:
            # This node has its own taxonomic label - start fresh numbering
            if own_label not in label_counters:
                label_counters[own_label] = 0
            counter = label_counters[own_label]
            label_counters[own_label] += 1
            width = label_widths[own_label]
            clade.name = f"{own_label}Anc{counter:0{width}d}"
            current_label = own_label
        elif inherited_label:
            # Inherit parent's label with incremented counter
            counter = label_counters[inherited_label]
            label_counters[inherited_label] += 1
            width = label_widths[inherited_label]
            clade.name = f"{inherited_label}Anc{counter:0{width}d}"
            current_label = inherited_label
        else:
            current_label = None

        if clade.name:
            sys.stderr.write(f'Labeling ancestor: {clade.name} ({len(clade.get_terminals())} leaves)\n')

        # Recurse to children
        for child in clade.clades:
            if not child.is_terminal():
                assign_names(child, current_label)

    # Start from root
    assign_names(tree.root, None)

    # Fallback pass: label any remaining unlabeled ancestors
    unlabeled_nodes = [c for c in tree.get_nonterminals() if not c.name]
    if unlabeled_nodes and fallback_label:
        fallback_width = len(str(len(unlabeled_nodes) - 1)) if len(unlabeled_nodes) > 1 else 2
        fallback_counter = 0

        def assign_fallback_names(clade):
            nonlocal fallback_counter
            if not clade.is_terminal() and not clade.name:
                clade.name = f"{fallback_label}{fallback_counter:0{fallback_width}d}"
                sys.stderr.write(f'Fallback labeling ancestor: {clade.name} ({len(clade.get_terminals())} leaves)\n')
                fallback_counter += 1
            for child in clade.clades:
                assign_fallback_names(child)

        assign_fallback_names(tree.root)

def main(command_line=None):
    parser = argparse.ArgumentParser(
        description='Label ancestral nodes in a Newick tree using taxonomic information')
    parser.add_argument('--tree', required=True,
                        help='tree to read in newick format')
    parser.add_argument('--table', required=True,
                        help='VGP table with taxonomic information, ex tables/VGPPhase1-freeze-1.0.tsv')
    parser.add_argument('--columns', type=str,
                        default='Lineage,Superorder,Orders Scientific Name (inferred >50 MYA divergence times),Family Scientific Name',
                        help='Comma-separated list of columns to use for labeling (ordered general to specific)')
    parser.add_argument('--fallback-label', type=str, default='Vertebrates',
                        help='Fallback label prefix for ancestors that cannot be labeled from the table (default: Vertebrates)')
    parser.add_argument('--no-fallback', action='store_true',
                        help='Disable fallback labeling (leave unlabeled ancestors without names)')
    parser.add_argument('--accession-column', type=str, default=None,
                        help='Column in table to match leaf names against. If not specified, searches standard accession columns.')
    parser.add_argument('--combine', action='append', metavar='LABEL:COLUMN:VALUES',
                        help='Define a combined label for multiple category values. '
                             'Format: "Label:Column:Value1,Value2,Value3". '
                             'Example: "Amniota:Lineage:Mammals,Birds,Reptiles". '
                             'Can be specified multiple times.')

    args = parser.parse_args(command_line)

    # Read tree
    tree = Phylo.read(args.tree, "newick")

    # Validate tree is binary
    for anc_clade in tree.get_nonterminals():
        if len(anc_clade) != 2:
            sys.stderr.write('Error: input tree is not binary!\n')
            return 1

    # Read table
    table = pd.read_csv(args.table, sep='\t')

    # Build accession map (leaf name -> accession for lookup)
    # By default, assume leaf names ARE the accessions
    accession_map = {}
    for leaf_clade in tree.get_terminals():
        # If leaf name contains a suffix (e.g., "GCA_123-Homo_sapiens"), extract the accession
        accession = leaf_clade.name.split('-')[0]
        accession_map[leaf_clade.name] = accession

    # Parse columns
    columns = [col.strip() for col in args.columns.split(',')]

    # Parse combined labels
    combined_labels = []
    if args.combine:
        for combine_spec in args.combine:
            parts = combine_spec.split(':')
            if len(parts) != 3:
                sys.stderr.write(f'Error: --combine format must be "Label:Column:Value1,Value2,...", got "{combine_spec}"\n')
                return 1
            label, col, values_str = parts
            values = set(v.strip() for v in values_str.split(','))
            combined_labels.append((label, col, values))
            sys.stderr.write(f'Combined label: {label} = {col} in {{{", ".join(sorted(values))}}}\n')

        # Sort by set size (smallest first = most specific)
        combined_labels.sort(key=lambda x: len(x[2]))

    # Determine fallback label
    fallback_label = None if args.no_fallback else args.fallback_label

    # Label the tree
    label_tree(tree, table, columns, accession_map, fallback_label, combined_labels)

    # Output the labeled tree
    tree_stream = io.StringIO()
    Phylo.write(tree, tree_stream, 'newick')
    tree_stream.seek(0)
    tree_string = tree_stream.read().strip()

    # Remove any branch length from the root
    m = re.search(r':[0-9,.]*;', tree_string)
    if m:
        tree_string = tree_string.replace(m.group(), ';')

    print(tree_string)

    return 0

if __name__ == '__main__':
    sys.exit(main())
