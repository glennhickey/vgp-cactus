#!/usr/bin/env python3
"""
Calculate phylogenetic distances between species in a Newick tree
"""

from Bio import Phylo
from io import StringIO
import sys
import argparse

def get_distance(tree, species1, species2):
    """Calculate the phylogenetic distance between two species"""
    return tree.distance(species1, species2)

def main():
    parser = argparse.ArgumentParser(
        description='Calculate phylogenetic distance between two species in a Newick tree'
    )
    parser.add_argument('--tree', '-t',
                        default='mammals.tre',
                        help='Path to Newick tree file (default: mammals.tre)')
    parser.add_argument('--species1', '-s1',
                        default='Homo_sapiens',
                        help='First species name (default: Homo_sapiens)')
    parser.add_argument('--species2', '-s2',
                        default='Mus_musculus',
                        help='Second species name (default: Mus_musculus)')

    args = parser.parse_args()

    # Read the tree file
    try:
        with open(args.tree, 'r') as f:
            tree_str = f.read()
    except FileNotFoundError:
        print(f"Error: Tree file '{args.tree}' not found", file=sys.stderr)
        sys.exit(1)

    # Parse the tree
    tree = Phylo.read(StringIO(tree_str), 'newick')
    print(args)

    # Calculate distance
    print(f"Calculating phylogenetic distance in {args.tree}:\n")

    try:
        distance = get_distance(tree, args.species1, args.species2)
        print(f"Distance between {args.species1} and {args.species2}: {distance:.6f}")
    except Exception as e:
        print(f"Error calculating distance: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
