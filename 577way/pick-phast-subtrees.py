#!/usr/bin/env python3
"""Print the named-ancestor lineage (with leaf counts) of one or more leaf
accessions in the 577-way tree. Use the output to hand-pick `--subtree`
arguments for cactus-phast: aim for two meaningful clades per reference
(one broader, one more local) that are not too small.

    ./pick-phast-subtrees.py vgp-577way.nwk GCF_037176765.1 GCF_011100685.1
"""
import sys
from Bio import Phylo


def main(tree_path, accessions):
    tree = Phylo.read(tree_path, 'newick')
    parents = {c: p for p in tree.find_clades() for c in p.clades}
    sizes = {c.name: len(c.get_terminals()) for c in tree.find_clades() if c.name}
    leaves_by_name = {l.name: l for l in tree.get_terminals()}

    for acc in accessions:
        leaf = leaves_by_name.get(acc)
        if leaf is None:
            print(f"{acc}: not found in tree", file=sys.stderr)
            continue
        print(f"\n{acc}:")
        node = leaf
        while node in parents:
            node = parents[node]
            if node.name:
                print(f"  {node.name}\t{sizes[node.name]} leaves")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2:])
