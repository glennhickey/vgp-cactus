#!/usr/bin/env python3
"""
Make a Cactus Seqfile from a Newick (sub)tree using the vgp data tables

Also make the --chromInfo cactus input using the big vgp table

Note: For labeling ancestral nodes, use the separate label-tree.py script.
"""

import sys
import argparse
from Bio import Phylo
import pandas as pd
import io
import re

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Download some tables necessary to make the Cactus input')
    parser.add_argument('--tree', required=True,
                        help='tree to read in newick format')
    parser.add_argument('--urls', required=True,
                        help='URL download table, ex tables/URL.download.table.tsv')
    parser.add_argument('--table',
                        help='Big vgp table for sex chromosome information, ex tables/VGPPhase1-freeze-1.0.tsv')
    parser.add_argument('--chrom-info', type=str,
                        help='Write sex chrom info to this file (requires --table)')
    parser.add_argument('--suffix-only', action='store_true',
                        help='Only keep suffix added by tree-extract.py --suffix-category')

    args = parser.parse_args(command_line)

    tree = Phylo.read(args.tree, "newick")

    # strip out confidence values and numeric ancestor names (but preserve labels from label-tree.py)
    for anc_clade in tree.get_nonterminals():
        anc_clade.confidence = None
        # Only clear numeric names (e.g., "100", "0.95") - preserve text labels
        if anc_clade.name and anc_clade.name.replace('.', '').replace('-', '').isdigit():
            anc_clade.name = None
        if len(anc_clade) != 2:
            sys.stderr.write('Error: input tree is not binary!\n')
            return 1

    # strip accessions but remember them
    name_to_accession = {}
    for leaf_clade in tree.get_terminals():
        accession = leaf_clade.name.split('-')[0]
        name = leaf_clade.name.split('-')[-1] if args.suffix_only else leaf_clade.name
        name_to_accession[name] = accession
        leaf_clade.name = name

    # make sure the ancestor doesn't have a distance
    tree_stream = io.StringIO()
    Phylo.write(tree, tree_stream, 'newick')
    tree_stream.seek(0)
    tree_string = tree_stream.read().strip()
    m = re.search(f':[0-9,.]*;', tree_string)
    tree_string = tree_string.replace(m.group(), ';')
    print(tree_string)
    print('')
    
    url_table = pd.read_csv(args.urls, sep='\t', index_col='# accession')

    for name, accession in name_to_accession.items():
        # in case we added a suffix
        url = url_table.at[accession, 'download URL']
        print(f'{name}\t{url}')

    if args.chrom_info:
        if not args.chrom_info:
            sys.stderr.write('need --table with --chrom-info')
            return 1        
        table = pd.read_csv(args.table, sep='\t', index_col='UCSC Browser main haplotype')
        with open(args.chrom_info, 'w') as chrom_file:
            for name, accession in name_to_accession.items():
                try:
                    sex_chroms = table.at[accession, 'Sex chromosomes main haploptype']
                except:
                    sys.stderr.write(f'Warning: {accession} not found in table\n')
                chrom_list = []
                for tok in str(sex_chroms).replace(',', '').split():
                    if len(tok) in [1,2] and tok[0].isupper() and (len(tok) == 1 or tok[1].isnumeric()):
                        chrom_list.append(tok)
                if chrom_list:
                    chrom_file.write('{}\t{}\n'.format(name, ','.join(chrom_list)))
                    
    
if __name__ == '__main__':
    main()
