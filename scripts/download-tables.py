#!/usr/bin/env python3
"""
Download some tables necessary to make the Cactus input
"""

import os, sys
import subprocess
import argparse

fasta_table='https://raw.githubusercontent.com/VGP/vgp-phase1/c5252eebc52f8cddbfdf90e4e12b43a60aa638a5/URL.download.table.tsv'
main_table='https://raw.githubusercontent.com/VGP/vgp-phase1/c5252eebc52f8cddbfdf90e4e12b43a60aa638a5/VGPPhase1-freeze-1.0.tsv'
newick_tree='https://raw.githubusercontent.com/VGP/vgp-trees/45fa9abcf6c9aa99cfa1b110b096457888e2a44a/phase-1/roadies_v1.1.16b.nwk'
annotations='https://raw.githubusercontent.com/VGP/vgp-trees/45fa9abcf6c9aa99cfa1b110b096457888e2a44a/phase-1/annotations.tsv'

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Download some tables necessary to make the Cactus input')
    parser.add_argument('--out-dir', required=True,
                        help='tree to read in newick format')
    
    args = parser.parse_args(command_line)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    for input_url in [fasta_table, newick_tree, annotations, main_table]:
        subprocess.check_call(['wget', '-q', input_url, '-O', os.path.join(args.out_dir, os.path.basename(input_url))])
    
if __name__ == '__main__':
    main()
