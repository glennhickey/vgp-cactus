#!/usr/bin/env python3
"""
Add a new genome to a seqfile, updating the tree an url list

Note that the genome to be added will be a direct sister to a genome already in the tree, and should be very similar (ie same species)
"""

import os, sys
import subprocess
import argparse
from Bio import Phylo
import io
import re

mash_docker = 'quay.io/comparative-genomics-toolkit/cactus:v2.9.9'
def main(command_line=None):                     
    parser = argparse.ArgumentParser('Add a new genome to a seqfile, updating the tree an url list')
    parser.add_argument('--seqfile', required=True,
                        help='seqfile to update')
    parser.add_argument('--target', required=True,
                        help='name of target species already in the seqfile to add a sister to')
    parser.add_argument('--name', required=True,
                        help='name of new species to add')
    parser.add_argument('--url', required=True,
                        help='url of new species to add')
    parser.add_argument('--chroms',
                        help='comma-separated list of chromosomes')
    parser.add_argument('--chrom-info',
                        help='chrom-info file to update with chromosomes (output will have .new suffix)')
                               
    args = parser.parse_args(command_line)
    if args.chroms is None != args.chrom_info is None:
        sys.stderr.write('--chroms / --chrom-info must be used together\n')
        return 1

    # parse the seqfile
    tree_string = None
    tree = None
    fasta_map = {}
    order = []
    with open(args.seqfile, 'r') as seqfile:
        for line in seqfile:
            if '(' in line:
                assert not tree
                tree_string = line.strip()
                tree_handle = io.StringIO(line)
                tree = Phylo.read(tree_handle, 'newick')
            else:
                toks = line.split()
                if len(toks) not in [0,2]:
                    print(toks)
                assert len(toks) in [0,2]
                if len(toks) == 2:
                    fasta_map[toks[0]] = toks[1].strip()
                    order.append(toks[0])

    if args.target not in fasta_map:
        sys.stderr.write(f'target genome {args.target} not in seqfile\n')
        return 1
    if len([t for t in tree.find_elements(name=args.target)]) == 0:
        sys.stderr.write(f'target genome {args.target} not in tree\n')
    
    # download the genomes
    fa_paths = []
    for download_url in [args.url, fasta_map[args.target]]:
        if not os.path.isfile(os.path.basename(download_url)):
            sys.stderr.write(f'downloading {download_url}\n')
            subprocess.check_call(['wget', '-q', download_url])
        else:
            sys.stderr.write(f'skipping download of {download_url}\n')
        fa_paths.append(os.path.basename(download_url))
        
    # compute the mash distance
    sys.stderr.write(f'computing mash distance between {args.target} and {args.name}\n')
    dist = subprocess.check_output(['docker', 'run', '-it', '--rm', '-v', os.getcwd() + ':/data',
                                    '-u', f'{os.getuid()}:{os.getgid()}',
                                    mash_docker, 'mash', 'dist',
                                    '-p', '8', os.path.join('/data', fa_paths[0]),
                                    os.path.join('/data/', fa_paths[1])]).decode().strip()
    dist = float(dist.split()[-3])
    sys.stderr.write(f'distance is {dist}\n')

    # update the tree
    # can't figure out how to update the tree using the biopython api, just do the string instead...
    
    # find the node in the tree
    m = re.search(f'{args.target}:[0-9,.]*', tree_string)
    tok = m.group()
    # pull out its distance
    tree_dist = float(tok[tok.index(':') + 1:])
    assert dist < tree_dist
    # make a cherry (pair of nodes) to replace it with.
    # the distance to the ancestor is the same for each new node
    new_cherry_string = f'({args.target}:{dist/2},{args.name}:{dist/2}):{tree_dist - dist/2}'
    new_tree_string = tree_string.replace(tok, new_cherry_string)

    # print the new seqfile
    print(new_tree_string)
    print()
    for genome in order:
        print(f'{genome}\t{fasta_map[genome]}')
    print(f'{args.name}\t{args.url}')

    # update the chrom-info file
    with open(args.chrom_info, 'r') as in_chrom_file, open(args.chrom_info + '.new', 'w') as out_chrom_file:
        for line in in_chrom_file:
            out_chrom_file.write(line)
        out_chrom_file.write(f'{args.name}\t{args.chroms}\n')
    

                    
if __name__ == '__main__':
    main()
