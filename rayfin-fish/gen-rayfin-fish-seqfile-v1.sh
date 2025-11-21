#!/usr/bin/bash
set -ex
BASE_DIR=..
BIN=${BASE_DIR}/scripts
TAB=${BASE_DIR}/tables

# download the input tables and tree from github
${BIN}/download-tables.py --out-dir ${TAB}

# extract a subtree, and add human-readable suffixes
# we don't select amphibian outgroups (--ignore) due to concerns about their genome sizes
${BIN}/tree-extract.py --tree ${TAB}/roadies_v1.1.16b.nwk --table ${TAB}/VGPPhase1-freeze-1.0.tsv --category Superorder --value "Actinopterygii (Ray-finned Fishes)" --outgroups 2 --suffix-category ScientificName  > rayfin-fish-v1.nwk

# convert the subtree to a seqfile, keeping only the suffixes as names (ie dropping accessions)
${BIN}/tree2seqfile.py --tree rayfin-fish-v1.nwk --urls ${TAB}/URL.download.table.tsv --table ${TAB}/VGPPhase1-freeze-1.0.tsv --chrom-info rayfin-fish-v1.chrom-info --suffix-only > rayfin-fish-v1.seqfile

