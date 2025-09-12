#!/usr/bin/bash

BASE_DIR=..
BIN=${BASE_DIR}/scripts
TAB=${BASE_DIR}/tables

# download the input tables and tree from github
${BIN}/download-tables.py --out-dir ${TAB}

# extract a subtree, and add human-readable suffixes
${BIN}/tree-extract.py --tree ${TAB}/roadies_v1.1.4.nwk --table ${TAB}/annotations.tsv --category Lineage --value Mammals --outgroups 1 --suffix-category ScientificName > mammals-v1.nwk

# convert the subtree to a seqfile
${BIN}/tree2seqfile.py --tree mammals-v1.nwk --urls ${TAB}/URL.download.table.tsv --table ${TAB}/VGPPhase1-freeze-1.0.tsv --chrom-info mammals-v1.chrom-info > mammals-v1.seqfile

