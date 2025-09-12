#!/usr/bin/bash
set -ex
BASE_DIR=..
BIN=${BASE_DIR}/scripts
TAB=${BASE_DIR}/tables

# download the input tables and tree from github
${BIN}/download-tables.py --out-dir ${TAB}

# extract a subtree, and add human-readable suffixes
${BIN}/tree-extract.py --tree ${TAB}/roadies_v1.1.4.nwk --table ${TAB}/annotations.tsv --category Lineage --value Mammals --outgroups 1 --suffix-category ScientificName > mammals-v1.nwk

# convert the subtree to a seqfile
${BIN}/tree2seqfile.py --tree mammals-v1.nwk --urls ${TAB}/URL.download.table.tsv --table ${TAB}/VGPPhase1-freeze-1.0.tsv --chrom-info mammals-v1.chrom-info > mammals-v1.seqfile

# add mm19
${BIN}/add-genome-to-seqfile.py --seqfile mammals-v1.seqfile --target GCA_949316315.1-Mus_musculus --name mm39 --url https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz --chroms X,Y --chrom-info mammals-v1.chrom-info > mammals-v1.seqfile.new
mv mammals-v1.chrom-info.new mammals-v1.chrom-info
mv mammals-v1.seqfile.new mammals-v1.seqfile

# add hg38
${BIN}/add-genome-to-seqfile.py --seqfile mammals-v1.seqfile --target GCA_009914755.4-Homo_sapiens --name hg38 --url https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz --chroms X,Y --chrom-info mammals-v1.chrom-info > mammals-v1.seqfile.new
mv mammals-v1.chrom-info.new mammals-v1.chrom-info
mv mammals-v1.seqfile.new mammals-v1.seqfile

