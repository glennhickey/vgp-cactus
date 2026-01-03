#!/usr/bin/bash
set -ex
BASE_DIR=..
BIN=${BASE_DIR}/scripts
TAB=${BASE_DIR}/tables

# download the input tables and tree from github
${BIN}/download-tables.py --out-dir ${TAB}

# extract a subtree, leave names as accessions
${BIN}/tree-extract.py --tree ${TAB}/roadies_v1.1.16b.nwk --table ${TAB}/VGPPhase1-freeze-1.0.tsv --category Lineage --value Mammals --value Reptiles --value Birds --outgroups 2 > amniotes-v1.nwk

# convert the subtree to a seqfile
${BIN}/tree2seqfile.py --tree amniotes-v1.nwk --urls ${TAB}/URL.download.table.tsv --table ${TAB}/VGPPhase1-freeze-1.0.tsv --chrom-info amniotes-v1.chrom-info --suffix-only > amniotes-v1.seqfile

# rename Homo_sapiens to hs1 (in the tree and 1st column, but not second column)
sed -i amniotes-v1.seqfile -e 's/GCA_009914755.4:/hs1:/g' -e 's/^GCA_009914755.4/hs1/g'

# add mm39
${BIN}/add-genome-to-seqfile.py --seqfile amniotes-v1.seqfile --target GCA_949316315.1 --name mm39 --url https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz --chroms X,Y --chrom-info amniotes-v1.chrom-info > amniotes-v1.seqfile.new
mv amniotes-v1.chrom-info.new amniotes-v1.chrom-info
mv amniotes-v1.seqfile.new amniotes-v1.seqfile

# add hg38
${BIN}/add-genome-to-seqfile.py --seqfile amniotes-v1.seqfile --target hs1 --name hg38 --url https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz --chroms X,Y --chrom-info amniotes-v1.chrom-info > amniotes-v1.seqfile.new
mv amniotes-v1.chrom-info.new amniotes-v1.chrom-info
mv amniotes-v1.seqfile.new amniotes-v1.seqfile

