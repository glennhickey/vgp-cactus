#!/usr/bin/bash
set -ex
BASE_DIR=..
BIN=${BASE_DIR}/scripts
TAB=${BASE_DIR}/tables

# download the input tables and tree from github
${BIN}/download-tables.py --out-dir ${TAB}

# filter out samples with big genomes
${BIN}/tree-extract.py  --tree ${TAB}/roadies_v1.1.16b.nwk --table ${TAB}/VGPPhase1-freeze-1.0.tsv --max-genome-size 12000000000 --write-accession-table vgp-577way-v1-scinames.tsv > vgp-575way-v1.raw.nwk

# convert the tree to a seqfile
${BIN}/tree2seqfile.py --tree vgp-575way-v1.raw.nwk --urls ${TAB}/URL.download.table.tsv --table ${TAB}/VGPPhase1-freeze-1.0.tsv --chrom-info vgp-577way-v1.chrom-info --suffix-only > vgp-577way-v1.raw.seqfile

# add mm39
${BIN}/add-genome-to-seqfile.py --seqfile vgp-577way-v1.raw.seqfile --target GCA_949316315.1 --name GCA_000001635.9 --url https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz --chroms X,Y --chrom-info vgp-577way-v1.chrom-info > vgp-577way-v1.raw.seqfile.new
mv vgp-577way-v1.chrom-info.new vgp-577way-v1.chrom-info
mv vgp-577way-v1.raw.seqfile.new vgp-577way-v1.raw.seqfile
printf "GCA_000001635.9\tmm39\n" >> vgp-577way-v1-scinames.tsv

# add hg38 (same as used in HPRC)
${BIN}/add-genome-to-seqfile.py --seqfile vgp-577way-v1.raw.seqfile --target GCA_009914755.4 --name GCA_000001405.15 --url https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz  --chroms X,Y --chrom-info vgp-577way-v1.chrom-info > vgp-577way-v1.raw.seqfile.new
mv vgp-577way-v1.chrom-info.new vgp-577way-v1.chrom-info
mv vgp-577way-v1.raw.seqfile.new vgp-577way-v1.raw.seqfile
printf "GCA_000001405.15\thg38\n" >> vgp-577way-v1-scinames.tsv

# label the internal nodes using classification info from the table
head -1 vgp-577way-v1.raw.seqfile > vgp-577way-v1.raw.nwk
${BIN}/label-tree.py  --tree vgp-577way-v1.raw.nwk --table ${TAB}/VGPPhase1-freeze-1.0.tsv  --combine 'Vertebrates:Lineage:Mammals,Birds,Reptiles,Amphibians,Fishes' --fallback-label Anc > vgp-577way.nwk
{ head -1 vgp-577way.nwk; tail -n +2 vgp-577way-v1.raw.seqfile; } >  vgp-577way-v1.seqfile




