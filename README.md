# vgp-cactus
Notes about making vgp phase-1 alignments

Note that the `./gen-xxx-seqfile-vx.sh` scripts below require pandas and biopython.  You can install them with

```
virtualenv -p python3 venv-vgp
. venv-vgp/bin/activate
pip install -U pandas biopython
```

## 577-way

This an alignment of the `v1.1.16b` ROADIES tree plus `hg38` and `mm39` minus genomes that Cactus can't run in under 2TB:

```
common newt (GCA_964263255.1): 24,226,223,864 bp
palmate newt (GCA_964261635.1): 23,170,028,842 bp
warty newt (GCA_964204655.1): 22,324,616,549 bp
Iberian ribbed newt (GCA_026652325.1): 20,300,798,037 bp
axolotl (GCF_040938575.1): 29,117,848,270 bp
West African lungfish (GCA_040939525.1): 40,524,157,013 bp
```

and correcting the band-tailed pigeon assembly from `GCA_036971685.2` to `GCF_037038585.1`.

Create the seqfile. All commands run in `./577way`.

```
./gen-577way-seqfile-v1.sh
```

Create the cactus workflow script (used cactus v3.1.3).  Note that the ROADIES branch lengths tend to be much shorter than we expect, so add `--branchScale 2` to double them up. 
```
cactus-prepare vgp-577way-v1.seqfile --outDir vgp-577way-v1-prep --chromInfo vgp-577way-v1.chrom-info --branchScale 2 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 2000Gi' --script --outHal vgp-577way-v1.hal --preprocessBatchSize 1000 > vgp-577way-v1.sh
chmod +x vgp-577way-v1.sh
```

**Note**: I switched to `v3.1.4` and reduced `--maxMemory to 1700Gi` part way through.  These changes don't affect any logic, but make paf i/o a bit faster and help with scheduling. 

Export the MAFs with

```
cactus-hal2maf jobstore/js-maf-hg38 vgp-577way-v1.hal vgp-577way-v1-hg38.maf.gz --refGenome GCA_000001405.15 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-hg38.maf.gz.log &
cactus-hal2maf jobstore/js-maf-mm39 vgp-577way-v1.hal vgp-577way-v1-mm39.maf.gz --refGenome GCA_000001635.9 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mm39.maf.gz.log &
cactus-hal2maf jobstore/js-maf-mouse vgp-577way-v1.hal vgp-577way-v1-mouse.maf.gz --refGenome GCA_949316315.1 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mouse.maf.gz.log &
cactus-hal2maf jobstore/js-maf-hs1 vgp-577way-v1.hal vgp-577way-v1-hs1.maf.gz --refGenome GCA_009914755.4 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-hs1.maf.gz.log &
cactus-hal2maf jobstore/js-maf-chicken vgp-577way-v1.hal vgp-577way-v1-chicken.maf.gz --refGenome GCF_016700215.2 --outType norm single --index --coverage --coverageSexChroms NC_059535.1 NC_059536.1 --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-chicken.maf.gz.log &
cactus-hal2maf jobstore/js-maf-zebrafish vgp-577way-v1.hal vgp-577way-v1-zebrafish.maf.gz --refGenome GCA_944039275.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-zebrafish.maf.gz.log &
cactus-hal2maf jobstore/js-maf-spotted_gar vgp-577way-v1.hal vgp-577way-v1-spotted_gar.maf.gz --refGenome GCF_040954835.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-spotted_gar.maf.gz.log &
cactus-hal2maf jobstore/js-maf-catshark vgp-577way-v1.hal vgp-577way-v1-catshark.maf.gz --refGenome GCF_902713615.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-catshark.maf.gz.log &
cactus-hal2maf jobstore/js-maf-clawed_frog vgp-577way-v1.hal vgp-577way-v1-clawed_frog.maf.gz --refGenome GCA_038501925.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-clawed_frog.maf.gz.log &
cactus-hal2maf jobstore/js-maf-brown_anole vgp-577way-v1.hal vgp-577way-v1-brown_anole.maf.gz --refGenome GCF_037176765.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-brown_anole.maf.gz.log &
cactus-hal2maf jobstore/js-maf-horseshoe_bat vgp-577way-v1.hal vgp-577way-v1-horseshoe_bat.maf.gz --refGenome GCF_004115265.2 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-horseshoe_bat.maf.gz.log &
cactus-hal2maf jobstore/js-maf-dog vgp-577way-v1.hal vgp-577way-v1-dog.maf.gz --refGenome GCF_011100685.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-dog.maf.gz.log &
cactus-hal2maf jobstore/js-maf-zebra_finch vgp-577way-v1.hal vgp-577way-v1-zebra_finch.maf.gz --refGenome GCA_048771995.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-zebra_finch.maf.gz.log &
cactus-hal2maf jobstore/js-maf-emu vgp-577way-v1.hal vgp-577way-v1-emu.maf.gz --refGenome GCF_036370855.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-emu.maf.gz.log &
cactus-hal2maf jobstore/js-maf-european_eel vgp-577way-v1.hal vgp-577way-v1-european_eel.maf.gz --refGenome GCF_013347855.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-european_eel.maf.gz.log &
cactus-hal2maf jobstore/js-maf-eastern_happy vgp-577way-v1.hal vgp-577way-v1-eastern_happy.maf.gz --refGenome GCA_964374335.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-eastern_happy.maf.gz.log &
cactus-hal2maf jobstore/js-maf-three_spined_stickleback vgp-577way-v1.hal vgp-577way-v1-three_spined_stickleback.maf.gz --refGenome GCA_964276395.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-three_spined_stickleback.maf.gz.log &
cactus-hal2maf jobstore/js-maf-green_sea_turtle vgp-577way-v1.hal vgp-577way-v1-green_sea_turtle.maf.gz --refGenome GCF_015237465.2 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-green_sea_turtle.maf.gz.log &
cactus-hal2maf jobstore/js-maf-gray_short_tailed_opossum vgp-577way-v1.hal vgp-577way-v1-gray_short_tailed_opossum.maf.gz --refGenome GCF_027887165.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-gray_short_tailed_opossum.maf.gz.log &
cactus-hal2maf jobstore/js-maf-european_river_lamprey vgp-577way-v1.hal vgp-577way-v1-european_river_lamprey.maf.gz --refGenome GCA_964198595.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-european_river_lamprey.maf.gz.log &
cactus-hal2maf jobstore/js-maf-tiny_cayenne_caecilian vgp-577way-v1.hal vgp-577way-v1-tiny_cayenne_caecilian.maf.gz --refGenome GCA_901765095.2 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-tiny_cayenne_caecilian.maf.gz.log &
cactus-hal2maf jobstore/js-maf-coelacanth vgp-577way-v1.hal vgp-577way-v1-coelacanth.maf.gz --refGenome GCF_037176945.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-coelacanth.maf.gz.log &
cactus-hal2maf jobstore/js-maf-mexican_tetra vgp-577way-v1.hal vgp-577way-v1-mexican_tetra.maf.gz --refGenome GCF_023375975.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mexican_tetra.maf.gz.log &
cactus-hal2maf jobstore/js-maf-coastal_tailed_frog vgp-577way-v1.hal vgp-577way-v1-coastal_tailed_frog.maf.gz --refGenome GCA_040206685.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-coastal_tailed_frog.maf.gz.log &

wait
```

PhyloP tracks (with `cactus-phast`, cactus v3.2.0).  We run phast for any reference where we can get a gene annotation.  For `GCF_` accessions we pull annotations from NCBI by browsing https://www.ncbi.nlm.nih.gov/datasets/genome/<ACCESSION>/ and grabbing the GFF URL from the FTP path (`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../<asm>_genomic.gff.gz`) -- the `<asm>` segment varies per accession and needs to be checked by hand.  For the standard browser genomes (hg38, mm39, hs1) we use UCSC's gene tracks directly.

For each MAF we could generate three phyloP tracks via `--subtree`: the full tree root (`Anc0`), a broader lineage clade (e.g. Mammals, Birds), and a more local clade (e.g. Carnivora, Galloanserformes).  The picks below were chosen to be biologically meaningful without being too small; you can audit/replay them with the helper:

```
./pick-phast-subtrees.py vgp-577way.nwk GCF_037176765.1
```

```
cactus-phast jobstore/js-phast-hg38 vgp-577way-v1-hg38.single.maf.gz vgp-577way-v1.hal GCA_000001405.15 conservation/GCA_000001405.15/vgp-577way-v1-hg38.single.phyloP.wig.gz --geneAnnotation https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-hg38.phyloP.log &

slecactus-phast jobstore/js-phast-mm39 vgp-577way-v1-mm39.single.maf.gz vgp-577way-v1.hal GCA_000001635.9 conservation/GCA_000001635.9/vgp-577way-v1-mm39.single.phyloP.wig.gz --geneAnnotation https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mm39.phyloP.log &

cactus-phast jobstore/js-phast-hs1 vgp-577way-v1-hs1.single.maf.gz vgp-577way-v1.hal GCA_009914755.4 conservation/GCA_009914755.4/vgp-577way-v1-hs1.single.phyloP.wig.gz --geneAnnotation https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-hs1.phyloP.log &

cactus-phast jobstore/js-phast-chicken vgp-577way-v1-chicken.single.maf.gz vgp-577way-v1.hal GCF_016700215.2 conservation/GCF_016700215.2/vgp-577way-v1-chicken.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/700/215/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-chicken.phyloP.log &

cactus-phast jobstore/js-phast-emu vgp-577way-v1-emu.single.maf.gz vgp-577way-v1.hal GCF_036370855.1 conservation/GCF_036370855.1/vgp-577way-v1-emu.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/370/855/GCF_036370855.1_bDroNov1.hap1/GCF_036370855.1_bDroNov1.hap1_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-emu.phyloP.log &

cactus-phast jobstore/js-phast-spotted_gar vgp-577way-v1-spotted_gar.single.maf.gz vgp-577way-v1.hal GCF_040954835.1 conservation/GCF_040954835.1/vgp-577way-v1-spotted_gar.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/954/835/GCF_040954835.1_fLepOcu1.hap2/GCF_040954835.1_fLepOcu1.hap2_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-spotted_gar.phyloP.log &

cactus-phast jobstore/js-phast-european_eel vgp-577way-v1-european_eel.single.maf.gz vgp-577way-v1.hal GCF_013347855.1 conservation/GCF_013347855.1/vgp-577way-v1-european_eel.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/347/855/GCF_013347855.1_fAngAng1.pri/GCF_013347855.1_fAngAng1.pri_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-european_eel.phyloP.log &

cactus-phast jobstore/js-phast-zebrafish vgp-577way-v1-zebrafish.single.maf.gz vgp-577way-v1.hal GCA_944039275.1 conservation/GCA_944039275.1/vgp-577way-v1-zebrafish.single.phyloP.wig.gz --geneAnnotation https://genome-test.gi.ucsc.edu/~hiram/VGP/liftedGenes/GCA_944039275.1.ncbiRefSeq.gtf.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-zebrafish.phyloP.log &

cactus-phast jobstore/js-phast-catshark vgp-577way-v1-catshark.single.maf.gz vgp-577way-v1.hal GCF_902713615.1 conservation/GCF_902713615.1/vgp-577way-v1-catshark.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/713/615/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-catshark.phyloP.log &

cactus-phast jobstore/js-phast-brown_anole vgp-577way-v1-brown_anole.single.maf.gz vgp-577way-v1.hal GCF_037176765.1 conservation/GCF_037176765.1/vgp-577way-v1-brown_anole.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/765/GCF_037176765.1_rAnoSag1.mat/GCF_037176765.1_rAnoSag1.mat_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-brown_anole.phyloP.log &

cactus-phast jobstore/js-phast-horseshoe_bat vgp-577way-v1-horseshoe_bat.single.maf.gz vgp-577way-v1.hal GCF_004115265.2 conservation/GCF_004115265.2/vgp-577way-v1-horseshoe_bat.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/115/265/GCF_004115265.2_mRhiFer1_v1.p/GCF_004115265.2_mRhiFer1_v1.p_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-horseshoe_bat.phyloP.log &

cactus-phast jobstore/js-phast-dog vgp-577way-v1-dog.single.maf.gz vgp-577way-v1.hal GCF_011100685.1 conservation/GCF_011100685.1/vgp-577way-v1-dog.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/685/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-dog.phyloP.log &

cactus-phast jobstore/js-phast-gray_short_tailed_opossum vgp-577way-v1-gray_short_tailed_opossum.single.maf.gz vgp-577way-v1.hal GCF_027887165.1 conservation/GCF_027887165.1/vgp-577way-v1-gray_short_tailed_opossum.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/887/165/GCF_027887165.1_mMonDom1.pri/GCF_027887165.1_mMonDom1.pri_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-gray_short_tailed_opossum.phyloP.log &

cactus-phast jobstore/js-phast-green_sea_turtle vgp-577way-v1-green_sea_turtle.single.maf.gz vgp-577way-v1.hal GCF_015237465.2 conservation/GCF_015237465.2/vgp-577way-v1-green_sea_turtle.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/237/465/GCF_015237465.2_rCheMyd1.pri.v2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-green_sea_turtle.phyloP.log &

cactus-phast jobstore/js-phast-mexican_tetra vgp-577way-v1-mexican_tetra.single.maf.gz vgp-577way-v1.hal GCF_023375975.1 conservation/GCF_023375975.1/vgp-577way-v1-mexican_tetra.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/375/975/GCF_023375975.1_AstMex3_surface/GCF_023375975.1_AstMex3_surface_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mexican_tetra.phyloP.log &

cactus-phast jobstore/js-phast-coelacanth vgp-577way-v1-coelacanth.single.maf.gz vgp-577way-v1.hal GCF_037176945.1 conservation/GCF_037176945.1/vgp-577way-v1-coelacanth.single.phyloP.wig.gz --geneAnnotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/176/945/GCF_037176945.1_fLatCha1.pri/GCF_037176945.1_fLatCha1.pri_genomic.gff.gz --substMod REV --modFreqs --precision HIGH --bigwig --batchSystem slurm --chunkCores 32 --phyloFitCores 64 --slurmPartition medium --doubleMem true --slurmTime 11:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-coelacanth.phyloP.log &

```

Chains
```
sbatch --partition=long --time=10-00:00:00 --mem=1800G --cpus-per-task=32 --job-name=chains-hg38 --wrap="cactus-hal2chains /data/tmp/js-chain-hg38 vgp-577way-v1.hal vgp-577way-v1-hg38-chains --targetGenomes GCA_000001405.15  --bigChain --retryCount 5 --maxCores 32 --maxMemory 1800Gi --symlinkImports=False" --error vgp-577way-v1-prep/logs/vgp-577way-v1-hg38-chains.log

sbatch --partition=long --time=10-00:00:00 --mem=1800G --cpus-per-task=32 --job-name=chains-catshark --wrap="cactus-hal2chains /data/tmp/js-chain-catshark vgp-577way-v1.hal vgp-577way-v1-catshark-chains --targetGenomes GCF_902713615.1  --bigChain --retryCount 5 --maxCores 32 --maxMemory 1800Gi --symlinkImports=False" --error vgp-577way-v1-prep/logs/vgp-577way-v1-catshark-chains.log

sbatch --partition=long --time=10-00:00:00 --mem=1800G --cpus-per-task=32 --job-name=chains-chicken --wrap="cactus-hal2chains /data/tmp/js-chain-chicken vgp-577way-v1.hal vgp-577way-v1-chicken-chains --targetGenomes GCF_016700215.2 --bigChain --retryCount 5 --maxCores 32 --maxMemory 1800Gi --symlinkImports=False" --error vgp-577way-v1-prep/logs/vgp-577way-v1-chicken-chains.log

sbatch --partition=long --time=10-00:00:00 --mem=1800G --cpus-per-task=32 --job-name=chains-hg38-13way-demo --wrap="cactus-hal2chains /data/tmp/js-chain-hg38-13way-demo vgp-577way-v1.hal vgp-577way-v1-hg38-13way-demo-chains --targetGenomes GCA_000001405.15 --queryGenomes GCA_949316315.1 GCA_009914755.4 GCF_016700215.2 GCA_944039275.1 GCF_040954835.1 GCF_902713615.1 GCA_038501925.1 GCF_037176765.1 GCF_004115265.2 GCF_011100685.1 GCA_048771995.1 GCF_036370855.1 GCF_013347855.1 --bigChain --retryCount 5 --maxCores 32 --maxMemory 1800Gi --symlinkImports=False" --error vgp-577way-v1-hg38-13way-demo-chains.log

```


### Coverage Plots

The MAF export commands above produce a `*.maf.gz.cov.tsv` file alongside each MAF (via `--coverage`). For every reference, four bar-chart PNGs are generated from that table using `scripts/visualize_coverage.R`. Each plot shows per-species alignment coverage (in Mbp), with bars colored by Lineage.

The script depends on the `annotations.tsv` table (lineage / family / english-name lookup) which lives in the per-alignment coverage directory.

Reproduce the plots for a single reference (using hs1 as example):

```
cd 577way/coverage
Rscript ../../scripts/visualize_coverage.R vgp-577way-v1-hs1.maf.gz.cov.tsv annotations.tsv Total
```

This produces four PNGs in the working directory:

- `vgp-577way-v1-hs1.maf.gz.cov.tsv_coverage_total.png` — linear-scale total alignment coverage; light bars = total aligned bases of hs1 covered by the query, dark overlay = identical matches.
- `vgp-577way-v1-hs1.maf.gz.cov.tsv_coverage_total_log.png` — same as above with a log-scaled y-axis (long tail of distantly related species becomes visible).
- `vgp-577way-v1-hs1.maf.gz.cov.tsv_coverage_total_1to1.png` — dark overlay now denotes 1:1 (single-copy, no duplications) coverage instead of identical matches.
- `vgp-577way-v1-hs1.maf.gz.cov.tsv_coverage_total_1to1_log.png` — same as above with a log-scaled y-axis.

To regenerate plots for every reference in the 577way coverage directory:

```
cd 577way/coverage
for f in vgp-577way-v1-*.maf.gz.cov.tsv; do
  [[ "$f" == *.single.* ]] && continue
  Rscript ../../scripts/visualize_coverage.R "$f" annotations.tsv Total
done
```

Optional 4th argument: a regex matched against the `Orders` column to restrict bars to a specific clade (when active, bars are colored by Family rather than Lineage). For example, to show only bats on the horseshoe bat reference:

```
Rscript ../../scripts/visualize_coverage.R vgp-577way-v1-horseshoe_bat.maf.gz.cov.tsv annotations.tsv Total Chiroptera
```


## Older alignments

### Amniotes

Create the amniotes (mammal/birds/reptile) seqfile. Commands run in `./amniotes`. Note that the branches in the ROADIES tree seem short, so this script scales them by 2 to help cactus select appropriate lastz parameters. 

```
./gen-amniotes-seqfile-v1.sh
```

Create the amniotes workflow script (used cactus v3.1.2).  Note that the ROADIES branch lengths tend to be much shorter than we expect, so add `--branchScale 2` to double them up. 
```
cactus-prepare amniotes-v1.seqfile --outDir amniotes-v1-prep --chromInfo amniotes-v1.chrom-info --branchScale 2 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1600Gi' --script --outHal vgp-amniotes-v1.hal --preprocessBatchSize 1000 > amniotes-v1.sh
chmod +x amniotes-v1.sh
```

Don't align or include the outgroups (which are under Anc002 and Anc000)
```
sed -i amniotes-v1.sh -e 's/Anc000.hal.append/Anc001.hal.append/g'
sed -i amniotes-v1.sh -e '/^cactus-halAppendSubtrees/s|amniotes-v1-prep/Anc000\.hal||g' -e '/^cactus-halAppendSubtrees/s|amniotes-v1-prep/Anc002\.hal||g'
sed -i amniotes-v1.sh -e '/Anc000/d' -e '/Anc002/d'
sed -i amniotes-v1.sh -e '/^pids=()$/{N;s/\npids+=(\$!)$//}'
```

Run the workflow
```
./amniotes-v1.sh
```

Make a backup
```
cp vgp-amniotes-v1.hal vgp-amniotes-v1.hal.bak
```

Rename the Ancestors in perparation of future merge.
```
printf "Anc001\tAmniotesAnc000\n" > amniotes-rename.tsv
for i in `seq 3 360`; do printf "Anc%03d\tAmniotesAnc%03d\n" $i $((i-2)); done >> amniotes-rename.tsv
halRenameGenomes vgp-amniotes-v1.hal amniotes-rename.tsv
```

Export the MAFs with

```
cactus-hal2maf jobstore/js-maf-hg38 vgp-amniotes-v1.hal vgp-amniotes-v1-hg38.maf.gz --refGenome hg38 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-hg38.maf.gz.log & 
cactus-hal2maf jobstore/js-maf-mm39 vgp-amniotes-v1.hal vgp-amniotes-v1-mm39.maf.gz --refGenome mm39 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-mm39.maf.gz.log &
cactus-hal2maf jobstore/js-maf-hs1 vgp-amniotes-v1.hal vgp-amniotes-v1-hs1.maf.gz --refGenome GCA_009914755.4 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-hs1.maf.gz.log &
cactus-hal2maf jobstore/js-maf-chicken vgp-amniotes-v1.hal vgp-amniotes-v1-chicken.maf.gz --refGenome GCF_016700215.2 --outType norm single --index --coverage --coverageSexChroms NC_059535.1 NC_059536.1 --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-chicken.maf.gz.log &

wait
```

And the chains with (note: these have to be submitted via slurm!)
```
sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=160 --job-name=chains-hg38 --wrap="cactus-hal2chains /data/tmp/js-chain-hg381 vgp-amniotes-fish-v1.hal vgp-amniotes-v1-hg38-chains --targetGenomes hg38 --bigChain --retryCount 5 --maxCores 160 --maxMemory 1Ti --symlinkImports=False" --error amniotes-fish-v1-prep/logs/vgp-amniotes-v1-hg38-chains.log
sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=160 --job-name=chains-mm39 --wrap="cactus-hal2chains /data/tmp/js-chain-mm391 vgp-amniotes-fish-v1.hal vgp-amniotes-v1-mm39-chains --targetGenomes mm39 --bigChain --retryCount 5 --maxCores 160 --maxMemory 1Ti --symlinkImports=False" --error amniotes-fish-v1-prep/logs/vgp-amniotes-v1-mm39-chains.log
sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=160 --job-name=chains-hs1 --wrap="cactus-hal2chains /data/tmp/js-chain-hs11 vgp-amniotes-fish-v1.hal vgp-amniotes-v1-hs1-chains --targetGenomes GCA_009914755.4 --bigChain --retryCount 5 --maxCores 160 --maxMemory 1Ti --symlinkImports=False" --error amniotes-fish-v1-prep/logs/vgp-amniotes-v1-hs1-chains.log
sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=160 --job-name=chains-chicken --wrap="cactus-hal2chains /data/tmp/js-chain-chicken1 vgp-amniotes-fish-v1.hal vgp-amniotes-v1-chicken-chains --targetGenomes GCF_016700215.2 --bigChain --retryCount 5 --maxCores 160 --maxMemory 1Ti --symlinkImports=False" --error amniotes-fish-v1-prep/logs/vgp-amniotes-v1-chicken-chains.log

```


### Rayfin-Fish

Create the rayfin-fish seqfile. Commands run in `./rayfin-fish`:

```
./gen-rayfin-fish-seqfile-v1.sh
```

Create the rayfin-fish workflow script (used cactus v3.1.2). We scale branches by 4 (twice as much as amniotes) because in addition to having shorter-than-expected branches, the fish are so difficult to align.  

```
cactus-prepare rayfin-fish-v1.seqfile --outDir rayfin-fish-v1-prep --chromInfo rayfin-fish-v1.chrom-info --branchScale 4 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1600Gi' --script --outHal rayfin-fish-v1.hal  --preprocessBatchSize 1000 > rayfin-fish-v1.sh
```

Don't align or include the outgroups (which are under Anc000 and Anc001)
```
sed -i rayfin-fish-v1.sh -e 's/Anc000.hal.append/Anc002.hal.append/g'
sed -i rayfin-fish-v1.sh -e '/^cactus-halAppendSubtrees/s|rayfin-fish-v1-prep/Anc000\.hal||g' -e '/^cactus-halAppendSubtrees/s|rayfin-fish-v1-prep/Anc001\.hal||g'
sed -i rayfin-fish-v1.sh -e '/Anc000/d' -e '/Anc001/d'
sed -i rayfin-fish-v1.sh -e '/^pids=()$/{N;s/\npids+=(\$!)$//}'
```

Rename the Ancestors in perparation of future merge.
```
for i in `seq 2 157`; do printf "Anc%03d\tRayfinAnc%03d\n" $i $((i-2)); done >> rayfin-rename.tsv
halRenameGenomes vgp-rayfin-v1.hal rayfin-rename.tsv
```

Export the MAFs with

```
cactus-hal2maf jobstore/js-maf-zebrafish vgp-rayfin-v1.hal vgp-rayfin-v1-zebrafish.maf.gz --refGenome GCA_944039275.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --logFile rayfin-fish-v1-prep/logs/vgp-rayfin-v1-zebrafish.maf.gz.log

cactus-hal2maf jobstore/js-maf-fugu vgp-rayfin-v1.hal vgp-rayfin-v1-fugu.maf.gz --refGenome GCF_901000725.2 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --logFile rayfin-fish-v1-prep/logs/vgp-rayfin-v1-fugu.maf.gz.log

cactus-hal2maf jobstore/js-maf-cichlid vgp-rayfin-v1.hal vgp-rayfin-v1-cichlid.maf.gz --refGenome GCA_964374335.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --logFile rayfin-fish-v1-prep/logs/vgp-rayfin-v1-cichlid.maf.gz.log

cactus-hal2maf jobstore/js-maf-spotted_gar vgp-rayfin-fish-v1.hal vgp-rayfin-v1-spotted_gar.maf.gz --refGenome GCF_040954835.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --logFile rayfin-fish-v1-prep/logs/vgp-rayfin-v1-spotted_gar.maf.gz.log
```

And do the chains.  Note: using `sbatch` and `--symlinkImports=False` to make sure that only one copy of the hal is ever made.
```
sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=96 --job-name=chains-zebrafish --wrap="cactus-hal2chains /data/tmp/js-chain-zebrafish1 vgp-rayfin-fish-v1.hal vgp-rayfin-v1-zebrafish-chains --targetGenomes GCA_944039275.1 --bigChain --retryCount 5 --maxCores 96 --maxMemory 1Ti --symlinkImports=False" --error rayfin-fish-v1-prep/logs/vgp-rayfin-v1-zebrafish-chains.log

sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=96 --job-name=chains-fugu --wrap="cactus-hal2chains /data/tmp/js-chain-fugu vgp-rayfin-fish-v1.hal vgp-rayfin-v1-fugu-chains --targetGenomes GCF_901000725.2 --bigChain --retryCount 5 --maxCores 96 --maxMemory 1Ti --symlinkImports=False" --error rayfin-fish-v1-prep/logs/vgp-rayfin-v1-fugu-chains.log

sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=96 --job-name=chains-cichlid --wrap="cactus-hal2chains /data/tmp/js-chain-cichlid vgp-rayfin-fish-v1.hal vgp-rayfin-v1-cichlid-chains --targetGenomes GCA_964374335.1 --bigChain --retryCount 5 --maxCores 96 --maxMemory 1Ti --symlinkImports=False" --error rayfin-fish-v1-prep/logs/vgp-rayfin-v1-cichlid-chains.log

sbatch --partition=long --time=5-00:00:00 --mem=1T --cpus-per-task=96 --job-name=chains-spotted_gar --wrap="cactus-hal2chains /data/tmp/js-chain-spotted_gar vgp-rayfin-fish-v1.hal vgp-rayfin-v1-spotted_gar-chains --targetGenomes GCF_040954835.1 --bigChain --retryCount 5 --maxCores 96 --maxMemory 1Ti --symlinkImports=False" --error rayfin-fish-v1-prep/logs/vgp-rayfin-v1-spotted_gar-chains.log

```



