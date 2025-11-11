# vgp-cactus
Notes about making vgp phase-1 alignments

Note that the scripts below require pandas and biopython.  You can install them with

```
virtualenv -p python3 venv-vgp
. venv-vgp/bin/activate
pip install -U pandas biopython
```

## Mammals

Create the mammals seqfile (used mammals-v1.0 tag of this repo). Commands run in `./mammals`:

```
./gen-mammals-seqfile-v1.sh
```

Create the mammals workflow script (used cactus v3.0.0)

```
cactus-prepare mammals-v1.seqfile --outDir mammals-v1-prep --chromInfo mammals-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal vgp-mammals-v1.hal > mammals-v1.sh
```

Don't put the outgroups in final tree
```
sed -i mammals-v1.sh -e 's#mammals-v1-prep/Anc000.hal mammals-v1-prep/Anc001.hal ##g'
```

Note: needed to make the following manual changes:
* update to Toil `740508407ea34388f1a6f2d412964b67483ae9e8`
* add `--skipProgress` and `--maxMemory 1Ti` to cactus commands
* Toggle off FastGA options for: `Anc002`, `Anc039`, `Anc065`, `Anc102`, `Anc131`, `Anc145`
* Correction for accidentially using the wrong hg38 assembly the first time (re-align entire hg38-to-root path)

Rename the Ancestors in perparation of future merge (note we start at 2 because of outgroups).
There's also a bug in the names where table formatting issues led to trailing underscores
Rename Homo_sapiens to hs1 to be a little more clear about it (since hg38 is also in the alignment)
```
for i in `seq 2 168`; do printf "Anc%03d\tMammalsAnc%03d\n" $i $((i-2)); done > mammals-rename.tsv
for i in Molossus_alvarezi Sciurus_carolinensis Pongo_abelii Cnephaeus_nilssonii; do printf "${i}_\t${i}\n"; done >> mammals-rename.tsv
printf "Homo_sapiens\ths1\n" >> mammals-rename.tsv
halRenameGenomes vgp-mammals-v1.hal mammals-rename.tsv
```

Export the TAFs with

```
cactus-hal2maf jobstore/js-maf-hg38 vgp-mammals-v1.hal vgp-mammals-v1-hg38.taf.gz --refGenome hg38 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile mammals-v1-prep/logs/vgp-mammals-v1-hg38.taf.gz.log & 
cactus-hal2maf jobstore/js-maf-mm39 vgp-mammals-v1.hal vgp-mammals-v1-mm39.taf.gz --refGenome mm39 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile mammals-v1-prep/logs/vgp-mammals-v1-mm39.taf.gz.log &
cactus-hal2maf jobstore/js-maf-hs1 vgp-mammals-v1.hal vgp-mammals-v1-hs1.taf.gz --refGenome hs1 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile mammals-v1-prep/logs/vgp-mammals-v1-hs1.taf.gz.log &
wait
```

And the chains with (note: these have to be submitted via slurm!)
```
cactus-hal2chains jobstore/js-chain-hs1 vgp-mammals-v1.hal vgp-mammals-v1-chains-hs1 --targetGenomes hs1 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-hs1.taf.gz.log
cactus-hal2chains jobstore/js-chain-hg38 vgp-mammals-v1.hal vgp-mammals-v1-chains-hg38 --targetGenomes hg38 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-hg38.taf.gz.log
cactus-hal2chains jobstore/js-chain-mm39 vgp-mammals-v1.hal vgp-mammals-v1-chains-mm39 --targetGenomes mm39 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-mm39.taf.gz.log 
```

Make note the accessions for Hiram
```
cat mammals-v1.seqfile  | grep -v Anc | sed -e 's#mammals-v1-prep/##g' -e 's/.fa.gz//g' | grep -v mm39 | grep -v hg38 > vgp-mammals-v1.accessions
printf "mm39\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz\n" >> vgp-mammals-v1.accessions
printf "hg38\thttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz\n" >> vgp-mammals-v1.accessions
for i in Molossus_alvarezi Sciurus_carolinensis Pongo_abelii Cnephaeus_nilssonii; do sed -i vgp-mammals-v1.accessions -e "s/${i}_/${i}/g"; done
```

## Birds

Create the birds seqfile (used birds-v1.0 tag of this repo). Commands run in `./birds`:

```
./gen-birds-seqfile-v1.sh
```

Create the birds workflow script (used cactus v3.0.1)

```
cactus-prepare birds-v1.seqfile --outDir birds-v1-prep --chromInfo birds-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal birds-v1.hal > birds-v1.sh
```

We don't want the outgroups, so pull them out of the final alignment
```
sed -i birds-v1.sh -e 's#birds-v1-prep/Anc000.hal birds-v1-prep/Anc001.hal ##g'
```



