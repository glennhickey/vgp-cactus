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
cactus-hal2maf jobstore/js-maf-mm39 vgp-577way-v1.hal vgp-577way-v1-mm39.maf.gz --refGenome GCA_949316315.1 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-mm39.maf.gz.log &
cactus-hal2maf jobstore/js-maf-hs1 vgp-577way-v1.hal vgp-577way-v1-hs1.maf.gz --refGenome GCA_009914755.4 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-hs1.maf.gz.log &
cactus-hal2maf jobstore/js-maf-chicken vgp-577way-v1.hal vgp-577way-v1-chicken.maf.gz --refGenome GCF_016700215.2 --outType norm single --index --coverage --coverageSexChroms NC_059535.1 NC_059536.1 --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-chicken.maf.gz.log &
cactus-hal2maf jobstore/js-maf-zebrafish vgp-577way-v1.hal vgp-577way-v1-zebrafish.maf.gz --refGenome GCA_944039275.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-zebrafish.maf.gz.log &
cactus-hal2maf jobstore/js-maf-spotted_gar vgp-577way-v1.hal vgp-577way-v1-spotted_gar.maf.gz --refGenome GCF_040954835.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-spotted_gar.maf.gz.log &
cactus-hal2maf jobstore/js-maf-catshark vgp-577way-v1.hal vgp-577way-v1-catshark.maf.gz --refGenome GCF_902713615.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-catshark.maf.gz.log &
cactus-hal2maf jobstore/js-maf-clawed_frog vgp-577way-v1.hal vgp-577way-v1-clawed_frog.maf.gz --refGenome GCA_038501925.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-clawed_frog.maf.gz.log &
cactus-hal2maf jobstore/js-maf-brown_anole vgp-577way-v1.hal vgp-577way-v1-brown_anole.maf.gz --refGenome GCF_037176765.1 --outType norm single --index --coverage  --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile vgp-577way-v1-prep/logs/vgp-577way-v1-brown_anole.maf.gz.log &
wait
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



