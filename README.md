# vgp-cactus
Notes about making vgp phase-1 alignments

Note that the `./gen-xxx-seqfile-vx.sh` scripts below require pandas and biopython.  You can install them with

```
virtualenv -p python3 venv-vgp
. venv-vgp/bin/activate
pip install -U pandas biopython
```

## Amniotes

Create the amniotes (mammal/birds/reptile) seqfile. Commands run in `./amniotes`. Note that the branches in the ROADIES tree seem short, so this script scales them by 2 to help cactus select appropriate lastz parameters. 

```
./gen-amniotes-seqfile-v1.sh
```

Create the amniotes workflow script (used cactus v3.1.0)
```
cactus-prepare amniotes-v1.seqfile --outDir amniotes-v1-prep --chromInfo amniotes-v1.chrom-info --gpu 4 --lastzCores 40 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal vgp-amniotes-v1.hal --remask --preprocessBatchSize 1000 > amniotes-v1.sh
chmod +x amniotes-v1.sh
sed -i amniotes-v1.sh -e 's/--gpu/--binariesMode docker --gpu/g'
```

Don't align or include the outgroups (which are under Anc002 and Anc000)
```
sed -i amniotes-v1.sh -e '/^cactus-halAppendSubtrees/s|amniotes-v1-prep/Anc000\.hal||g' -e '/^cactus-halAppendSubtrees/s|amniotes-v1-prep/Anc002\.hal||g'
sed -i amniotes-v1.sh -e '/Anc000/d' -e '/Anc002/d'
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
Rename Homo_sapiens to hs1 to be a little more clear about it (since hg38 is also in the alignment)
```
printf "Anc001\tAmniotesAnc000\n" > amniotes-rename.tsv
for i in `seq 3 360`; do printf "Anc%03d\tAmniotesAnc%03d\n" $i $((i-2)); done >> amniotes-rename.tsv
printf "Homo_sapiens\ths1\n" >> amniotes-rename.tsv
halRenameGenomes vgp-amniotes-v1.hal amniotes-rename.tsv
```

Export the TAFs with

```
cactus-hal2maf jobstore/js-maf-hg38 vgp-amniotes-v1.hal vgp-amniotes-v1-hg38.taf.gz --refGenome hg38 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-hg38.taf.gz.log & 
cactus-hal2maf jobstore/js-maf-mm39 vgp-amniotes-v1.hal vgp-amniotes-v1-mm39.taf.gz --refGenome mm39 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-mm39.taf.gz.log &
cactus-hal2maf jobstore/js-maf-hs1 vgp-amniotes-v1.hal vgp-amniotes-v1-hs1.taf.gz --refGenome hs1 --outType norm single --index --coverage --coverageSexChroms chrX chrY --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-hs1.taf.gz.log &
cactus-hal2maf jobstore/js-maf-chicken vgp-amniotes-v1.hal vgp-amniotes-v1-chicken.taf.gz --refGenome Gallus_gallus --outType norm single --index --coverage --coverageSexChroms NC_059535.1 NC_059536.1 --chunkSize 100000 --batchCores 96 --batchCount 32 --noAncestors --batchParallelTaf 32 --batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti --disableProgress --logFile amniotes-v1-prep/logs/vgp-amniotes-v1-chicken.taf.gz.log &

wait
```

Make MAFs using the accessions for the browser
```
cat amniotes-v1.seqfile  | grep -v Anc | sed -e 's#amniotes-v1-prep/##g' -e 's/.fa.gz//g' | grep -v mm39 | grep -v hg38 | grep -v Homo_sapiens | grep -v hs1 | grep -v '^$' | awk -F'\t' '{if ($2 ~ /^https?:\/\//) {n=split($2,a,"/"); $2=a[n]} print}' OFS='\t'> vgp-amniotes-v1.accessions

taffy view -i vgp-amniotes-v1-hg38.single.taf.gz -n vgp-amniotes-v1.accessions -m | mafFilter -m - -l 2 | bgzip > vgp-amniotes-v1-hg38.single.maf.gz &
taffy view -i vgp-amniotes-v1-hs1.single.taf.gz -n vgp-amniotes-v1.accessions -m | mafFilter -m - -l 2 | bgzip > vgp-amniotes-v1-hs1.single.maf.gz &
taffy view -i vgp-amniotes-v1-mm39.single.taf.gz -n vgp-amniotes-v1.accessions -m | mafFilter -m - -l 2 | bgzip > vgp-amniotes-v1-mm39.single.maf.gz &
taffy view -i vgp-amniotes-v1-chicken.single.taf.gz -n vgp-amniotes-v1.accessions -m | mafFilter -m - -l 2 | bgzip > vgp-amniotes-v1-chicken.single.maf.gz &
wait
```

And the chains with (note: these have to be submitted via slurm!)
```
cactus-hal2chains jobstore/js-chain-hs1 vgp-mammals-v1.hal vgp-mammals-v1-chains-hs1 --targetGenomes hs1 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-hs1.taf.gz.log
cactus-hal2chains jobstore/js-chain-hg38 vgp-mammals-v1.hal vgp-mammals-v1-chains-hg38 --targetGenomes hg38 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-hg38.taf.gz.log
cactus-hal2chains jobstore/js-chain-mm39 vgp-mammals-v1.hal vgp-mammals-v1-chains-mm39 --targetGenomes mm39 --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-mm39.taf.gz.log
cactus-hal2chains jobstore/js-chain-chicken vgp-mammals-v1.hal vgp-mammals-v1-chains-chicken --targetGenomes Gallus_gallus --bigChain --retryCount 5 --maxMemory 1Ti --logFile mammals-v1-prep/logs/vgp-mammals-v1-chicken.taf.gz.log 
```


## Rayfin-Fish

Create the rayfin-fish seqfile. Commands run in `./rayfin-fish`:

```
./gen-rayfin-fish-seqfile-v1.sh
```

Create the rayfin-fish workflow script (used cactus v3.0.1)

```
cactus-prepare rayfin-fish-v1.seqfile --outDir rayfin-fish-v1-prep --chromInfo rayfin-fish-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal rayfin-fish-v1.hal > rayfin-fish-v1.sh
```

We don't want the outgroups, so pull them out of the final alignment
```
sed -i rayfin-fish-v1.sh -e 's#rayfin-fish-v1-prep/Anc000.hal rayfin-fish-v1-prep/Anc001.hal ##g'
```



