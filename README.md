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
cactus-prepare mammals-v1.seqfile --outDir mammals-v1-prep --chromInfo mammals-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal mammals-v1.hal > mammals-v1.sh
```

Note: needed to make the following manual changes:
* update to Toil `740508407ea34388f1a6f2d412964b67483ae9e8`
* add `--skipProgress` and `--maxMemory 1Ti` to cactus commands
* Toggle off FastGA options for: `Anc002`, `Anc039`, `Anc065`, `Anc102`, `Anc131`, `Anc145`

## Birds

Create the birds seqfile (used birds-v1.0 tag of this repo). Commands run in `./birds`:

```
./gen-birds-seqfile-v1.sh
```

Create the mammals workflow script (used cactus v3.0.1)

```
cactus-prepare birds-v1.seqfile --outDir birds-v1-prep --chromInfo birds-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5 --maxMemory 1Ti' --script --outHal birds-v1.hal > birds-v1.sh
```



