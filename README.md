# vgp-cactus
Notes about making vgp phase-1 alignments

## Mammals

Note that the scripts below require pandas and biopython.  You can install them with

```
virtualenv -p python3 venv-vgp
. venv-vgp/bin/activate
pip install -U pandas biopython
```

Create the mammals seqfile (used mammals-v1.0 tag of this repo). Commands run in `./mammals`:

```
./gen-mammals-seqfile-v1.sh
```

Create the mammals workflow script (used cactus v3.0.0)

```
cactus-prepare mammals-v1.seqfile --outDir mammals-v1-prep --chromInfo mammals-v1.chrom-info --fastga --lastzCores 8 --alignCores 64 --cactusOptions '--batchSystem slurm --doubleMem true --slurmTime 100:00:00 --retryCount 5' --script --outHal mammals-v1.hal > mammals-v1.sh
```

