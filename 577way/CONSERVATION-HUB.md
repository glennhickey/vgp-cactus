# 577-way phyloP track hub

This document explains how the
[VGP 577-way phyloP conservation hub](https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/conservation/hub.txt)
is built and how to regenerate / extend it.

## What's in the hub

For each reference where we had a CDS gene annotation, we ran
`cactus-phast` on the single-coverage MAF anchored at that reference
(see the [cactus-phast recipes in ../README.md](../README.md#577-way))
and uploaded the resulting bigwig (plus the trained model, 4d-site SS
file, name-map sed script, and a gzipped wig copy for non-browser
tooling) under

```
s3://genomeark/downstream_analyses/genome_alignments/cactus/577way/conservation/<ref>/
```

`<ref>` is the **UCSC database name** for assemblies that have a long-
established browser presence (`hg38`, `hs1`, `mm39`) and the **GenArk
accession** otherwise (`GCF_902713615.1`, `GCA_944039275.1`, …). The
hub's `genome <db>` declaration in `genomes.txt` uses the same
identifier so the UCSC browser lands the tracks in the right assembly
context — for non-UCSC-native dbs, the GenArk auto-resolver fetches the
assembly hub from
[hgdownload.soe.ucsc.edu/hubs/](https://hgdownload.soe.ucsc.edu/hubs/).

## Track design

- One **composite track** per genome (`track vgp577wayPhyloP`), with
  one subtrack per `cactus-phast` invocation. A genome may produce
  multiple invocations:
  - **Full-tree** (no `--root` / `--subtree`): scores against all 577
    species. On by default in the browser.
  - **`--root <CladeAnc>`**: restricts the entire pipeline (model
    training and scoring) to the named clade. Output filename gets a
    `.r<root>` token (`vgp-577way-v1-zebra_finch.single.phyloP.rBirdsAnc0.bw`).
  - **`--subtree <CladeAnc>`**: trains on the full tree but performs
    the per-base LRT only on the named clade. Output gets a `.s<subtree>`
    token. Both flags can be combined.
- **Color**: `25,25,95` (dark blue) for positive scores =
  conservation, `230,170,40` (orange) for negative = acceleration.
  Matches UCSC's standard phyloP track convention.
- **`viewLimits -5.0:10.0`** as the default zoom; **`viewLimitsMax
  -20:20`** caps the user-adjustable slider. Phast hard-clips the
  `-log10(p)` it emits to ±20 so the slider exactly covers the value
  domain.
- **`windowingFunction mean+whiskers`** so zoomed-out views show a
  min/max envelope plus mean — the natural choice for signed-signal
  tracks.

## How to (re)build the hub

The script that generates `hub.txt` / `genomes.txt` / per-genome
`trackDb.txt` + `phyloP.html` is at
[../scripts/build_hub.py](../scripts/build_hub.py). It works by:

1. Parsing the `cactus-phast` invocations out of `../README.md`. Each
   command's positional args + `--root` / `--subtree` flags determine
   the bigwig basename it produces (mirroring `cactus-phast`'s own
   output-naming logic).
2. Walking `/mnt/cluster_phast/<accession>/` (the sshfs mount of the
   cluster's conservation output) to enumerate which bigwigs actually
   exist; if that's unavailable, falling back to the README-only
   listing.
3. Emitting one trackDb stanza per bigwig, grouped under one composite
   per genome. The first (alphabetically) is `on` by default; the
   rest are `off` so they don't auto-stack.

```
./scripts/build_hub.py                  # writes into 577way/hub/
./scripts/build_hub.py --hub-dir /tmp/x # custom output location
```

The script is **idempotent** — safe to rerun any time a new track
comes online. It preserves `577way/hub/hub.txt` if it already exists
(so manual edits to `shortLabel` / `email` / `descriptionUrl` aren't
overwritten); `genomes.txt` and the per-genome dirs are always
regenerated from scratch.

## How to add a new track

1. Run `cactus-phast` per the recipe in [../README.md](../README.md).
   Output lands at `conservation/<accession>/vgp-577way-v1-<ref>.single.phyloP.*`.
2. Add the `cactus-phast` invocation to `../README.md` so it's
   reproducible. (The hub builder reads this README, so the new track
   won't show up in the hub without this step.)
3. Sync the new data to s3:
   ```
   S3=s3://genomeark/downstream_analyses/genome_alignments/cactus/577way/conservation
   aws s3 sync /private/.../conservation/<accession>/ "$S3/<accession>/"
   ```
4. If `<accession>` is in `NATIVE_DB` (hg38 / hs1 / mm39), rename the
   s3 dir:
   ```
   aws s3 mv "$S3/<accession>/" "$S3/<db>/" --recursive --copy-props none
   ```
   `--copy-props none` skips the multipart-copy's
   `GetObjectTagging` call, which our IAM user can't make.
5. Regenerate the hub locally and resync just the trackDb files:
   ```
   ./scripts/build_hub.py
   aws s3 sync 577way/hub/ "$S3/" --exclude "build_hub.py"
   ```

The local hub regeneration step picks the new track up automatically
as long as its `cactus-phast` invocation is in the README and its
bigwig is non-empty under `/mnt/cluster_phast/` (or simply visible in
the README, via the fallback path).

## How to add a new variant (`--root` / `--subtree`) for an existing genome

Same as adding a new track, but you'll generate a second bigwig in the
same `conservation/<accession>/` dir; the builder picks both up and
turns the genome's trackDb into a composite with two (or more)
subtracks.

Example: zebra finch has both
`vgp-577way-v1-zebra_finch.single.phyloP.bw` (full tree) and
`vgp-577way-v1-zebra_finch.single.phyloP.rBirdsAnc0.bw` (Bird-clade
only). See the corresponding commands in [../README.md](../README.md).

## Special case: GCA-vs-GCF assemblies

If the alignment contains a GCA accession but the only published CDS
annotation is on the paired GCF version (same physical contigs,
different chrom-name conventions), use the renamer:

```
./scripts/rename-gff-refseq-to-genbank.sh <GCF_acc> <GCA_acc> <out.gff.gz>
```

It downloads both the GCF GFF and the GCA `assembly_report.txt`,
builds a RefSeq → GenBank mapping for the chromosomes, and rewrites
the GFF chrom column in a single pass. The output is suitable as
`--geneAnnotation` for the GCA-anchored `cactus-phast` run.

Example: zebra finch (`GCF_048771995.1` / `GCA_048771995.1`), used in
the recipe in [../README.md](../README.md#zebra-finch-special-case).

## Validation

After uploading, sanity-check that the hub is well-formed and every
bigwig is reachable:

```
URL=https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/conservation
hubCheck "$URL/hub.txt"                                    # exit 0
curl -sSI "$URL/hub.txt" | head -1                         # 200
for g in $(curl -sS "$URL/genomes.txt" | awk '/^genome / {print $2}'); do
    curl -sSI "$URL/$g/trackDb.txt"          | head -1 | sed "s|^|$g |"
done
```

Per-bigwig header sanity (header-only fetch over https; doesn't
download data sections):

```
bigWigInfo "$URL/hg38/vgp-577way-v1-hg38.single.phyloP.bw" \
    | grep -E "^(min|max|mean|std|chromCount|basesCovered):"
```

Expect `min: -20, max: 20` (phast cap), `mean ~ 0.3`, `std ~ 1.5–2.5`
for full-tree tracks. `--root <CladeAnc>` tracks generally have a
lower max (less power on smaller subtrees).

## Source-of-truth references

- Hub: <https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/conservation/hub.txt>
- Top-level README on s3: <https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/README.html>
- Cactus-phast docs: <https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#conservation-scores>
- Per-track command lines: [../README.md](../README.md)
- PHAST package: <http://compgen.cshl.edu/phast/>
