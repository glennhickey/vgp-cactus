#!/usr/bin/env python3
"""Generate a UCSC track hub for the vgp-577way phyloP conservation data set.

Walks a per-genome data directory (or, if it's unavailable, falls back to
parsing the cactus-phast invocations out of ../README.md), then emits the
hub.txt / genomes.txt / per-genome trackDb.txt + phyloP.html into a target
hub directory.

Each genome may have multiple .bw files — one per cactus-phast invocation,
distinguished by .r<root> / .s<subtree> tokens in the basename:
    vgp-577way-v1-<ref>.single.phyloP.bw                       full tree
    vgp-577way-v1-<ref>.single.phyloP.r<root>.bw               --root <root>
    vgp-577way-v1-<ref>.single.phyloP.s<subtree>.bw            --subtree <subtree>
    vgp-577way-v1-<ref>.single.phyloP.r<root>.s<subtree>.bw    both
All variants for a genome land as subtracks of one composite track.

Usage:
    build_hub.py [--hub-dir DIR] [--readme PATH] [--src-dir DIR]

Defaults assume the script lives in vgp-cactus/scripts/ and the hub
artifacts live in vgp-cactus/577way/hub/."""
import argparse
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_HUB  = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '577way', 'hub'))
DEFAULT_README = os.path.normpath(os.path.join(SCRIPT_DIR, '..', 'README.md'))
DEFAULT_SRC  = '/mnt/cluster_phast'

# UCSC-native database names used for both the `genome` declaration and
# the per-genome subdir for assemblies that have a long-established UCSC
# browser presence; everything else uses the GenArk accession directly.
NATIVE_DB = {
    'GCA_000001405.15': 'hg38',
    'GCA_009914755.4':  'hs1',
    'GCA_000001635.9':  'mm39',
}

# vgp-577way-v1-<commonName>.single.phyloP[.r<root>][.s<subtree>].bw
BW_RE = re.compile(
    r'^vgp-577way-v1-(?P<common>[A-Za-z0-9_]+)\.single\.phyloP'
    r'(?P<vtokens>(?:\.[rs][A-Za-z0-9_]+)*)\.bw$'
)


def parse_variant(bw_basename):
    """Return (common, root_or_None, subtree_or_None) for a cactus-phast
    bigwig basename, or None if the name doesn't fit."""
    m = BW_RE.match(bw_basename)
    if not m:
        return None
    common = m.group('common')
    root = subtree = None
    for tok in m.group('vtokens').split('.'):
        if not tok:
            continue
        if tok.startswith('r'):
            root = tok[1:]
        elif tok.startswith('s'):
            subtree = tok[1:]
    return common, root, subtree


def parse_phast_cmd(cmd):
    """Pull refGenome, output path, --root, --subtree values out of a
    cactus-phast command string."""
    toks = cmd.split()
    ref = toks[4] if len(toks) > 4 else None
    out = toks[5] if len(toks) > 5 else None
    root = None
    subtrees = []
    i = 0
    while i < len(toks):
        if toks[i] == '--root' and i + 1 < len(toks):
            root = toks[i + 1]
            i += 2
        elif toks[i] == '--subtree':
            i += 1
            while i < len(toks) and not toks[i].startswith('-'):
                subtrees.append(toks[i])
                i += 1
        else:
            i += 1
    return ref, out, root, subtrees


def cmd_produced_bws(cmd):
    """Given a cactus-phast command, return the basenames of every bigwig
    it would emit (mirrors cactus-phast's .r<root>/.s<subtree> output
    naming logic)."""
    ref, out, root, subtrees = parse_phast_cmd(cmd)
    if not out:
        return []
    base = os.path.basename(out)
    for ext in ('.wig.gz', '.wig', '.mod'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    if root:
        base = base + '.r' + root
    if subtrees:
        return [base + '.s' + s + '.bw' for s in subtrees]
    return [base + '.bw']


def read_readme_cmds(readme_path):
    """{bw_basename: cactus-phast-cmd}"""
    out = {}
    with open(readme_path) as f:
        for line in f:
            if not line.startswith('cactus-phast '):
                continue
            cleaned = line.strip().rstrip('&').rstrip()
            for bw in cmd_produced_bws(cleaned):
                out[bw] = cleaned
    return out


def enumerate_from_dir(src_dir, cmd_by_bw):
    """Walk src_dir/<accession>/ and collect every non-empty .bw."""
    genomes = []
    for acc in sorted(os.listdir(src_dir)):
        full = os.path.join(src_dir, acc)
        if not os.path.isdir(full):
            continue
        tracks = []
        for bw in sorted(os.listdir(full)):
            if not bw.endswith('.bw'):
                continue
            if os.path.getsize(os.path.join(full, bw)) == 0:
                continue
            info = parse_variant(bw)
            if not info:
                continue
            tracks.append((bw, info[1], info[2], cmd_by_bw.get(bw, '')))
        if not tracks:
            continue
        m = BW_RE.match(tracks[0][0])
        common = m.group('common') if m else acc
        genomes.append((acc, common, tracks))
    return genomes


def enumerate_from_readme(cmd_by_bw):
    """Build the genome list from the README cactus-phast lines alone
    (when no local data dir is available)."""
    by_acc = {}
    for bw, cmd in cmd_by_bw.items():
        ref, _, _, _ = parse_phast_cmd(cmd)
        info = parse_variant(bw)
        if not info or not ref:
            continue
        common = info[0]
        by_acc.setdefault(ref, (common, []))[1].append((bw, info[1], info[2], cmd))
    out = []
    for acc in sorted(by_acc):
        common, tracks = by_acc[acc]
        tracks.sort(key=lambda t: t[0])
        out.append((acc, common, tracks))
    return out


def variant_label(root, subtree):
    parts = []
    if root:    parts.append('clade ' + root)
    if subtree: parts.append('lineage ' + subtree)
    return ' (' + '; '.join(parts) + ')' if parts else ''


def variant_short(root, subtree):
    parts = []
    if root:    parts.append(root)
    if subtree: parts.append('s' + subtree)
    return '_'.join(parts) if parts else 'full'


def variant_track_id(root, subtree):
    s = ''
    if root:    s += 'R' + root
    if subtree: s += 'S' + subtree
    return s or 'Full'


def write_genome(hub_dir, acc, common, tracks):
    db = NATIVE_DB.get(acc, acc)
    gdir = os.path.join(hub_dir, db)
    os.makedirs(gdir, exist_ok=True)

    with open(os.path.join(gdir, 'trackDb.txt'), 'w') as f:
        f.write(
            'track vgp577wayPhyloP\n'
            'compositeTrack on\n'
            'shortLabel VGP 577-way phyloP\n'
            'longLabel VGP 577-way phyloP conservation (anchored at {common})\n'
            'type bigWig\n'
            'visibility full\n'
            'maxHeightPixels 100:50:11\n'
            'viewLimits -5.0:10.0\n'
            'viewLimitsMax -20:20\n'
            'autoScale off\n'
            'windowingFunction mean+whiskers\n'
            'color 25,25,95\n'
            'altColor 230,170,40\n'
            'html phyloP\n\n'.format(common=common))
        for i, (bw, root, subtree, _) in enumerate(tracks):
            tid = 'vgp577wayPhyloP{}'.format(variant_track_id(root, subtree))
            on_off = 'on' if i == 0 else 'off'
            short = 'full-tree' if not (root or subtree) else variant_short(root, subtree)
            f.write(
                '    track {tid}\n'
                '    parent vgp577wayPhyloP {on_off}\n'
                '    shortLabel {short}\n'
                '    longLabel VGP 577-way phyloP{lab} — {common}\n'
                '    bigDataUrl {bw}\n'
                '    type bigWig\n\n'.format(
                    tid=tid, on_off=on_off, short=short,
                    lab=variant_label(root, subtree), common=common, bw=bw))

    sections = []
    for bw, root, subtree, cmd in tracks:
        if root or subtree:
            heading = 'Variant: {}'.format(variant_label(root, subtree).strip(' ()'))
        else:
            heading = 'Variant: full 577-way tree'
        cmd_block = cmd if cmd else '(command not recorded; see https://github.com/glennhickey/vgp-cactus)'
        sections.append(
            '<h3>{h}</h3>\n<p>Bigwig: <code>{bw}</code></p>\n<pre>\n{c}\n</pre>\n'
            .format(h=heading, bw=bw, c=cmd_block))

    html = """<h2>VGP 577-way PhyloP conservation — {common} ({acc})</h2>

<p>Per-base PhyloP conservation scores from the
<b>vgp-577way-v1</b> Progressive Cactus alignment of 577 vertebrate
genomes, anchored at <code>{common}</code> ({acc}).</p>

<h3>Methods</h3>

<p>Computed with <code>cactus-phast</code> from
<a href="https://github.com/ComparativeGenomicsToolkit/cactus">Cactus</a>
v3.2.0. The reference's CDS annotation was used to extract 4-fold-degenerate
sites for neutral-model training; the trained model is the <code>REV</code>
substitution model fit by <code>phyloFit --EM --precision HIGH</code> with
<code>modFreqs</code>-symmetrized background frequencies; per-base scores
are <code>-log10(p)</code> from <code>phyloP --method LRT --mode CONACC</code>.
A <code>--root &lt;CladeAnc&gt;</code> variant restricts the entire pipeline
(model training and scoring) to the named clade; a
<code>--subtree &lt;CladeAnc&gt;</code> variant trains on the full tree but
performs the per-base LRT only on the named clade.</p>

{sections}
<p>Source: <a href="https://github.com/glennhickey/vgp-cactus">github.com/glennhickey/vgp-cactus</a> /
<a href="https://github.com/ComparativeGenomicsToolkit/cactus">github.com/ComparativeGenomicsToolkit/cactus</a></p>

<h3>Contact</h3>

<p>Glenn Hickey &lt;<a href="mailto:glenn.hickey@gmail.com">glenn.hickey@gmail.com</a>&gt;</p>
""".format(common=common, acc=acc, sections='\n'.join(sections))
    with open(os.path.join(gdir, 'phyloP.html'), 'w') as f:
        f.write(html)


def write_hub_root(hub_dir, genomes):
    # hub.txt only emitted if missing — preserves manual edits to email /
    # shortLabel / descriptionUrl across regenerations.
    hub_txt = os.path.join(hub_dir, 'hub.txt')
    if not os.path.exists(hub_txt):
        with open(hub_txt, 'w') as f:
            f.write(
                'hub vgp577wayConservation\n'
                'shortLabel Cactus VGP 577-way phyloP\n'
                'longLabel PhyloP Conservation for the vgp-577way-v1 Progressive Cactus alignment\n'
                'genomesFile genomes.txt\n'
                'descriptionUrl description.html\n'
                'email glenn.hickey@gmail.com\n')

    # genomes.txt is always regenerated from the current track list.
    with open(os.path.join(hub_dir, 'genomes.txt'), 'w') as f:
        for acc, _, _ in genomes:
            db = NATIVE_DB.get(acc, acc)
            f.write('genome {}\n'.format(db))
            f.write('trackDb {}/trackDb.txt\n\n'.format(db))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hub-dir', default=DEFAULT_HUB,
                    help='Output hub directory (default: ../577way/hub/)')
    ap.add_argument('--readme', default=DEFAULT_README,
                    help='vgp-cactus README.md with cactus-phast commands')
    ap.add_argument('--src-dir', default=DEFAULT_SRC,
                    help='Per-genome data dir (uses README-only fallback when absent)')
    args = ap.parse_args()

    cmd_by_bw = read_readme_cmds(args.readme)

    if os.path.isdir(args.src_dir) and os.listdir(args.src_dir):
        genomes = enumerate_from_dir(args.src_dir, cmd_by_bw)
    else:
        print('NOTE: {} not available; enumerating from README cmds only.'.format(args.src_dir),
              file=sys.stderr)
        genomes = enumerate_from_readme(cmd_by_bw)

    print('Producing hub for {} genomes:'.format(len(genomes)))
    for acc, common, tracks in genomes:
        print('  {} ({}) — {} track(s)'.format(acc, common, len(tracks)))
        for bw, root, subtree, _ in tracks:
            tag = []
            if root:    tag.append('--root ' + root)
            if subtree: tag.append('--subtree ' + subtree)
            print('      {}{}'.format(bw, '  [' + ', '.join(tag) + ']' if tag else ''))

    os.makedirs(args.hub_dir, exist_ok=True)
    write_hub_root(args.hub_dir, genomes)
    for acc, common, tracks in genomes:
        write_genome(args.hub_dir, acc, common, tracks)

    print('\nGenerated {} per-genome subdirs under {}'.format(len(genomes), args.hub_dir))


if __name__ == '__main__':
    main()
