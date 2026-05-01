#!/usr/bin/env python3
"""Summarize cactus pipeline logs (blast/align/halAppend/preprocess) into a
high-level breakdown by command category, plus per-phase timings for
cactus_consolidated.

Usage:
    ./analyze_cactus_logs.py <logs_dir> [--out-dir DIR] [--prefix STR]

Produces:
    <out>/<prefix>command_summary.png   bar chart: count / total walltime / peak memory per category
    <out>/<prefix>command_summary.tsv   raw aggregated numbers
    <out>/<prefix>consolidated_phases.png   stacked bars of caf/bar/reference/hal duration per genome
    <out>/<prefix>consolidated_phases.tsv   raw per-genome phase durations
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Regex for the standard "Successfully ran" line. Two flavors:
#   Successfully ran: "<cmd>" in <secs> seconds and <mem> <unit> memory
#   Successfully ran cactus_consolidated(<event>): "<cmd>" in <secs> seconds and <mem> <unit> memory ...
SUCC_RE = re.compile(
    r'Successfully ran(?:\((?P<scope1>[^)]+)\))?(?: (?P<scope2>cactus_consolidated\([^)]+\)))?: "(?P<cmd>.*?)" in (?P<secs>[0-9.]+) seconds and (?P<mem>[0-9.]+) (?P<unit>[KMGTP]i?) memory'
)

# Phase marker line for cactus_consolidated:
#   cactus_consolidated(<event>): <description>, <secs> seconds have elapsed
PHASE_RE = re.compile(
    r'cactus_consolidated\((?P<event>[^)]+)\): (?P<desc>.*?), (?P<secs>[0-9]+) seconds have elapsed'
)

UNIT_TO_BYTES = {
    "K": 1e3, "Ki": 2 ** 10,
    "M": 1e6, "Mi": 2 ** 20,
    "G": 1e9, "Gi": 2 ** 30,
    "T": 1e12, "Ti": 2 ** 40,
    "P": 1e15, "Pi": 2 ** 50,
}


# Match the two FASTA arguments at the start of a lastz command line:
#   lastz <target_genome>_<chunk>.fa[...] <query_genome>_<chunk>.fa[...] ...
# Genome name may itself contain underscores (e.g. GCA_028022735.1).
LASTZ_GENOMES_RE = re.compile(
    r"^lastz\s+'?(\S+?)_\d+\.fa\S*?'?\s+'?(\S+?)_\d+\.fa"
)


def lastz_genomes(cmd):
    """Return (target, query) genome names for a lastz command, or None."""
    m = LASTZ_GENOMES_RE.match(cmd)
    if not m:
        return None
    return m.group(1), m.group(2)


# Differentiating params (after dropping --queryhspbest=100000 and
# --ambiguous=iupac,100,100 which are common to every class) for each
# lastz parameter class defined in cactus_progressive_config.xml.
LASTZ_CLASSES = {
    "one":     "--step=2 --ydrop=3000 --notransition",
    "two":     "--step=5 --ydrop=3000",
    "three":   "--step=4 --ydrop=3500 --hspthresh=2800",
    "four":    "--step=3 --ydrop=3500 --hspthresh=2600 --gappedthresh=2800",
    "five":    "--step=2 --ydrop=4000 --hspthresh=2400 --gappedthresh=2600",
    "default": "--step=1 --ydrop=4000 --hspthresh=2200 --gappedthresh=2400",
}
LASTZ_CLASS_ORDER = ["one", "two", "three", "four", "five", "default"]
LASTZ_CLASS_COLORS = {
    "one":     "#1f77b4",
    "two":     "#2ca02c",
    "three":  "#ff7f0e",
    "four":    "#d62728",
    "five":    "#9467bd",
    "default": "#7f7f7f",
}


def classify_lastz(cmd):
    """Map a lastz command line to one of the 6 cactus parameter classes
    by inspecting --step / --ydrop / --hspthresh / --notransition.
    """
    step = re.search(r"--step=(\d+)", cmd)
    ydrop = re.search(r"--ydrop=(\d+)", cmd)
    hsp = re.search(r"--hspthresh=(\d+)", cmd)
    has_notrans = "--notransition" in cmd
    s = int(step.group(1)) if step else None
    y = int(ydrop.group(1)) if ydrop else None
    h = int(hsp.group(1)) if hsp else None
    if s == 2 and has_notrans:
        return "one"
    if s == 5:
        return "two"
    if s == 4:
        return "three"
    if s == 3:
        return "four"
    if s == 2 and not has_notrans:
        return "five"
    if s == 1:
        return "default"
    return "unknown"


def _shell_head(cmd):
    """Extract the dominant program name from a `bash -c '...'` command line.
    Strips the leading shell preamble (set -eo pipefail &&), then returns the
    first verb of the resulting command. Words inside `pipefail` etc. should
    NOT be matched as the verb."""
    # Drop the bash -c '...' wrapper if present
    m = re.match(r"^bash -c '(.*)'$", cmd, flags=re.DOTALL)
    inner = m.group(1) if m else cmd
    # Drop a leading shell preamble like `set -eo pipefail &&` or just `set -e &&`
    inner = re.sub(r"^\s*set -[a-z]+(?:\s+[a-z]+)?\s*&&\s*", "", inner)
    # Take the first whitespace-separated token
    tokens = inner.split()
    return tokens[0] if tokens else ""


def categorize(cmd):
    """Bucket a command-line into one of our coarse categories."""
    head = cmd.lstrip()
    if head.startswith("bash -c"):
        head = _shell_head(head)
    first = head.split()[0] if head.split() else ""
    base = os.path.basename(first)

    if base.startswith("lastz"):
        return "lastz"
    if base.startswith("cactus_consolidated"):
        return "cactus_consolidated"
    if base == "halAppendSubtree":
        return "halAppendSubtree"
    if base.startswith("paffy"):
        return "paffy"
    if base == "Red" or base.startswith("cactus_red"):
        return "Red"
    return "Other"


def parse_log(path):
    """Extract (command, secs, mem_bytes, category) tuples and per-event
    phase elapsed dicts from one log file.

    Returns (commands, phases) where:
      commands: list of dicts {cmd, category, secs, mem_bytes}
      phases:   {event: [(desc, secs), ...]} for cactus_consolidated phases
    """
    commands = []
    phases = defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            m = SUCC_RE.search(line)
            if m:
                cmd = m.group("cmd")
                secs = float(m.group("secs"))
                mem_bytes = float(m.group("mem")) * UNIT_TO_BYTES[m.group("unit")]
                # If this is the cactus_consolidated wrap line, force the
                # category regardless of what the inner command looks like.
                if m.group("scope2") and m.group("scope2").startswith("cactus_consolidated"):
                    cat = "cactus_consolidated"
                else:
                    cat = categorize(cmd)
                # Extract --threads N from the cmd; default to 1 (most cactus
                # commands run single-core; cactus_consolidated is multi-core).
                tm = re.search(r"--threads\s+(\d+)", cmd)
                threads = int(tm.group(1)) if tm else 1
                commands.append({"cmd": cmd, "category": cat,
                                 "secs": secs, "mem_bytes": mem_bytes,
                                 "threads": threads})
                continue
            mp = PHASE_RE.search(line)
            if mp:
                phases[mp.group("event")].append(
                    (mp.group("desc"), int(mp.group("secs")))
                )
    return commands, phases


# Phase grouping. Each entry: (group, marker_substring) where marker_substring
# matches the description ending the group (the cumulative elapsed value at
# this marker is the end-of-phase time).
PHASE_GROUPS = [
    ("caf",       "Ran cactus caf"),
    ("bar",       "Ran cactus bar"),
    ("reference", "Ran cactus make reference top down coordinates"),
    ("hal",       "Cactus consolidated is done!"),
]


def phase_durations(markers):
    """Convert a list of (desc, cumulative_secs) markers into per-phase
    durations. Returns dict {caf, bar, reference, hal} -> seconds. Missing
    markers leave entries at 0.
    """
    # Build cumulative end-of-phase times by scanning markers in order.
    ends = {}
    for desc, secs in markers:
        for group, needle in PHASE_GROUPS:
            if needle in desc:
                ends[group] = secs
    durations = {}
    prev = 0
    for group, _ in PHASE_GROUPS:
        if group in ends:
            durations[group] = ends[group] - prev
            prev = ends[group]
        else:
            durations[group] = 0
    return durations


def fmt_hours(secs):
    return secs / 3600.0


def fmt_gib(b):
    return b / (2 ** 30)


def plot_lastz_classes(lastz_by_class, out_dir, prefix, log_scale=True,
                       suffix=""):
    """Box-plot of lastz job runtimes per parameter class, with a legend
    showing the differentiating arguments for each class.
    """
    classes_present = [c for c in LASTZ_CLASS_ORDER if lastz_by_class.get(c)]
    if "unknown" in lastz_by_class:
        classes_present.append("unknown")
    data = [np.array(lastz_by_class[c]) / 60.0 for c in classes_present]  # minutes
    counts = [len(lastz_by_class[c]) for c in classes_present]
    colors = [LASTZ_CLASS_COLORS.get(c, "#bbbbbb") for c in classes_present]

    fig, (ax_box, ax_legend) = plt.subplots(
        1, 2, figsize=(13, 5.5),
        gridspec_kw={"width_ratios": [3, 2]})

    positions = np.arange(1, len(data) + 1)
    if log_scale:
        # Compute violins in log-space so the KDE shapes are sensible across
        # 4+ decades; we'll relabel the axis in original minutes below.
        plot_data = [np.log10(np.clip(arr, 1e-3, None)) for arr in data]
    else:
        plot_data = data
    vp = ax_box.violinplot(plot_data, positions=positions, showmeans=False,
                           showmedians=True, showextrema=True, widths=0.8)
    for body, color in zip(vp["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(0.85)
    for part_name in ("cbars", "cmins", "cmaxes", "cmedians"):
        if part_name in vp:
            vp[part_name].set_color("black")
            vp[part_name].set_linewidth(0.8)
    ax_box.set_xticks(positions)
    ax_box.set_xticklabels(classes_present)
    if log_scale:
        # Y-axis is log10(minutes); show ticks at decade boundaries.
        all_log = np.concatenate(plot_data)
        lo = np.floor(all_log.min())
        hi = np.ceil(all_log.max())
        ticks = np.arange(lo, hi + 1)
        ax_box.set_yticks(ticks)
        ax_box.set_yticklabels(
            [f"{10 ** t:g}" if t >= 0 else f"{10 ** t:.3g}" for t in ticks]
        )
        ax_box.set_ylabel("walltime (minutes, log)")
    else:
        ax_box.set_ylabel("walltime (minutes)")
        # On a linear scale extreme outliers crush the violins; clip the y-axis
        # to a robust upper bound based on the slowest class' tail.
        upper = max(np.percentile(arr, 99) for arr in data) * 1.1
        ax_box.set_ylim(0, upper)
    ax_box.set_xlabel("lastz parameter class")
    total = sum(counts)
    ax_box.set_title(f"lastz runtime by parameter class "
                     f"(n={total:,} jobs)")
    # n labels above each box
    for i, n in enumerate(counts, 1):
        ax_box.text(i, ax_box.get_ylim()[1], f"n={n:,}",
                    ha="center", va="top", fontsize=8, color="#444")

    # Legend panel: class name + differentiating params
    ax_legend.axis("off")
    ax_legend.set_title("class → differentiating params\n"
                        "(common to all: --ambiguous=iupac,100,100 "
                        "--queryhspbest=100000)",
                        fontsize=10, loc="left")
    y = 0.92
    for c in classes_present:
        if c == "unknown":
            params = "(could not classify)"
        else:
            params = LASTZ_CLASSES[c]
        ax_legend.add_patch(plt.Rectangle((0.02, y - 0.02), 0.05, 0.04,
                                          facecolor=LASTZ_CLASS_COLORS.get(c, "#bbbbbb"),
                                          transform=ax_legend.transAxes,
                                          clip_on=False))
        ax_legend.text(0.10, y, c, fontsize=11, fontweight="bold",
                       transform=ax_legend.transAxes, va="center")
        ax_legend.text(0.22, y, params, fontsize=9, family="monospace",
                       transform=ax_legend.transAxes, va="center")
        y -= 0.13

    fig.tight_layout()
    out = os.path.join(out_dir, f"{prefix}lastz_by_class{suffix}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}", file=sys.stderr)

    # Also write a TSV with summary stats per class
    tsv = os.path.join(out_dir, f"{prefix}lastz_by_class{suffix}.tsv")
    with open(tsv, "w") as fh:
        fh.write("class\tcount\tmean_sec\tmedian_sec\tp95_sec\tmax_sec\ttotal_h\n")
        for c, vals in [(c, lastz_by_class[c]) for c in classes_present]:
            arr = np.array(vals)
            fh.write(f"{c}\t{len(arr)}\t{arr.mean():.2f}\t"
                     f"{np.median(arr):.2f}\t{np.percentile(arr, 95):.2f}\t"
                     f"{arr.max():.2f}\t{arr.sum() / 3600.0:.2f}\n")
    print(f"Wrote {tsv}", file=sys.stderr)


def load_name_map(annotations_path, extra_scinames_paths):
    """Build {accession: display_name}, preferring EnglishName over
    ScientificName. Reads a TSV whose first row is a header (with `# ` prefix
    or not) containing 'ScientificName' and optionally 'EnglishName' columns.
    Then layers any extra accession\\tScientificName files on top, only
    filling in entries not already present.
    """
    name_map = {}

    def add(accession, scientific=None, english=None):
        if not accession:
            return
        if english:
            name_map.setdefault(accession, english)
        elif scientific:
            name_map.setdefault(accession, scientific.replace("_", " ")
                                .replace(" ", "_"))

    if annotations_path and os.path.isfile(annotations_path):
        with open(annotations_path) as fh:
            header = fh.readline().lstrip("# ").rstrip("\n").split("\t")
            try:
                acc_idx = header.index("accession")
            except ValueError:
                acc_idx = 0
            sci_idx = header.index("ScientificName") if "ScientificName" in header else None
            eng_idx = header.index("EnglishName") if "EnglishName" in header else None
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= acc_idx:
                    continue
                acc = cols[acc_idx]
                sci = cols[sci_idx] if sci_idx is not None and sci_idx < len(cols) else ""
                eng = cols[eng_idx] if eng_idx is not None and eng_idx < len(cols) else ""
                add(acc, sci or None, eng or None)

    for path in extra_scinames_paths:
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    add(cols[0], scientific=cols[1])

    return name_map


def display_name(label, name_map):
    """Return the English/scientific name if `label` looks like a GC[AF]_
    accession that we have a mapping for; otherwise return the label
    unchanged (ancestor names, anything else)."""
    if label in name_map:
        return name_map[label]
    return label


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs_dir")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--prefix", default="")
    p.add_argument("--top", type=int, default=60,
                   help="show top-N genomes (by walltime) in per-genome phase plot; "
                        "0 = show all")
    p.add_argument("--annotations",
                   default=os.path.join(
                       os.path.dirname(os.path.abspath(__file__)),
                       "..", "mammals", "coverage-2", "annotations.tsv"),
                   help="annotations TSV (e.g. mammals/coverage-2/annotations.tsv) "
                        "with at least accession+ScientificName, ideally also EnglishName")
    p.add_argument("--scinames", action="append", default=[],
                   help="extra accession\\tScientificName TSV; can be repeated")
    p.add_argument("--exclude-stages", default="vgp",
                   help="comma-separated log filename prefixes to skip "
                        "(default: 'vgp', which excludes the chains and "
                        "maf-export logs)")
    args = p.parse_args()
    excluded_stages = set(s for s in args.exclude_stages.split(",") if s)

    name_map = load_name_map(args.annotations, args.scinames)
    print(f"Loaded {len(name_map)} accession name mappings", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)

    # Collect log files. Group by their prefix (blast/align/halAppend/preprocess
    # /vgp/hal2fasta) so we can report per-stage totals as well.
    log_files = sorted(f for f in os.listdir(args.logs_dir)
                       if f.endswith(".log"))
    by_stage = defaultdict(list)
    for f in log_files:
        m = re.match(r"([a-zA-Z0-9]+?)-", f)
        stage = m.group(1) if m else "other"
        if stage in excluded_stages:
            continue
        by_stage[stage].append(os.path.join(args.logs_dir, f))

    print(f"Using {sum(len(v) for v in by_stage.values())} of {len(log_files)} "
          f"logs across stages: { {k: len(v) for k, v in by_stage.items()} }"
          + (f" (excluded: {sorted(excluded_stages)})" if excluded_stages else ""),
          file=sys.stderr)

    # Per-(stage, category) aggregates
    agg = defaultdict(lambda: {"count": 0, "secs": 0.0, "cpu_secs": 0.0,
                               "peak_mem_b": 0.0})
    # Per-genome cactus_consolidated phase durations (only meaningful in align logs)
    phase_table = {}
    # Per-genome cactus_consolidated wall-time and peak memory
    consolidated_meta = {}
    # Per-class lastz runtimes (seconds), pooled across all stages
    lastz_by_class = defaultdict(list)
    # Per-genome total lastz seconds (each lastz job contributes its runtime
    # to BOTH the target and query genomes it operated on)
    lastz_by_genome = defaultdict(float)

    for stage, paths in by_stage.items():
        for path in paths:
            commands, phases = parse_log(path)
            for c in commands:
                key = (stage, c["category"])
                agg[key]["count"] += 1
                agg[key]["secs"] += c["secs"]
                agg[key]["cpu_secs"] += c["secs"] * c["threads"]
                agg[key]["peak_mem_b"] = max(agg[key]["peak_mem_b"], c["mem_bytes"])
                if c["category"] == "lastz":
                    lastz_by_class[classify_lastz(c["cmd"])].append(c["secs"])
                    pair = lastz_genomes(c["cmd"])
                    if pair:
                        lastz_by_genome[pair[0]] += c["secs"]
                        if pair[1] != pair[0]:
                            lastz_by_genome[pair[1]] += c["secs"]
            # Only treat phases from align- logs as belonging to that genome run
            if stage == "align":
                for event, markers in phases.items():
                    phase_table[event] = phase_durations(markers)
                # Recover this align log's event name and find its
                # cactus_consolidated runtime + memory by scanning the entries.
                event_name = re.match(r"align-(.+)\.log$",
                                      os.path.basename(path))
                if event_name:
                    ev = event_name.group(1)
                    for c in commands:
                        if c["category"] == "cactus_consolidated":
                            consolidated_meta[ev] = {
                                "secs": c["secs"],
                                "mem_bytes": c["mem_bytes"],
                            }
                            break

    # ---- Write command summary TSV ----
    summary_tsv = os.path.join(args.out_dir, f"{args.prefix}command_summary.tsv")
    categories = ["lastz", "cactus_consolidated", "halAppendSubtree",
                  "paffy", "Red", "Other"]
    stages = sorted(by_stage.keys())
    with open(summary_tsv, "w") as fh:
        fh.write("stage\tcategory\tcount\ttotal_walltime_h\ttotal_cpu_h\t"
                 "peak_mem_GiB\n")
        for stage in stages:
            for cat in categories:
                a = agg.get((stage, cat))
                if not a or a["count"] == 0:
                    continue
                fh.write(f"{stage}\t{cat}\t{a['count']}\t"
                         f"{fmt_hours(a['secs']):.3f}\t"
                         f"{fmt_hours(a['cpu_secs']):.3f}\t"
                         f"{fmt_gib(a['peak_mem_b']):.2f}\n")
    print(f"Wrote {summary_tsv}", file=sys.stderr)

    # ---- Plot command summary figure ----
    # We'll show 3 panels: count, total walltime (hours), peak memory (GiB),
    # with bars for each category, summed across stages.
    cat_totals = {c: {"count": 0, "secs": 0.0, "cpu_secs": 0.0,
                      "peak_mem_b": 0.0}
                  for c in categories}
    for (stage, cat), a in agg.items():
        cat_totals[cat]["count"] += a["count"]
        cat_totals[cat]["secs"] += a["secs"]
        cat_totals[cat]["cpu_secs"] += a["cpu_secs"]
        cat_totals[cat]["peak_mem_b"] = max(cat_totals[cat]["peak_mem_b"],
                                            a["peak_mem_b"])

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes = axes.flatten()
    cats_present = [c for c in categories if cat_totals[c]["count"] > 0]
    counts = [cat_totals[c]["count"] for c in cats_present]
    hours = [fmt_hours(cat_totals[c]["secs"]) for c in cats_present]
    cpu_hours = [fmt_hours(cat_totals[c]["cpu_secs"]) for c in cats_present]
    mems = [fmt_gib(cat_totals[c]["peak_mem_b"]) for c in cats_present]

    color_map = {
        "lastz": "#1f77b4",
        "cactus_consolidated": "#d62728",
        "halAppendSubtree": "#9467bd",
        "paffy": "#2ca02c",
        "Red": "#ff7f0e",
        "Other": "#7f7f7f",
    }
    colors = [color_map[c] for c in cats_present]

    for ax, vals, title, ylabel, log in [
        (axes[0], counts, "Command count", "count", True),
        (axes[1], hours, "Total walltime", "hours (sum of per-job walltime)", True),
        (axes[2], cpu_hours, "Total CPU time", "core-hours (walltime × threads)", True),
        (axes[3], mems, "Peak memory (single job)", "GiB", False),
    ]:
        bars = ax.bar(cats_present, vals, color=colors)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        if log and max(vals) / max(min(vals), 1e-9) > 100:
            ax.set_yscale("log")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                    f"{v:.1f}" if v < 1000 else f"{v:.0f}",
                    ha="center", va="bottom", fontsize=9)
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")

    n_used = sum(len(v) for v in by_stage.values())
    excl_note = (f"; excluded stages: {sorted(excluded_stages)}"
                 if excluded_stages else "")
    fig.suptitle(f"Cactus log summary ({sum(counts):,} commands across "
                 f"{n_used} log files{excl_note})", fontsize=12)
    fig.tight_layout()
    summary_png = os.path.join(args.out_dir, f"{args.prefix}command_summary.png")
    fig.savefig(summary_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {summary_png}", file=sys.stderr)

    # ---- lastz runtime distribution by parameter class ----
    if lastz_by_class:
        plot_lastz_classes(lastz_by_class, args.out_dir, args.prefix,
                           log_scale=True, suffix="")
        plot_lastz_classes(lastz_by_class, args.out_dir, args.prefix,
                           log_scale=False, suffix="_linear")

    # ---- per-genome lastz total walltime ----
    if lastz_by_genome:
        plot_lastz_per_genome(lastz_by_genome, args.out_dir, args.prefix,
                              args.top, name_map)

    # ---- cactus_consolidated phase analysis ----
    if not phase_table:
        print("No cactus_consolidated phase markers found; skipping phase plot",
              file=sys.stderr)
        return

    phases_tsv = os.path.join(args.out_dir, f"{args.prefix}consolidated_phases.tsv")
    phase_names = [g for g, _ in PHASE_GROUPS]
    with open(phases_tsv, "w") as fh:
        fh.write("event\t" + "\t".join(f"{p}_secs" for p in phase_names)
                 + "\ttotal_secs\n")
        for event in sorted(phase_table.keys()):
            durs = phase_table[event]
            row = [event] + [str(durs[p]) for p in phase_names] \
                + [str(sum(durs.values()))]
            fh.write("\t".join(row) + "\n")
    print(f"Wrote {phases_tsv}", file=sys.stderr)

    # Stacked bar plot of phase durations per genome, sorted by total time.
    events_sorted = sorted(phase_table.keys(),
                           key=lambda e: -sum(phase_table[e].values()))
    n_total = len(events_sorted)
    if args.top and args.top > 0:
        events_plot = events_sorted[:args.top]
    else:
        events_plot = events_sorted
    n = len(events_plot)
    phase_colors = {"caf": "#4c72b0", "bar": "#dd8452",
                    "reference": "#55a467", "hal": "#c44e52"}

    fig, ax = plt.subplots(figsize=(max(8, n * 0.22), 6))
    bottoms = np.zeros(n)
    for phase in phase_names:
        vals = np.array([fmt_hours(phase_table[e][phase]) for e in events_plot])
        ax.bar(np.arange(n), vals, bottom=bottoms,
               color=phase_colors[phase], label=phase, width=0.9)
        bottoms += vals
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([display_name(e, name_map) for e in events_plot],
                       rotation=90, fontsize=8)
    ax.set_ylabel("hours")
    title = f"cactus_consolidated phase breakdown (top {n} of {n_total} ancestors by walltime)" \
            if n < n_total else \
            f"cactus_consolidated phase breakdown ({n} ancestors)"
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.tight_layout()
    phases_png = os.path.join(args.out_dir, f"{args.prefix}consolidated_phases.png")
    fig.savefig(phases_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {phases_png}", file=sys.stderr)

    # Distribution of total runtime + per-phase violin for whole-tree view
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    totals_h = np.array([sum(phase_table[e].values()) / 3600.0
                         for e in events_sorted])
    axes[0].hist(totals_h, bins=40, color="#4c72b0", edgecolor="black")
    axes[0].set_xlabel("total walltime (hours)")
    axes[0].set_ylabel("# ancestors")
    axes[0].set_title(f"Distribution of cactus_consolidated walltime "
                      f"(n={n_total})")

    raw = [np.array([phase_table[e][p] / 3600.0 for e in events_sorted])
           for p in phase_names]
    # Compute violins in log-space so KDE shapes work across decades.
    plot_data = [np.log10(np.clip(arr, 1e-3, None)) for arr in raw]
    positions = np.arange(1, len(phase_names) + 1)
    vp = axes[1].violinplot(plot_data, positions=positions, showmeans=False,
                            showmedians=True, showextrema=True, widths=0.8)
    for body, p in zip(vp["bodies"], phase_names):
        body.set_facecolor(phase_colors[p])
        body.set_edgecolor("black")
        body.set_alpha(0.85)
    for part in ("cbars", "cmins", "cmaxes", "cmedians"):
        if part in vp:
            vp[part].set_color("black")
            vp[part].set_linewidth(0.8)
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels(phase_names)
    all_log = np.concatenate(plot_data)
    lo = np.floor(all_log.min())
    hi = np.ceil(all_log.max())
    ticks = np.arange(lo, hi + 1)
    axes[1].set_yticks(ticks)
    axes[1].set_yticklabels(
        [f"{10 ** t:g}" if t >= 0 else f"{10 ** t:.3g}" for t in ticks]
    )
    axes[1].set_ylabel("hours (log)")
    axes[1].set_title("Per-phase walltime across all ancestors")
    fig.tight_layout()
    dist_png = os.path.join(args.out_dir,
                            f"{args.prefix}consolidated_phases_dist.png")
    fig.savefig(dist_png, dpi=150)
    plt.close(fig)
    print(f"Wrote {dist_png}", file=sys.stderr)

    # ---- cactus_consolidated peak memory per ancestor ----
    if consolidated_meta:
        plot_consolidated_memory(consolidated_meta, args.out_dir,
                                 args.prefix, args.top, name_map)


def plot_lastz_per_genome(per_genome, out_dir, prefix, top, name_map):
    """Bar chart of per-genome cumulative lastz walltime (top N) plus a
    distribution histogram across all genomes. A "genome" here is anything
    that appeared as either side of a lastz pair, so each lastz job
    contributes its runtime to two genomes.
    """
    events_sorted = sorted(per_genome.keys(),
                           key=lambda e: -per_genome[e])
    n_total = len(events_sorted)
    if top and top > 0:
        events_plot = events_sorted[:top]
    else:
        events_plot = events_sorted
    n = len(events_plot)
    hours = [per_genome[e] / 3600.0 for e in events_plot]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, n * 0.22), 9),
                             gridspec_kw={"height_ratios": [3, 2]})

    ax = axes[0]
    ax.bar(np.arange(n), hours, color="#1f77b4", width=0.9)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([display_name(e, name_map) for e in events_plot],
                       rotation=90, fontsize=8)
    ax.set_ylabel("total lastz walltime (hours)")
    title = (f"Per-genome cumulative lastz walltime "
             f"(top {n} of {n_total} genomes)") \
        if n < n_total else \
        f"Per-genome cumulative lastz walltime ({n} genomes)"
    ax.set_title(title)

    ax2 = axes[1]
    all_h = [per_genome[e] / 3600.0 for e in events_sorted]
    ax2.hist(all_h, bins=40, color="#1f77b4", edgecolor="black")
    ax2.set_xlabel("total lastz walltime (hours)")
    ax2.set_ylabel("# genomes")
    ax2.set_title(f"Distribution of per-genome lastz walltime "
                  f"(n={n_total})")
    ax2.axvline(np.median(all_h), color="red", linestyle="--",
                label=f"median = {np.median(all_h):.1f} h")
    ax2.axvline(np.percentile(all_h, 95), color="orange", linestyle="--",
                label=f"p95 = {np.percentile(all_h, 95):.1f} h")
    ax2.legend()

    fig.tight_layout()
    out = os.path.join(out_dir, f"{prefix}lastz_by_genome.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}", file=sys.stderr)

    tsv = os.path.join(out_dir, f"{prefix}lastz_by_genome.tsv")
    with open(tsv, "w") as fh:
        fh.write("genome\tdisplay_name\ttotal_lastz_h\n")
        for e in events_sorted:
            fh.write(f"{e}\t{display_name(e, name_map)}\t"
                     f"{per_genome[e] / 3600.0:.3f}\n")
    print(f"Wrote {tsv}", file=sys.stderr)


def plot_consolidated_memory(meta, out_dir, prefix, top, name_map):
    """Bar chart of cactus_consolidated peak memory for the top-N ancestors,
    plus a histogram of all ancestors' memory usage.
    """
    events_sorted = sorted(meta.keys(),
                           key=lambda e: -meta[e]["mem_bytes"])
    n_total = len(events_sorted)
    if top and top > 0:
        events_plot = events_sorted[:top]
    else:
        events_plot = events_sorted
    n = len(events_plot)
    mems_gib = [meta[e]["mem_bytes"] / (2 ** 30) for e in events_plot]
    secs = [meta[e]["secs"] / 3600.0 for e in events_plot]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, n * 0.22), 9),
                             gridspec_kw={"height_ratios": [3, 2]})

    # Top: per-event peak memory bar; color encodes walltime so the user
    # can see whether memory and time scale together.
    ax = axes[0]
    norm = plt.Normalize(vmin=0, vmax=max(secs) if secs else 1)
    cmap = plt.get_cmap("viridis")
    bar_colors = [cmap(norm(s)) for s in secs]
    ax.bar(np.arange(n), mems_gib, color=bar_colors, width=0.9)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([display_name(e, name_map) for e in events_plot],
                       rotation=90, fontsize=8)
    ax.set_ylabel("peak memory (GiB)")
    title = (f"cactus_consolidated peak memory "
             f"(top {n} of {n_total} ancestors)") \
        if n < n_total else \
        f"cactus_consolidated peak memory ({n} ancestors)"
    ax.set_title(title)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label("walltime (h)")

    # Bottom: histogram of all ancestors' memory.
    ax2 = axes[1]
    all_mem = [meta[e]["mem_bytes"] / (2 ** 30) for e in events_sorted]
    ax2.hist(all_mem, bins=40, color="#4c72b0", edgecolor="black")
    ax2.set_xlabel("peak memory (GiB)")
    ax2.set_ylabel("# ancestors")
    ax2.set_title(f"Distribution of cactus_consolidated peak memory "
                  f"(n={n_total})")
    ax2.axvline(np.median(all_mem), color="red", linestyle="--",
                label=f"median = {np.median(all_mem):.1f} GiB")
    ax2.axvline(np.percentile(all_mem, 95), color="orange", linestyle="--",
                label=f"p95 = {np.percentile(all_mem, 95):.1f} GiB")
    ax2.legend()

    fig.tight_layout()
    out = os.path.join(out_dir, f"{prefix}consolidated_memory.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}", file=sys.stderr)

    tsv = os.path.join(out_dir, f"{prefix}consolidated_memory.tsv")
    with open(tsv, "w") as fh:
        fh.write("event\tpeak_mem_GiB\twalltime_h\n")
        for e in events_sorted:
            fh.write(f"{e}\t{meta[e]['mem_bytes'] / (2 ** 30):.2f}\t"
                     f"{meta[e]['secs'] / 3600.0:.3f}\n")
    print(f"Wrote {tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
