#!/bin/bash
# rename-gff-refseq-to-genbank.sh
#
# Download a RefSeq GFF for a GCF assembly and rewrite its chrom column
# from RefSeq names (NC_/NW_/NT_/...) to the GenBank names (CM_/scaffold-
# id) used by the paired GCA assembly. Useful when an alignment contains
# the GCA version of an assembly with the same physical contigs as a
# published GCF, and you want to feed the GCF's RefSeq annotation to
# cactus-phast (or any downstream tool) keyed on the GCA contig names.
#
# Usage:
#   rename-gff-refseq-to-genbank.sh <GCF_accession> <GCA_accession> [output_gff_path]
#
# Default output path: <GCA_accession>.refseq-renamed.gff.gz
#
# Example:
#   rename-gff-refseq-to-genbank.sh GCF_048771995.1 GCA_048771995.1

set -euo pipefail

usage() {
    sed -nE 's/^# ?//; 4,/^$/p' "$0" >&2
    exit 1
}

[ $# -lt 2 ] && usage

gcf="$1"
gca="$2"
out_gff="${3:-${gca}.refseq-renamed.gff.gz}"

[[ "${gcf}" == GCF_* ]] || { echo "ERROR: first arg must start with GCF_ (got ${gcf})" >&2; exit 1; }
[[ "${gca}" == GCA_* ]] || { echo "ERROR: second arg must start with GCA_ (got ${gca})" >&2; exit 1; }

# Resolve the full FTP subdirectory for an NCBI assembly accession. NCBI
# splits the path as /genomes/all/<GCx>/<3>/<3>/<3>/<dir>, where <dir>
# begins with the accession; the suffix after the accession is the
# assembly name and varies per release.
ncbi_dir() {
    local acc="$1"
    local prefix="${acc%_*}"        # GCF or GCA
    local digits="${acc#*_}"
    digits="${digits%%.*}"
    local parent="https://ftp.ncbi.nlm.nih.gov/genomes/all/${prefix}/${digits:0:3}/${digits:3:3}/${digits:6:3}/"
    local subdir
    subdir=$(curl -fsS "${parent}" \
        | sed -nE "s|.*href=\"(${acc}_[^/\"]+)/?\".*|\1|p" \
        | head -1)
    if [ -z "${subdir}" ]; then
        echo "ERROR: could not resolve NCBI subdir for ${acc} under ${parent}" >&2
        return 1
    fi
    echo "${parent}${subdir}"
}

work=$(mktemp -d)
trap "rm -rf '${work}'" EXIT

echo "Resolving FTP paths..." >&2
gcf_dir=$(ncbi_dir "${gcf}")
gca_dir=$(ncbi_dir "${gca}")
echo "  GCF: ${gcf_dir}" >&2
echo "  GCA: ${gca_dir}" >&2

# Either side's assembly_report.txt has both naming conventions; use the
# GCA one (closer to the target naming).
report="${work}/assembly_report.txt"
echo "Fetching assembly report..." >&2
curl -fsS "${gca_dir}/$(basename "${gca_dir}")_assembly_report.txt" -o "${report}"

# assembly_report.txt columns (after #-prefixed header lines):
#   1=Sequence-Name 2=Sequence-Role 3=Assigned-Molecule 4=...Type
#   5=GenBank-Accn 6=Relationship 7=RefSeq-Accn 8=Assembly-Unit
#   9=Sequence-Length 10=UCSC-style-name
mapping="${work}/rs2gb.tsv"
awk -F'\t' '!/^#/ && NF >= 7 && $7 != "na" && $7 != "" {print $7"\t"$5}' "${report}" > "${mapping}"
n_pairs=$(wc -l < "${mapping}")
echo "  ${n_pairs} RefSeq->GenBank pairs in assembly report" >&2
[ "${n_pairs}" -gt 0 ] || { echo "ERROR: no RefSeq->GenBank pairs in report" >&2; exit 1; }

# Pull the GFF.
gff_in="${work}/in.gff.gz"
echo "Fetching GFF..." >&2
curl -fsS "${gcf_dir}/$(basename "${gcf_dir}")_genomic.gff.gz" -o "${gff_in}"

# Rewrite chrom column in a single pass via awk's hash-table lookup.
echo "Rewriting chrom column..." >&2
zcat "${gff_in}" \
    | awk -F'\t' -v OFS='\t' -v MAP="${mapping}" '
        BEGIN {
            while ((getline line < MAP) > 0) {
                split(line, a, "\t")
                m[a[1]] = a[2]
            }
            close(MAP)
        }
        /^#/ {print; next}
        NF > 1 && ($1 in m) { $1 = m[$1] }
        {print}
    ' \
    | gzip > "${out_gff}"

# Sanity check: anything still bearing a RefSeq accession in col 1?
unmapped=$(zcat "${out_gff}" | awk -F'\t' '!/^#/ && $1 ~ /^(NC|NW|NT|NZ)_/' | wc -l)
total_features=$(zcat "${out_gff}" | awk -F'\t' '!/^#/' | wc -l)
distinct_chroms=$(zcat "${out_gff}" | awk -F'\t' '!/^#/ {print $1}' | sort -u | wc -l)

echo "Done." >&2
echo "  output:                 ${out_gff}" >&2
echo "  total feature lines:    ${total_features}" >&2
echo "  distinct chroms in out: ${distinct_chroms}" >&2
echo "  unmapped (RefSeq left): ${unmapped}" >&2
if [ "${unmapped}" -gt 0 ]; then
    echo "" >&2
    echo "  Unmapped lines correspond to RefSeq contigs with no GenBank" >&2
    echo "  pair in the assembly report (typically unplaced / patch / alt" >&2
    echo "  contigs that aren't in the GCA submission). cactus-phast will" >&2
    echo "  simply produce no 4d sites for those contigs; the rest of the" >&2
    echo "  annotation is fine." >&2
fi
