#!/usr/bin/env bash
set -euo pipefail

GENOME="${GENOME:-hg38}"
OUTDIR="${OUTDIR:-data/loci}"

mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# Helper: fetch locus + convert repeats to local coords
# ------------------------------------------------------------
fetch_and_localize () {
  local chrom="$1"
  local start="$2"
  local end="$3"
  local name="$4"

  echo "> Fetching ${GENOME} ${chrom}:${start}-${end} -> ${name}"
  python scripts/fetch_ucsc_locus.py \
    --genome "$GENOME" \
    --chrom "$chrom" --start "$start" --end "$end" \
    --name "$name" \
    --outdir "$OUTDIR"

  echo "> Converting repeats BED to local coords -> ${name}.repeats.local.bed"
  python scripts/bed_to_local.py \
    --in-bed  "${OUTDIR}/${name}.repeats.bed" \
    --out-bed "${OUTDIR}/${name}.repeats.local.bed" \
    --locus-start "$start" --locus-end "$end" \
    --chrom-filter "$chrom" \
    --contig "$name"
}

# ------------------------------------------------------------
# Locus 1: HBB +/- 10kb
# ------------------------------------------------------------
fetch_and_localize \
  chr11 5215464 5237071 \
  locus1_hg38_chr11_5215464_5237071_HBB

# ------------------------------------------------------------
# Locus 2: repeat-rich window
# ------------------------------------------------------------

SEED="${SEED:-5}"
echo "> Finding repeat-rich window (seed=${SEED})"
found="$(python scripts/find_repeat_window_ucsc.py --chrom chr8 --seed "$SEED" | tail -n 1)"
echo "$found"
# Expected format: FOUND   hg38    chr8:96243643-96263643  Alu=7   LTR=3
# gene overlapping window is MTERF3
GENE="${GENE:-MTERF3}"

coord="$(echo "$found" | awk '{print $3}')"
chrom="$(echo "$coord" | cut -d: -f1)"
range="$(echo "$coord" | cut -d: -f2)"
start="$(echo "$range" | cut -d- -f1)"
end="$(echo "$range" | cut -d- -f2)"
name="locus2_${GENOME}_${chrom}_${start}_${end}_${GENE}"
fetch_and_localize "$chrom" "$start" "$end" "$name"

echo "Done. Contents:"
ls -lh "$OUTDIR"
