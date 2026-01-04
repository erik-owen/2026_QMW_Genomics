#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# build_toy_data.sh
#
# Generate deterministic toy FASTQ + truth VCF for workshop loci.
#
# Usage:
#   ./scripts/build_toy_data.sh
#   ./scripts/build_toy_data.sh --force
# ------------------------------------------------------------

FORCE=0
if [[ "${1:-}" == "--force" ]]; then
  FORCE=1
fi

OUTDIR="data/loci"
mkdir -p "$OUTDIR"

die() { echo "ERROR: $*" >&2; exit 1; }

need_file() { [[ -f "$1" ]] || die "Missing file: $1"; }

maybe_rm() {
  local f="$1"
  if [[ -e "$f" ]]; then
    if [[ "$FORCE" -eq 1 ]]; then
      rm -f "$f"
    else
      die "Output exists: $f (re-run with --force to overwrite)"
    fi
  fi
}

# ------------------------------------------------------------
# Inputs (from build_loci.sh)
# ------------------------------------------------------------
HBB_FASTA="$OUTDIR/locus1_hg38_chr11_5215464_5237071_HBB.fa"
MTERF_FASTA="$OUTDIR/locus2_hg38_chr8_96243643_96263643_MTERF3.fa"
MTERF_REP_LOCAL="$OUTDIR/locus2_hg38_chr8_96243643_96263643_MTERF3.repeats.local.bed"

need_file "$HBB_FASTA"
need_file "$MTERF_FASTA"
need_file "$MTERF_REP_LOCAL"
need_file "scripts/make_toy_data.py"

# ------------------------------------------------------------
# Global sim settings 
#
# Feel free to adjust these if you'd like to make your own
# toy data! 
# ------------------------------------------------------------

# Illumina-like values
READLEN=100
Q_START=35
Q_END=18
Q_SD=3.0

# Thousands are reasonable for colab
NREADS_HBB_1=2471
NREADS_HBB_2=2732
NREADS_REPEAT=2943

# ------------------------------------------------------------
# Helper to run make_toy_data.py and protect outputs
# ------------------------------------------------------------
run_make() {
  local fasta="$1"
  local prefix="$2"
  local sample="$3"
  shift 3

  local tag="${prefix}.${sample}"
  local vcf="$OUTDIR/${tag}.truth.vcf"
  local r1="$OUTDIR/${tag}.reads_R1.fastq.gz"
  local r2="$OUTDIR/${tag}.reads_R2.fastq.gz"

  maybe_rm "$vcf"
  maybe_rm "$r1"
  maybe_rm "$r2"  # Only matters for paired end reads

  echo "> Generating: ${tag}"
  python scripts/make_toy_data.py \
    --fasta "$fasta" \
    --outdir "$OUTDIR" \
    --prefix "$prefix" \
    --sample "$sample" \
    --readlen "$READLEN" \
    --q-start "$Q_START" \
    --q-end "$Q_END" \
    --q-sd "$Q_SD" \
    "$@"

  [[ -f "$vcf" ]] || die "Expected VCF not created: $vcf"
  [[ -f "$r1"  ]] || die "Expected R1 not created:  $r1"
}

# ------------------------------------------------------------
# 1) HBB locus
#   - Individual 1: wild-type                 -> alt_fraction 0.0
#   - Individual 2: HbS carrier (heterozygous) -> alt_fraction 0.5 + add-hbb-hbs
# ------------------------------------------------------------
HBB_PREFIX="locus1_HBB"

run_make "$HBB_FASTA" "$HBB_PREFIX" "ind1_noHbS" \
  --seed 1 \
  --alt-fraction 0.0 \
  --n_reads "$NREADS_HBB_1"

run_make "$HBB_FASTA" "$HBB_PREFIX" "ind2_HbS_carrier" \
  --seed 2 \
  --add-hbb-hbs \
  --alt-fraction 0.5 \
  --n_reads "$NREADS_HBB_2"

# ------------------------------------------------------------
# 2) Repeat-rich locus (MTERF3 window)
#   - Put 2 SNPs in repeats + 2 outside repeats to demonstrate MAPQ/multimapping effects
# ------------------------------------------------------------
REP_PREFIX="locus2_MTERF3_repeat"

run_make "$MTERF_FASTA" "$REP_PREFIX" "ind1_repeat_mix" \
  --seed 3 \
  --repeats-local-bed "$MTERF_REP_LOCAL" \
  --n_snps 4 \
  --n_snps_in_repeats 2 \
  --alt-fraction 1.0 \
  --n_reads "$NREADS_REPEAT"

echo
echo "Done. Outputs in $OUTDIR:"
ls -lh "$OUTDIR" | sed 's/^/  /'
