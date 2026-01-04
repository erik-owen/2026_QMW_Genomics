#!/usr/bin/env bash
set -euo pipefail

echo "> Installing tools..."
sudo apt-get -qq update
sudo DEBIAN_FRONTEND=noninteractive apt-get -qq install -y \
  bowtie2 samtools bcftools tabix pigz \
  tree
echo "Tools installed."
echo
echo "Tool versions:"
bt2_ver="$(bowtie2 --version 2>&1 || true)"
echo "${bt2_ver%%$'\n'*}"          # first line, no SIGPIPE
samtools --version | head -n 1 || true
bcftools --version | head -n 1 || true

echo "> Downloading workshop data..."
mkdir -p data/loci

BASE="https://raw.githubusercontent.com/erik-owen/2026_QMW_Genomics/main/data/loci"
FILES=(
  locus1_hg38_chr11_5215464_5237071_HBB.fa
  locus1_hg38_chr11_5215464_5237071_HBB.repeats.local.bed
  locus1_HBB.ind1_noHbS.reads_R1.fastq.gz
  locus1_HBB.ind2_HbS_carrier.reads_R1.fastq.gz
  locus2_hg38_chr8_96243643_96263643_MTERF3.fa
  locus2_hg38_chr8_96243643_96263643_MTERF3.repeats.local.bed
  locus2_MTERF3_repeat.ind1_repeat_mix.reads_R1.fastq.gz
)

for f in "${FILES[@]}"; do
  url="$BASE/$f"
  out="data/loci/$f"
  echo "  - $f"
  curl -L --fail --retry 3 --retry-delay 1 -o "$out" "$url"
done

echo
echo "Downloaded files:"
ls -lh data/loci | sed 's/^/  /'