#!/usr/bin/env bash
set -euo pipefail

echo "Installing tools..."
sudo apt-get -qq update
sudo apt-get -qq install -y bowtie2 samtools bcftools tabix pigz

# echo "Downloading workshop data..."
# mkdir -p data/loci
# BASE="https://raw.githubusercontent.com/erik-owen/2026_QMW_Genomics/main/data/loci"
# for f in locus1.fa locus1.truth.vcf locus1.reads_R1.fastq.gz; do
#   curl -L -o "data/loci/$f" "$BASE/$f"
# done

echo "Tool versions:"
bowtie2 --version | head -n 1 || true
samtools --version | head -n 1 || true
bcftools --version | head -n 1 || true