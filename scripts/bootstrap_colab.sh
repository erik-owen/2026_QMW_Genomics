#!/usr/bin/env bash
set -euo pipefail

echo "> Installing tools..."
sudo apt-get -qq update
sudo DEBIAN_FRONTEND=noninteractive apt-get -qq install -y \
  bowtie2 samtools bcftools tabix pigz fastqc perl \
  tree

# cutadapt for TrimGalore
python3 -m pip -q install --upgrade pip
python3 -m pip -q install cutadapt

# Download TrimGalore from GitHub
TG_DIR="/opt/trim_galore"
mkdir -p "${TG_DIR}"
curl -L -o /tmp/trim_galore.tar.gz \
  https://github.com/FelixKrueger/TrimGalore/archive/refs/heads/master.tar.gz
tar -xzf /tmp/trim_galore.tar.gz -C "${TG_DIR}" --strip-components 1

chmod +x "${TG_DIR}/trim_galore"
ln -sf "${TG_DIR}/trim_galore" /usr/local/bin/trim_galore

echo "Tools installed."
echo
echo "Tool versions:"
bt2_ver="$(bowtie2 --version 2>&1 || true)"
echo "${bt2_ver%%$'\n'*}"          # first line, no SIGPIPE
samtools --version | head -n 1 || true
bcftools --version | head -n 1 || true
trim_galore --version | head -n 1 || true
cutadapt --version | head -n 1 || true
fastqc --version | head -n 1 || true


echo
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