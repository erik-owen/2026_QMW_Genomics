#!/usr/bin/env bash
set -euo pipefail

fail(){ echo "ERROR: $*" >&2; exit 1; }

for fa in data/loci/*.fa; do
  contig="$(grep -m1 '^>' "$fa" | sed 's/^>//' | awk '{print $1}')"
  # check any truth vcfs that correspond to this contig
  for vcf in data/loci/*.truth.vcf; do
    grep -q "^##contig=<ID=${contig}," "$vcf" || continue
    # ensure CHROM column only contains this contig
    bad=$(awk -v c="$contig" 'BEGIN{FS="\t"} $0 !~ /^#/ && $1!=c {print $1; exit 0}' "$vcf" || true)
    [[ -z "$bad" ]] || fail "$vcf has CHROM=$bad but expected $contig"
  done
done

echo "Toy bundle checks passed."
