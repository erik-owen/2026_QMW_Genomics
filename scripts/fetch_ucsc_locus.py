#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import requests

UCSC = "https://api.genome.ucsc.edu"

def _get_int(x: dict, *keys: str) -> int:
    for k in keys:
        if k in x and x[k] is not None:
            return int(x[k])
    raise KeyError(f"None of the keys {keys} found. Available keys: {sorted(x.keys())[:30]} ...")

def fetch_sequence(genome: str, chrom: str, start: int, end: int) -> str:
    r = requests.get(
        f"{UCSC}/getData/sequence",
        params={"genome": genome, "chrom": chrom, "start": start, "end": end},
        timeout=60,
    )
    r.raise_for_status()
    js = r.json()
    seq = js.get("dna", "")
    if not seq:
        raise RuntimeError(f"No sequence returned for {genome} {chrom}:{start}-{end}")
    return seq.upper()

def fetch_rmsk(genome: str, chrom: str, start: int, end: int) -> list[dict]:
    r = requests.get(
        f"{UCSC}/getData/track",
        params={"genome": genome, "track": "rmsk", "chrom": chrom, "start": start, "end": end},
        timeout=60,
    )
    r.raise_for_status()
    return r.json().get("rmsk", [])

def write_fasta(path: Path, name: str, seq: str, width: int = 60) -> None:
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i:i+width] + "\n")

def write_repeats_bed(path: Path, default_chrom: str, rows: list[dict]) -> None:
    """
    BED6 + extra cols:
      chrom start end name score strand repClass repFamily
    """
    with open(path, "w") as out:
        for x in rows:
            chrom = x.get("genoName") or x.get("chrom") or default_chrom
            start = _get_int(x, "genoStart", "chromStart", "start")
            end   = _get_int(x, "genoEnd",   "chromEnd",   "end")

            strand = x.get("strand", ".")
            rep = x.get("repName", "repeat")
            repClass = x.get("repClass", "")
            repFamily = x.get("repFamily", "")

            out.write(f"{chrom}\t{start}\t{end}\t{rep}\t0\t{strand}\t{repClass}\t{repFamily}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--name", required=True, help="Contig name to use in FASTA header")
    ap.add_argument("--outdir", default="data/loci")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seq = fetch_sequence(args.genome, args.chrom, args.start, args.end)
    rows = fetch_rmsk(args.genome, args.chrom, args.start, args.end)

    fa = outdir / f"{args.name}.fa"
    bed = outdir / f"{args.name}.repeats.bed"
    write_fasta(fa, args.name, seq)
    write_repeats_bed(bed, args.chrom, rows)

    # quick summary
    n_alu = sum(1 for r in rows if r.get("repClass") == "SINE" and ("Alu" in (r.get("repFamily","") + r.get("repName",""))))
    n_ltr = sum(1 for r in rows if r.get("repClass") == "LTR")
    print(f"Wrote {fa}")
    print(f"Wrote {bed}")
    print(f"Repeat summary: Alu={n_alu} LTR={n_ltr} total={len(rows)}")

if __name__ == "__main__":
    main()

# Example Usage:

# ❯ python scripts/fetch_ucsc_locus.py \
#   --genome hg38 \
#   --chrom chr11 --start 5215464 --end 5237071 \
#   --name locus1_hg38_chr11_5215464_5237071_HBB
# Wrote data/loci/locus1_hg38_chr11_5215464_5237071_HBB.fa
# Wrote data/loci/locus1_hg38_chr11_5215464_5237071_HBB.repeats.bed
# Repeat summary: Alu=4 LTR=0 total=28

# ❯ python scripts/fetch_ucsc_locus.py \
#   --chrom chr8 --start 96243643 --end 96263643 \
#   --name locus2_hg38_chr8_96243643_96263643_mterf3       
# Wrote data/loci/locus2_hg38_chr8_96243643_96263643_mterf3.fa
# Wrote data/loci/locus2_hg38_chr8_96243643_96263643_mterf3.repeats.bed
# Repeat summary: Alu=7 LTR=3 total=34