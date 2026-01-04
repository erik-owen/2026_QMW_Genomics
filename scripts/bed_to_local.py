#!/usr/bin/env python3
"""
Convert a genomic-coordinate BED (e.g., UCSC rmsk) into local coordinates for a
mini-contig locus FASTA.

Example:
  Input BED:  chr8  96245000  96245500  AluY  0  +
  locus_start=96243643 locus_end=96263643 contig=locus2_...
  Output BED: locus2_...  (96245000-96243643)=1357  ... etc

Rules:
- shift: local = genomic - locus_start
- clip intervals to [0, locus_len]
- drop intervals that don't overlap locus
- rewrite chrom to --contig
- preserve extra columns (name/score/strand/repClass/repFamily...) as-is
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import List, Optional


def parse_int(x: str) -> int:
    try:
        return int(x)
    except Exception as e:
        raise ValueError(f"Expected int, got {x!r}") from e


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-bed", required=True, help="Input BED (genomic coords).")
    ap.add_argument("--out-bed", required=True, help="Output BED (local coords).")
    ap.add_argument("--locus-start", type=int, required=True, help="Genomic start (0-based) of fetched locus window.")
    ap.add_argument("--locus-end", type=int, required=True, help="Genomic end (0-based, exclusive) of fetched locus window.")
    ap.add_argument("--contig", required=True, help="Chrom/contig name to write in output BED.")
    ap.add_argument("--chrom-filter", default=None, help="Only convert intervals with this chrom (e.g., chr8).")
    ap.add_argument("--keep-header", action="store_true", help="Copy through 'track'/'browser' lines as comments.")
    args = ap.parse_args()

    in_path = Path(args.in_bed)
    out_path = Path(args.out_bed)

    locus_len = args.locus_end - args.locus_start
    if locus_len <= 0:
        raise ValueError("locus_end must be > locus_start")

    n_in = 0
    n_written = 0
    n_dropped = 0

    with open(in_path) as inp, open(out_path, "w") as out:
        for line in inp:
            if not line.strip():
                continue
            if line.startswith(("track", "browser", "#")):
                if args.keep_header:
                    out.write("# " + line if not line.startswith("#") else line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            chrom = parts[0]
            if args.chrom_filter and chrom != args.chrom_filter:
                continue

            start = parse_int(parts[1])
            end = parse_int(parts[2])
            if end <= start:
                continue

            n_in += 1

            # shift to local
            ls = start - args.locus_start
            le = end - args.locus_start

            # drop if no overlap with [0, locus_len)
            if le <= 0 or ls >= locus_len:
                n_dropped += 1
                continue

            # clip to locus
            if ls < 0:
                ls = 0
            if le > locus_len:
                le = locus_len

            out_parts = [args.contig, str(ls), str(le)] + parts[3:]
            out.write("\t".join(out_parts) + "\n")
            n_written += 1

    print(f"Input intervals considered: {n_in}")
    print(f"Written: {n_written}")
    print(f"Dropped (outside locus): {n_dropped}")
    print(f"Output: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
