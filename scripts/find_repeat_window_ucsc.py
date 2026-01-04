#!/usr/bin/env python3
import argparse, random, requests

UCSC = "https://api.genome.ucsc.edu"

def get_rmsk(genome, chrom, start, end):
    url = f"{UCSC}/getData/track"
    params = {"genome": genome, "track": "rmsk", "chrom": chrom, "start": start, "end": end}
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json().get("rmsk", [])

def has_alu_and_ltr(rmsk_rows, min_alu=2):
    alu = [x for x in rmsk_rows if x.get("repClass") == "SINE" and "Alu" in (x.get("repFamily") or x.get("repName") or "")]
    ltr = [x for x in rmsk_rows if x.get("repClass") == "LTR"]
    return (len(alu) >= min_alu) and (len(ltr) >= 1), len(alu), len(ltr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--chrom", default="chr8")
    ap.add_argument("--window", type=int, default=20000)
    ap.add_argument("--tries", type=int, default=200)
    ap.add_argument("--chrom_len", type=int, default=145138636)  # chr8 size: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
    ap.add_argument("--min_alu", type=int, default=2)
    ap.add_argument("--seed", type=int, default=5)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    for _ in range(args.tries):
        start = rng.randint(0, args.chrom_len - args.window - 1)
        end = start + args.window
        rows = get_rmsk(args.genome, args.chrom, start, end)
        ok, nalu, nltr = has_alu_and_ltr(rows, min_alu=args.min_alu)
        if ok:
            print(f"FOUND\t{args.genome}\t{args.chrom}:{start}-{end}\tAlu={nalu}\tLTR={nltr}")
            return
    raise SystemExit("No suitable window found; increase --tries or change --chrom/--window.")

if __name__ == "__main__":
    main()
