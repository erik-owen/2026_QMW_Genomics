#!/usr/bin/env python3
"""
scripts/make_toy_data.py

Given a FASTA for a locus, generate Illumina-like FASTQ reads and a ground-truth VCF.

Features:
- Deterministic RNG via --seed
- Supports SNPs + small indels as VCF-style REF/ALT replacements
- Optional built-in pathogenic HBB HbS variant (rs334; HBB c.20A>T) when locus spans GRCh38 chr11:5227002
  (HbS is classically described as a codon change GAG -> GTG in HBB)
- Illumina-ish per-cycle quality decay + sampling from a distribution
- Substitution errors derived from Q via Phred definition p = 10^(-Q/10)
- Single-end or paired-end
- Optional mixture of reference + ALT haplotype reads via --alt-fraction (useful for heterozygous examples)

Outputs:
  <outdir>/<prefix>.<sample>.truth.vcf
  <outdir>/<prefix>.<sample>.reads_R1.fastq.gz
  <outdir>/<prefix>.<sample>.reads_R2.fastq.gz   (if --paired)

Notes:
- Reads are simulated from an "ALT haplotype" = reference locus with inserted truth variants.
- Alignment is performed later against the original reference locus FASTA (mini-contig).
"""

from __future__ import annotations
import argparse
import gzip
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

DNA = "ACGT"


# ----------------------------
# Basic sequence helpers
# ----------------------------

def rc(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def read_single_fasta(path: Path) -> tuple[str, str]:
    name = None
    parts: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].split()[0]
                else:
                    raise ValueError(f"{path} has multiple FASTA records; expected one.")
            else:
                parts.append(line.upper())
    if name is None:
        raise ValueError(f"{path} missing FASTA header.")
    seq = "".join(parts)
    if not seq:
        raise ValueError(f"{path} empty sequence.")
    return name, seq

def phred_char(q: int) -> str:
    # Phred+33
    return chr(q + 33)

def q_to_p_err(q: int) -> float:
    # Phred definition: Q = -10 log10(p)
    return 10 ** (-q / 10)

def clamp_int(x: int, lo: int, hi: int) -> int:
    return lo if x < lo else hi if x > hi else x

def mutate_base(refb: str, rng: random.Random) -> str:
    alts = [b for b in DNA if b != refb]
    return rng.choice(alts)


# ----------------------------
# BED intervals (local coords)
# ----------------------------

@dataclass(frozen=True)
class Interval:
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive

def parse_local_bed_intervals(path: Path, contig_filter: Optional[str] = None) -> List[Interval]:
    intervals: List[Interval] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if contig_filter and chrom != contig_filter:
                continue
            s = int(parts[1]); e = int(parts[2])
            if e > s:
                intervals.append(Interval(s, e))
    intervals.sort(key=lambda x: (x.start, x.end))
    return intervals

def in_any_interval(pos0: int, intervals: List[Interval]) -> bool:
    for itv in intervals:
        if itv.start <= pos0 < itv.end:
            return True
    return False


# ----------------------------
# Variants (VCF-style)
# ----------------------------

@dataclass
class Variant:
    pos1: int      # 1-based position on the mini-contig
    ref: str
    alt: str
    vid: str = "."
    info: str = "."

def apply_variant(seq: str, v: Variant) -> str:
    """
    Apply a VCF-style REF/ALT replacement at v.pos1.
    - pos1 points to the first base of REF on the reference.
    - REF must match reference sequence exactly.
    """
    pos0 = v.pos1 - 1
    if pos0 < 0 or pos0 >= len(seq):
        raise ValueError(f"Variant position out of range: {v.pos1}")
    ref = v.ref.upper()
    alt = v.alt.upper()
    if seq[pos0:pos0+len(ref)] != ref:
        ctx = seq[max(0, pos0-10):min(len(seq), pos0+len(ref)+10)]
        raise ValueError(
            f"REF mismatch at pos1={v.pos1}. Expected {ref}, saw {seq[pos0:pos0+len(ref)]}. "
            f"Context: {ctx}"
        )
    return seq[:pos0] + alt + seq[pos0+len(ref):]

def write_vcf(path: Path, contig: str, contig_len: int, vars: List[Variant]) -> None:
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={contig},length={contig_len}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in vars:
            out.write(f"{contig}\t{v.pos1}\t{v.vid}\t{v.ref}\t{v.alt}\t.\tPASS\t{v.info}\n")


# ----------------------------
# Infer genomic window from contig name
# e.g. locus1_hg38_chr11_5215464_5237071_HBB
# ----------------------------

def infer_genomic_window_from_name(contig: str) -> tuple[Optional[str], Optional[int], Optional[int]]:
    m = re.search(r"_(chr[0-9XYM]+)_(\d+)_(\d+)_", contig)
    if not m:
        return None, None, None
    chrom = m.group(1)
    start = int(m.group(2))  # 0-based locus start used in your file naming
    end = int(m.group(3))
    return chrom, start, end


# ----------------------------
# Built-in HBB HbS (rs334)
# GRCh38 chr11:5227002 ; HBB c.20A>T
# ----------------------------

def hbb_hbs_variant_if_present(contig: str, seq: str, chrom: Optional[str], gstart0: Optional[int]) -> Variant:
    if chrom != "chr11" or gstart0 is None:
        raise ValueError("HbS requested but genomic chrom/start are unavailable (need chr11 + genomic start).")

    # ClinVar HGVS includes: NC_000011.10:g.5227002T>A (GRCh38)
    pos1_genome = 5227002
    pos1_local = pos1_genome - gstart0  # gstart0 is 0-based locus start
    if pos1_local < 1 or pos1_local > len(seq):
        raise ValueError(
            f"HbS requested but site not within locus: genome pos1 {pos1_genome}, "
            f"locus spans {chrom}:{gstart0}-{gstart0+len(seq)} (0-based start)."
        )

    refb = seq[pos1_local - 1]
    if refb != "T":
        raise ValueError(
            f"Expected reference base T at HbS genomic site (g.5227002T>A), but saw {refb} in {contig} "
            f"(local pos1 {pos1_local}, genome pos1 {pos1_genome})."
        )

    # Genomic: T>A ; Transcript: c.20A>T (same biology, strand-aware representation)
    return Variant(
        pos1=pos1_local,
        ref="T",
        alt="A",
        vid="rs334",
        info="HBB_HbS_g.5227002T>A;c.20A>T"
    )

# ----------------------------
# Random variant generation (deterministic RNG)
# ----------------------------

def choose_positions(seq: str, n: int, rng: random.Random, margin: int,
                     avoid: Optional[List[Interval]] = None,
                     require: Optional[List[Interval]] = None) -> List[int]:
    L = len(seq)
    candidates: List[int] = []
    for pos0 in range(margin, L - margin):
        if seq[pos0] not in DNA:
            continue
        if avoid and in_any_interval(pos0, avoid):
            continue
        if require and (not in_any_interval(pos0, require)):
            continue
        candidates.append(pos0)

    if len(candidates) < n:
        raise ValueError(f"Not enough candidate positions (need {n}, have {len(candidates)})")
    return sorted(rng.sample(candidates, n))

def make_random_snps(seq: str, positions0: List[int], rng: random.Random, id_prefix: str) -> List[Variant]:
    vars: List[Variant] = []
    for i, pos0 in enumerate(positions0, 1):
        refb = seq[pos0]
        altb = mutate_base(refb, rng)
        vars.append(Variant(pos1=pos0+1, ref=refb, alt=altb, vid=f"{id_prefix}_snp{i}", info="."))
    return vars

def make_single_insertion(seq: str, pos1: int, ins: str, vid: str) -> Variant:
    # VCF insertion: REF = base at pos1; ALT = REF + inserted sequence
    refb = seq[pos1-1]
    return Variant(pos1=pos1, ref=refb, alt=refb + ins.upper(), vid=vid, info="INS")

def make_single_deletion(seq: str, pos1: int, del_len: int, vid: str) -> Variant:
    # VCF deletion: REF spans anchor + deleted bases; ALT = anchor base
    ref = seq[pos1-1:pos1-1+1+del_len]
    if len(ref) != 1 + del_len:
        raise ValueError("Deletion would run off end of contig")
    return Variant(pos1=pos1, ref=ref, alt=ref[0], vid=vid, info="DEL")


# ----------------------------
# Illumina-like read simulation
# ----------------------------

def simulate_read(seq: str, start0: int, readlen: int, rng: random.Random,
                  q_start: int, q_end: int, q_sd: float,
                  q_min: int, q_max: int) -> tuple[str, str]:
    frag = seq[start0:start0+readlen]
    bases = list(frag)
    quals: List[str] = []

    for i, b in enumerate(bases):
        if readlen == 1:
            q_mu = q_start
        else:
            q_mu = q_start + (q_end - q_start) * (i / (readlen - 1))

        q = int(round(rng.gauss(q_mu, q_sd)))
        q = clamp_int(q, q_min, q_max)

        # substitution errors based on Q
        if b in DNA and rng.random() < q_to_p_err(q):
            bases[i] = mutate_base(b, rng)

        quals.append(phred_char(q))

    return "".join(bases), "".join(quals)

def write_fastq_gz(path: Path, reads: List[tuple[str, str]], name_prefix: str) -> None:
    with gzip.open(path, "wt") as out:
        for i, (seq, qual) in enumerate(reads, 1):
            out.write(f"@{name_prefix}{i}\n{seq}\n+\n{qual}\n")


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser()

    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="data/loci")
    ap.add_argument("--prefix", default=None, help="Default = FASTA header name")
    ap.add_argument("--sample", default="ind1", help="Sample label used in filenames + read names")

    # Determinism
    ap.add_argument("--seed", type=int, default=1)

    # Genomic mapping (for HbS placement)
    ap.add_argument("--genomic-chrom", default=None, help="e.g., chr11 (optional; inferred from contig name if possible)")
    ap.add_argument("--genomic-start", type=int, default=None, help="0-based locus start (optional; inferred from contig name)")
    ap.add_argument("--add-hbb-hbs", action="store_true", help="Explicitly add HbS rs334 if locus supports it")

    # Repeats (local coords; for choosing where random SNPs go)
    ap.add_argument("--repeats-local-bed", default=None, help="Local BED in mini-contig coords (0..len)")
    ap.add_argument("--n_snps", type=int, default=0, help="Number of random SNPs to add (besides HbS/indels)")
    ap.add_argument("--n_snps_in_repeats", type=int, default=0, help="How many of the random SNPs should be inside repeats")
    ap.add_argument("--margin", type=int, default=150, help="Avoid placing variants near ends of locus")

    # Optional deterministic indels (repeatable)
    ap.add_argument("--ins", action="append", default=[], help="Insertion as POS1:SEQ (e.g., 500:AT)")
    ap.add_argument("--del", dest="dels", action="append", default=[], help="Deletion as POS1:LEN (e.g., 800:3)")

    # Zygosity-ish knob: fraction of reads drawn from ALT haplotype (0..1)
    ap.add_argument("--alt-fraction", type=float, default=1.0,
                    help="Fraction of reads from ALT haplotype. Use 0.5 to mimic heterozygous variants.")

    # Reads
    ap.add_argument("--paired", action="store_true")
    ap.add_argument("--readlen", type=int, default=100)
    ap.add_argument("--n_reads", type=int, default=2000, help="Reads (SE) or read-pairs (PE)")
    ap.add_argument("--insert-mean", type=int, default=350)
    ap.add_argument("--insert-sd", type=int, default=30)

    # Quality model
    ap.add_argument("--q-start", type=int, default=35)
    ap.add_argument("--q-end", type=int, default=18)
    ap.add_argument("--q-sd", type=float, default=3.0)
    ap.add_argument("--q-min", type=int, default=2)
    ap.add_argument("--q-max", type=int, default=40)

    args = ap.parse_args()

    if not (0.0 <= args.alt_fraction <= 1.0):
        raise ValueError("--alt-fraction must be between 0 and 1")

    rng = random.Random(args.seed)

    contig_id, ref = read_single_fasta(Path(args.fasta))
    out_prefix = args.prefix or contig_id
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Infer genomic window if not provided
    inferred_chrom, inferred_start, _ = infer_genomic_window_from_name(contig_id)
    gchrom = args.genomic_chrom or inferred_chrom
    gstart0 = args.genomic_start if args.genomic_start is not None else inferred_start

    # Load repeats (local coords)
    repeats: List[Interval] = []
    if args.repeats_local_bed:
        repeats = parse_local_bed_intervals(Path(args.repeats_local_bed), contig_filter=contig_id)

    # Build truth variants list
    truth_vars: List[Variant] = []

    if args.add_hbb_hbs:
        truth_vars.append(hbb_hbs_variant_if_present(contig_id, ref, gchrom, gstart0))

    for j, spec in enumerate(args.ins, 1):
        pos_s, seq_s = spec.split(":", 1)
        truth_vars.append(make_single_insertion(ref, int(pos_s), seq_s, vid=f"toyINS{j}"))

    for j, spec in enumerate(args.dels, 1):
        pos_s, len_s = spec.split(":", 1)
        truth_vars.append(make_single_deletion(ref, int(pos_s), int(len_s), vid=f"toyDEL{j}"))

    # Random SNPs: optionally push some into repeats
    if args.n_snps:
        n_in = max(0, min(args.n_snps_in_repeats, args.n_snps))
        n_out = args.n_snps - n_in

        occupied0 = {v.pos1 - 1 for v in truth_vars}

        pos_in0: List[int] = []
        if n_in:
            if not repeats:
                raise ValueError("--n-snps-in-repeats > 0 requires --repeats-local-bed")
            pos_in0 = choose_positions(ref, n_in, rng, margin=args.margin, require=repeats)

        pos_out0: List[int] = []
        if n_out:
            pos_out0 = choose_positions(ref, n_out, rng, margin=args.margin, avoid=repeats if repeats else None)

        chosen0: List[int] = []
        for p in (pos_in0 + pos_out0):
            if p not in occupied0:
                chosen0.append(p)
                occupied0.add(p)

        if len(chosen0) != args.n_snps:
            raise RuntimeError("Variant position collision; rerun with different --seed.")

        truth_vars.extend(make_random_snps(ref, chosen0, rng, id_prefix=f"{args.sample}"))

    # Apply variants onto ALT haplotype (ascending pos1)
    truth_vars.sort(key=lambda v: v.pos1)
    alt = ref
    for v in truth_vars:
        alt = apply_variant(alt, v)

    # Simulate reads (sample from REF vs ALT according to alt_fraction)
    L_ref = len(ref)
    L_alt = len(alt)
    if L_alt < args.readlen + 1 or L_ref < args.readlen + 1:
        raise ValueError("Locus shorter than read length.")

    reads_r1: List[tuple[str, str]] = []
    reads_r2: List[tuple[str, str]] = []

    def pick_template() -> str:
        return alt if rng.random() < args.alt_fraction else ref

    if not args.paired:
        for _ in range(args.n_reads):
            templ = pick_template()
            start0 = rng.randint(0, len(templ) - args.readlen)
            rseq, rqual = simulate_read(
                templ, start0, args.readlen, rng,
                q_start=args.q_start, q_end=args.q_end, q_sd=args.q_sd,
                q_min=args.q_min, q_max=args.q_max
            )
            reads_r1.append((rseq, rqual))
    else:
        min_insert = max(2 * args.readlen + 10, 2 * args.readlen)
        for _ in range(args.n_reads):
            templ = pick_template()

            insert = int(round(rng.gauss(args.insert_mean, args.insert_sd)))
            insert = max(min_insert, insert)
            insert = min(insert, len(templ))
            start0 = rng.randint(0, len(templ) - insert)
            frag = templ[start0:start0+insert]

            r1_seq, r1_qual = simulate_read(
                frag, 0, args.readlen, rng,
                q_start=args.q_start, q_end=args.q_end, q_sd=args.q_sd,
                q_min=args.q_min, q_max=args.q_max
            )
            r2_seq, r2_qual = simulate_read(
                rc(frag), 0, args.readlen, rng,
                q_start=args.q_start, q_end=args.q_end, q_sd=args.q_sd,
                q_min=args.q_min, q_max=args.q_max
            )

            reads_r1.append((r1_seq, r1_qual))
            reads_r2.append((r2_seq, r2_qual))

    # Write outputs
    tag = f"{out_prefix}.{args.sample}"
    vcf_path = outdir / f"{tag}.truth.vcf"
    r1_path  = outdir / f"{tag}.reads_R1.fastq.gz"
    write_vcf(vcf_path, contig=contig_id, contig_len=len(ref), vars=truth_vars)
    write_fastq_gz(r1_path, reads_r1, name_prefix=f"{args.sample}_R1_")

    if args.paired:
        r2_path = outdir / f"{tag}.reads_R2.fastq.gz"
        write_fastq_gz(r2_path, reads_r2, name_prefix=f"{args.sample}_R2_")

    # Summary
    print(f"FASTA:  {args.fasta} (contig={contig_id}, len(ref)={len(ref)}, len(alt)={len(alt)})")
    if gchrom and gstart0 is not None:
        print(f"Genomic mapping (inferred/used): {gchrom}:{gstart0}-{gstart0+len(ref)} (0-based start)")
    if args.repeats_local_bed:
        print(f"Repeats local BED: {args.repeats_local_bed}")
    print(f"ALT fraction: {args.alt_fraction}")
    print(f"Truth:  {vcf_path}")
    print(f"Reads:  {r1_path}" + (f" + {outdir / f'{tag}.reads_R2.fastq.gz'}" if args.paired else ""))
    for v in truth_vars:
        print(f"  VAR {v.vid}: pos1={v.pos1} {v.ref}>{v.alt} info={v.info}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
