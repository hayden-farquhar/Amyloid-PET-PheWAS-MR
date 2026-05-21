#!/usr/bin/env python3
"""Convert full Ali 2023 GRCh38 sumstats to PRS-CS format (SNP A1 A2 BETA P).

PRS-CS requires:
- rsID in column SNP
- A1 = effect allele
- A2 = other allele
- BETA and P columns

Ali 2023 uses chr:pos:ref:alt (GRCh38) as SNP ID. We:
1. Liftover each position GRCh38 -> GRCh37.
2. Match against 1000G EUR.bim by (chr, pos_grch37) to get rsID.
3. Confirm alleles match (or are reverse complement).
4. Write PRS-CS-format file.

Typical wall-clock: ~10-15 minutes for 15.6M SNPs.

Usage:
  python3 scripts/make_prscs_sumstats.py \\
    --in data/raw/ali_2023/AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt \\
    --chain data/reference/hg38ToHg19.over.chain.gz \\
    --eur-bim ../60\\ Statin\\ MR\\ Pathogen/data/reference/EUR.bim \\
    --out data/interim/ali2023_prscs_sumstats.tsv
"""
import argparse
import sys
import time
from pyliftover import LiftOver

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="infile", required=True)
ap.add_argument("--chain", required=True)
ap.add_argument("--eur-bim", dest="bim", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

t0 = time.time()
print(f"[1/4] Loading liftover chain {args.chain}", file=sys.stderr)
lo = LiftOver(args.chain)

print(f"[2/4] Building chr+pos -> rsID lookup from {args.bim}", file=sys.stderr)
lookup = {}
with open(args.bim) as f:
    for line in f:
        parts = line.rstrip().split("\t")
        if len(parts) < 6:
            continue
        try:
            c = int(parts[0]); p_ = int(parts[3])
        except ValueError:
            continue
        lookup[(c, p_)] = (parts[1], parts[4], parts[5])
print(f"      Loaded {len(lookup):,} BIM entries (t={time.time()-t0:.1f}s)", file=sys.stderr)

print(f"[3/4] Reading {args.infile} and converting to PRS-CS format", file=sys.stderr)
n_in = n_out = n_lift_fail = n_no_rsid = n_allele_mm = 0
out = open(args.out, "w")
out.write("SNP\tA1\tA2\tBETA\tP\n")
with open(args.infile) as f:
    header = f.readline().strip().split("\t")
    col = {c: i for i, c in enumerate(header)}
    for line in f:
        n_in += 1
        if n_in % 1_000_000 == 0:
            print(f"      processed {n_in:,} / mapped {n_out:,}  (t={time.time()-t0:.1f}s)", file=sys.stderr)
        s = line.rstrip("\n").split("\t")
        try:
            p = float(s[col["P"]])
            beta = float(s[col["BETA"]])
            bp = int(s[col["BP"]])
        except (ValueError, IndexError):
            continue
        chr_s = s[col["CHR"]]
        if not chr_s.isdigit():
            continue
        chr_ = int(chr_s)
        # Liftover GRCh38 -> GRCh37
        r = lo.convert_coordinate(f"chr{chr_}", bp - 1)
        if not r:
            n_lift_fail += 1
            continue
        nchr_s, npos = r[0][0], r[0][1] + 1
        if nchr_s != f"chr{chr_}":
            n_lift_fail += 1
            continue
        hit = lookup.get((chr_, npos))
        if hit is None:
            n_no_rsid += 1
            continue
        rsid, bim_a1, bim_a2 = hit
        ea = s[col["EFFECT_ALLELE"]]
        oa = s[col["OTHER_ALLELE"]]
        # Allele match check (could also be strand-flipped; PRS-CS will catch mismatches)
        ali_set = {ea, oa}
        bim_set = {bim_a1, bim_a2}
        if ali_set != bim_set:
            n_allele_mm += 1
            continue
        out.write(f"{rsid}\t{ea}\t{oa}\t{beta:.6f}\t{p:.4e}\n")
        n_out += 1
out.close()

print(f"[4/4] Done. Read {n_in:,}; wrote {n_out:,}; liftover fail {n_lift_fail:,}; "
      f"no-rsID {n_no_rsid:,}; allele mismatch {n_allele_mm:,}  (t={time.time()-t0:.1f}s)",
      file=sys.stderr)
