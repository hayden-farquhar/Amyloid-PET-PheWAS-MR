#!/usr/bin/env python3
"""Build the wide PGS-MR non-APOE IV list with GRCh38 positions.

PRS-CS uses the 1000G EUR.bim positions (GRCh37) from P60. FinnGen R12 is
GRCh38. We lift PRS-CS positions GRCh37→GRCh38 so the sweep_finngen_r12.sh
tabix queries hit the right coordinates.

The wide-APOE LD block (chr19:44.4-46.5 Mb in GRCh38) is excluded BEFORE
liftover, using a GRCh37-equivalent boundary (chr19:44.9-47 Mb in GRCh37 —
roughly +500 kb offset). The exclusion is generous (2.1 Mb wide); the
build mismatch (~500 kb at this locus) doesn't change which SNPs are
excluded in practice.

Output: data/interim/instruments/pgs_wide_iv_list_noAPOE.tsv
Columns: rsid, chr, hg38_pos, A1, A2, beta_pgs, hg37_pos (for audit)
"""
import sys, os
from pathlib import Path
from pyliftover import LiftOver
import csv

PROJ_ROOT = Path(os.environ.get("AMYLOID_PHEWASMR_ROOT", os.getcwd()))
os.chdir(PROJ_ROOT)

N_TOP = int(sys.argv[1]) if len(sys.argv) >= 2 else 50

# Load all PRS-CS chr files
prscs = []
for chr_file in sorted((PROJ_ROOT / "data/interim/prscs").glob("amyloid_pgs_pst_eff_*_chr*.txt")):
    with open(chr_file) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            chrom, rsid, bp, A1, A2, beta_pgs = parts[:6]
            try:
                prscs.append({
                    "chr": int(chrom), "rsid": rsid, "bp_hg37": int(bp),
                    "A1": A1, "A2": A2, "beta_pgs": float(beta_pgs),
                })
            except ValueError:
                continue

print(f"Total PRS-CS SNPs loaded: {len(prscs):,}", file=sys.stderr)

# Exclude wide-APOE LD block (GRCh37 coords: chr19:44,907,000-47,007,000
# approximates the GRCh38 chr19:44.4-46.5 Mb block after +500 kb offset)
def in_apoe_grch37(d):
    return d["chr"] == 19 and 44907000 <= d["bp_hg37"] <= 47007000

prscs_noapoe = [d for d in prscs if not in_apoe_grch37(d)]
print(f"After wide-APOE exclusion (GRCh37 chr19:44.9-47.0 Mb): {len(prscs_noapoe):,}", file=sys.stderr)

# Top N by absolute beta_pgs
prscs_noapoe.sort(key=lambda d: -abs(d["beta_pgs"]))
top_n = prscs_noapoe[:N_TOP]

# Liftover GRCh37 → GRCh38
print(f"Liftover GRCh37→GRCh38 for top {N_TOP}...", file=sys.stderr)
lo = LiftOver("data/reference/hg19ToHg38.over.chain.gz")
for d in top_n:
    r = lo.convert_coordinate(f"chr{d['chr']}", d["bp_hg37"] - 1)  # 0-based
    if r and len(r) > 0 and r[0][0] == f"chr{d['chr']}":
        d["hg38_pos"] = r[0][1] + 1  # back to 1-based
    else:
        d["hg38_pos"] = None

n_lifted = sum(1 for d in top_n if d["hg38_pos"] is not None)
print(f"Lifted: {n_lifted} / {N_TOP}", file=sys.stderr)
print(f"  Weight range: {min(abs(d['beta_pgs']) for d in top_n):.2e} to "
      f"{max(abs(d['beta_pgs']) for d in top_n):.2e}", file=sys.stderr)

# Filter to successfully lifted + sort by chr, pos
top_n = [d for d in top_n if d["hg38_pos"] is not None]
top_n.sort(key=lambda d: (d["chr"], d["hg38_pos"]))

# Write
out_fp = PROJ_ROOT / "data/interim/instruments/pgs_wide_iv_list_noAPOE.tsv"
with open(out_fp, "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["rsid", "chr", "hg38_pos", "A1", "A2", "beta_pgs", "hg37_pos"])
    for d in top_n:
        w.writerow([d["rsid"], d["chr"], d["hg38_pos"],
                    d["A1"], d["A2"], f"{d['beta_pgs']:.7g}", d["bp_hg37"]])

print(f"Wrote {out_fp} ({len(top_n)} SNPs across "
      f"{len(set(d['chr'] for d in top_n))} chromosomes)", file=sys.stderr)
print("First 5 rows:", file=sys.stderr)
for d in top_n[:5]:
    print(f"  {d['rsid']} chr{d['chr']}:{d['hg38_pos']} (hg37 {d['bp_hg37']})  beta={d['beta_pgs']:+.4f}", file=sys.stderr)
