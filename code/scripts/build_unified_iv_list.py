#!/usr/bin/env python3
"""Build unified instrument-variant list across all 5 IV panels.

Output: data/interim/instruments/unified_iv_list.tsv
Schema: rsid, chr, hg38_pos, hg19_pos, effect_allele, other_allele, beta, se, p, eaf, f_stat, panels

The `panels` column is a comma-separated list of which panels include this rsID.
Used by scripts/sweep_finngen_r12.sh to filter FinnGen R12 sumstats.
Used by Phase 2 MR runs to select IVs per panel.
"""
from pathlib import Path

INSTR = Path("data/interim/instruments")
PANELS = {
    "tier1_apoe_inclusive":    "tier1_apoe_inclusive_index_snps.tsv",
    "tier1_apoe_narrow_excl":  "tier1_apoe_narrow_excl_index_snps.tsv",
    "tier1_apoe_wide_excl":    "tier1_apoe_wide_excl_index_snps.tsv",
    "tier2_apoe_narrow_excl":  "tier2_apoe_narrow_excl_index_snps.tsv",
    "tier2_apoe_wide_excl":    "tier2_apoe_wide_excl_index_snps.tsv",
}

rows = {}      # rsid -> first-seen row dict
membership = {}  # rsid -> set of panel names

for panel, fn in PANELS.items():
    fp = INSTR / fn
    with open(fp) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < len(header):
                continue
            row = dict(zip(header, cells))
            rsid = row["rsid"]
            if rsid not in rows:
                rows[rsid] = row
            membership.setdefault(rsid, set()).add(panel)

out_fp = INSTR / "unified_iv_list.tsv"
cols = ["rsid", "chr", "hg38_pos", "hg19_pos", "effect_allele", "other_allele",
        "beta", "se", "p", "eaf", "f_stat", "panels"]
with open(out_fp, "w") as out:
    out.write("\t".join(cols) + "\n")
    # Sort by chr (numeric), then bp (numeric)
    def sortkey(rsid):
        r = rows[rsid]
        try:
            return (int(r["chr"]), int(r["hg38_pos"]))
        except ValueError:
            return (99, 0)
    for rsid in sorted(rows.keys(), key=sortkey):
        r = rows[rsid]
        panels_str = ",".join(sorted(membership[rsid]))
        out.write("\t".join([
            r["rsid"], r["chr"], r["hg38_pos"], r["hg19_pos"],
            r["effect_allele"], r["other_allele"],
            r["beta"], r["se"], r["p"], r["eaf"], r["f_stat"],
            panels_str,
        ]) + "\n")

print(f"Wrote {out_fp} with {len(rows)} unique rsIDs across {len(PANELS)} panels.")
print()
print("Panel membership summary:")
from collections import Counter
counts = Counter()
for s in membership.values():
    counts[frozenset(s)] += 1
for panels, n in sorted(counts.items(), key=lambda x: -x[1]):
    print(f"  {n:3d}  {','.join(sorted(panels))}")
