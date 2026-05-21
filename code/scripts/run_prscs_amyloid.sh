#!/usr/bin/env bash
# Run PRS-CS to derive the amyloid-PET polygenic score from Ali 2023 sumstats.
# Generates the Tier 3 PGS-MR instrument (now PRIMARY under Pivot A activation).
#
# PREREQUISITES (to be set up in next session):
#   1. Download PRS-CS 1KG EUR LD reference panel (~6 GB):
#      mkdir -p data/reference/PRScs_ld_ref
#      cd data/reference/PRScs_ld_ref
#      wget https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/1KG/ldblk_1kg_eur.tar.gz
#      tar -xzf ldblk_1kg_eur.tar.gz
#      cd ../../..
#   2. Sumstats in PRS-CS format: SNP A1 A2 BETA P (rsid effect_allele other_allele beta p_value)
#      Generated below from the full Ali 2023 sumstats via the make_prscs_sumstats.py helper.
#   3. BIM file from a validation cohort (use the 1000G EUR.bim as proxy validation target).
#
# OUTPUT: per-SNP posterior effect sizes (one file per chromosome, named like
#         pst_eff_a1_b0.5_phiauto_chr*.txt), which collectively constitute the
#         PGS instrument. Combine into a single PGS via PLINK --score.

set -euo pipefail

REF_DIR="data/reference/PRScs_ld_ref/ldblk_1kg_eur"
SST_FILE="data/interim/ali2023_prscs_sumstats.tsv"
BIM_PREFIX="${EUR_LD_REFERENCE_PREFIX:-./data/reference/EUR}"
OUT_DIR="data/interim/prscs"
N_GWAS=13409

if [[ ! -d "$REF_DIR" ]]; then
  echo "ERROR: PRS-CS LD reference not found at $REF_DIR"
  echo "Download instructions in script header."
  exit 1
fi

mkdir -p "$OUT_DIR"

# Generate PRS-CS-format sumstats: SNP A1 A2 BETA P
# Need to convert the FULL Ali 2023 file (15.6M SNPs) to this format,
# not just the significant ones - PRS-CS uses ALL SNPs to estimate posteriors.
echo "Generating PRS-CS-format sumstats from full Ali 2023 file..."
.venv/bin/python3 scripts/make_prscs_sumstats.py \
  --in data/raw/ali_2023/AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt \
  --chain data/reference/hg38ToHg19.over.chain.gz \
  --eur-bim "$BIM_PREFIX.bim" \
  --out "$SST_FILE"

# Run PRS-CS per chromosome (parallelizable via xargs/parallel)
for chr in {1..22}; do
  .venv/bin/python3 data/reference/PRScs/PRScs.py \
    --ref_dir="$REF_DIR" \
    --bim_prefix="$BIM_PREFIX" \
    --sst_file="$SST_FILE" \
    --n_gwas=$N_GWAS \
    --chrom=$chr \
    --phi=1e-2 \
    --out_dir="$OUT_DIR/amyloid_pgs" \
    --seed=81
done

echo "PRS-CS complete. Combine posterior effect sizes into a single PGS via:"
echo "  cat $OUT_DIR/amyloid_pgs*_pst_eff_a1_b0.5_phi*_chr*.txt > $OUT_DIR/amyloid_pgs_combined.txt"
echo "Then use as PLINK --score input to compute per-individual PGS in any validation cohort."
