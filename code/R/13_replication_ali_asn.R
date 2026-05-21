# =============================================================================
# Script 13: Phase 5 — Cross-ancestry replication (Ali 2023 ASN arm)
# =============================================================================
#
# Per Amendment 3 (osf/amendment_log.md, pushed to OSF 2026-05-20), the H4
# hypothesis is now tested against the Ali 2023 ASN arm (N=1,494, East Asian
# amyloid-PET multi-ethnic meta-analysis subset) instead of BBJ (BBJ does not
# host any AD/dementia/amyloid GWAS suitable for the replication).
#
# H4 (revised): "≥1 non-APOE hit replicates in the Ali 2023 East Asian arm at
# nominal p < 0.05, OR in the Kim 2025 cross-ancestry GWAS if/when sumstats
# become available from the Won lab."
#
# This script performs the **instrument-portability** replication: for each of
# the EUR-identified IVs in our 16-IV unified panel, look up its association
# with amyloid-PET in the Ali ASN arm, and test for direction concordance +
# per-IV p-value.
#
# Note: a phecode-MR cross-ancestry replication would require EAS outcome
# GWASes for the FinnGen phecodes that show signal, which don't exist in the
# major repositories (see §15.3 — MR-MEGA deferred). The instrument-portability
# test is the realisable form of H4.
#
# Prerequisites:
#   1. Re-download the Ali 2023 ASN arm from the WashU Box folder:
#      https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l
#      File: AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt
#      Save to: data/raw/ali_2023/
#      (The file was deleted locally on 2026-05-20 for disk space; SHA256 is
#       preserved in data/raw/ali_2023/SHA256SUMS for integrity verification.)
#
#   2. Phase 2 must be complete (we need the EUR hit list to know which IVs
#      to spotlight, although the script also reports on all 16 unified IVs).
#
# Input:
#   - data/raw/ali_2023/AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt
#   - data/interim/instruments/unified_iv_list.tsv
#
# Output:
#   - data/processed/replication_ali_asn.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")

ALI_ASN_FP <- file.path(proj_root, "data/raw/ali_2023",
  "AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt")

if (!file.exists(ALI_ASN_FP)) {
  cat("\n>>> Ali 2023 ASN arm not found at:\n>>> ", ALI_ASN_FP, "\n\n")
  cat("Re-download required:\n")
  cat("  1. Open in browser: https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l\n")
  cat("  2. Download: AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt\n")
  cat("  3. Save to: data/raw/ali_2023/\n")
  cat("  4. Verify SHA256 against data/raw/ali_2023/SHA256SUMS\n")
  cat("  5. Re-run this script.\n\n")
  stop("ASN arm missing")
}

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(dir_interim, "instruments/unified_iv_list.tsv"))
cat(sprintf("[%s] Unified IVs (EUR-identified): %d\n", Sys.time(), nrow(iv_unified)))

# Read Ali ASN sumstats
# Expected schema (matches Ali 2023 multi-ethnic file format):
#   CHR, BP, EFFECT_ALLELE, OTHER_ALLELE, BETA, SE, P, ...
cat(sprintf("[%s] Reading Ali ASN sumstats...\n", Sys.time()))
asn <- fread(ALI_ASN_FP)
cat(sprintf("  ASN sumstats: %d SNPs\n", nrow(asn)))

# Extract just the IV positions
ivpos <- iv_unified[, .(chr = as.integer(chr), pos = as.integer(hg38_pos),
                        rsid, ea_eur = effect_allele, oa_eur = other_allele,
                        beta_eur = beta, se_eur = se, p_eur = p, eaf_eur = eaf,
                        f_stat_eur = f_stat, panels = panels)]

# Match by chr+pos
m <- merge(ivpos,
           asn[, .(CHR, BP, EA_ASN = EFFECT_ALLELE, OA_ASN = OTHER_ALLELE,
                   beta_asn = BETA, se_asn = SE, p_asn = P,
                   eaf_asn = if ("EAF" %in% names(asn)) EAF else NA_real_)],
           by.x = c("chr", "pos"), by.y = c("CHR", "BP"))

cat(sprintf("[%s] Matched %d / %d unified IVs in Ali ASN\n",
            Sys.time(), nrow(m), nrow(iv_unified)))

if (nrow(m) == 0) stop("No IVs matched in Ali ASN sumstats — check schema")

# --- Harmonise effect alleles -----------------------------------------------

# Flip ASN beta if alleles are swapped vs EUR
m[, flip := FALSE]
m[toupper(ea_eur) == toupper(OA_ASN) & toupper(oa_eur) == toupper(EA_ASN), flip := TRUE]
m[flip == TRUE, beta_asn := -beta_asn]

# Drop SNPs where alleles mismatch entirely (neither orientation works)
m_clean <- m[(toupper(ea_eur) == toupper(EA_ASN) & toupper(oa_eur) == toupper(OA_ASN)) |
             (toupper(ea_eur) == toupper(OA_ASN) & toupper(oa_eur) == toupper(EA_ASN))]
cat(sprintf("[%s] Post-harmonise: %d / %d IVs (dropped %d for allele mismatch)\n",
            Sys.time(), nrow(m_clean), nrow(m), nrow(m) - nrow(m_clean)))

# --- Concordance + Z statistic -----------------------------------------------

m_clean[, concordant_direction := sign(beta_eur) == sign(beta_asn)]
m_clean[, replication_z := beta_asn / se_asn]
m_clean[, replication_pval := 2 * pnorm(-abs(replication_z))]

# A directional replication is significant if:
#   - Same direction as EUR
#   - Replication p < 0.05
m_clean[, h4_replicates := concordant_direction & replication_pval < 0.05]

# --- Per-panel summary ------------------------------------------------------

panels <- c("tier1_apoe_inclusive", "tier1_apoe_narrow_excl",
            "tier1_apoe_wide_excl", "tier2_apoe_narrow_excl", "tier2_apoe_wide_excl")

panel_summary <- list()
for (pn in panels) {
  pn_ivs <- m_clean[grepl(pn, panels, fixed = TRUE)]
  panel_summary[[length(panel_summary) + 1]] <- data.table(
    panel = pn,
    n_iv_eur = nrow(pn_ivs),
    n_concordant_direction = sum(pn_ivs$concordant_direction, na.rm = TRUE),
    frac_concordant = mean(pn_ivs$concordant_direction, na.rm = TRUE),
    n_h4_replicates = sum(pn_ivs$h4_replicates, na.rm = TRUE),
    min_repl_pval = min(pn_ivs$replication_pval, na.rm = TRUE)
  )
}
panel_summary_dt <- rbindlist(panel_summary)

# --- Write -----------------------------------------------------------------

write_parquet(m_clean,
              file.path(dir_processed, "replication_ali_asn_per_iv.parquet"))
write_parquet(panel_summary_dt,
              file.path(dir_processed, "replication_ali_asn_per_panel.parquet"))
cat(sprintf("[%s] Wrote replication_ali_asn_*.parquet\n", Sys.time()))

# --- Report ----------------------------------------------------------------

cat("\n=== Per-IV replication (EUR vs Ali ASN N=1,494) ===\n")
print(m_clean[order(replication_pval),
              .(rsid, chr, pos, beta_eur = round(beta_eur, 3),
                beta_asn = round(beta_asn, 3),
                concordant = concordant_direction,
                p_eur = signif(p_eur, 3),
                p_asn = signif(replication_pval, 3),
                h4_replicates)])

cat("\n=== Per-panel replication summary ===\n")
print(panel_summary_dt)

cat(sprintf("\n=== H4 replication outcome (Amendment 3) ===\n"))
n_total_replicates <- sum(m_clean$h4_replicates, na.rm = TRUE)
cat(sprintf("Total IVs replicating in Ali ASN at p<0.05 + same direction: %d / %d\n",
            n_total_replicates, nrow(m_clean)))
non_apoe_replicates <- m_clean[grepl("wide_excl", panels) & h4_replicates == TRUE]
cat(sprintf("Non-APOE IVs (in any wide-excl panel) replicating: %d\n",
            nrow(non_apoe_replicates)))
if (nrow(non_apoe_replicates) > 0) {
  cat("\nH4 HYPOTHESIS SUPPORTED — non-APOE IVs replicate cross-ancestry.\n")
  print(non_apoe_replicates[, .(rsid, beta_eur, beta_asn,
                                  p_asn = replication_pval)])
} else {
  cat("\nH4 hypothesis NOT supported — no non-APOE IVs replicate at p<0.05 in EAS.\n")
  cat("(Power-limited: ASN N=1,494 is much smaller than EUR N=13,409.)\n")
}
