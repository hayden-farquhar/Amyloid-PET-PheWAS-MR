# =============================================================================
# Script 13b: Amendment 4 Component A — Kim Korean cross-ancestry replication
# =============================================================================
#
# Per Amendment 4 Component A: replaces Ali ASN (Amendment 3) with Kim Korean
# (GCST90483382, N=3,387 EAS) as the primary East Asian cross-ancestry
# replication reference. Ali ASN is retained as secondary sensitivity (handled
# by R/13_replication_ali_asn.R).
#
# H4 (revised): "≥1 non-APOE hit replicates in the Kim 2025 Korean amyloid-PET
# cohort (GCST90483382, N=3,387) at nominal p < 0.05 in concordant direction
# with the EUR discovery effect."
#
# Method: instrument-portability test — for each of the 16 unified IVs from
# the Ali EUR primary, look up the corresponding effect in the Kim Korean
# sumstats and test for direction concordance + per-IV replication p-value.
#
# This script is materially smaller than R/15 because it doesn't run a full
# locus-MR sweep — it just performs the per-IV portability table that
# determines H4 outcome.
#
# Input:
#   - data/processed/kim_korean_iv_lookup.parquet (produced by R/15)
#     OR data/reference/kim_2025/GCST90483382.tsv (if lookup parquet absent)
#   - data/interim/instruments/unified_iv_list.tsv
#
# Output:
#   - data/processed/replication_kim_korean.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
KIM_FP <- file.path(proj_root, "data/reference/kim_2025/GCST90483382.tsv")

# --- Get the per-IV Kim Korean lookup ----------------------------------------

lookup_fp <- file.path(dir_processed, "kim_korean_iv_lookup.parquet")
if (file.exists(lookup_fp)) {
  cat(sprintf("[%s] Reading pre-built IV lookup from %s\n",
              Sys.time(), basename(lookup_fp)))
  m <- as.data.table(read_parquet(lookup_fp))
} else {
  cat(sprintf("[%s] No lookup parquet — building inline from Kim Korean sumstats\n",
              Sys.time()))
  iv_unified <- fread(file.path(dir_interim, "instruments/unified_iv_list.tsv"))
  kim <- fread(KIM_FP)
  setnames(kim, "rs_id", "rsid")
  m <- merge(iv_unified[, .(rsid, chr, hg38_pos, ea_eur = effect_allele,
                             oa_eur = other_allele, beta_eur = beta,
                             se_eur = se, p_eur = p, panels)],
             kim[, .(rsid, chr_kim = chromosome, pos_kim = base_pair_location,
                     ea_kim = effect_allele, oa_kim = other_allele,
                     beta_kim = beta, se_kim = standard_error,
                     eaf_kim = effect_allele_frequency, p_kim = p_value)],
             by = "rsid")
  m[, flip := FALSE]
  m[toupper(ea_eur) == toupper(oa_kim) & toupper(oa_eur) == toupper(ea_kim), flip := TRUE]
  m[flip == TRUE, beta_kim := -beta_kim]
  m <- m[(toupper(ea_eur) == toupper(ea_kim) & toupper(oa_eur) == toupper(oa_kim)) |
          (toupper(ea_eur) == toupper(oa_kim) & toupper(oa_eur) == toupper(ea_kim))]
}
cat(sprintf("  IVs available: %d\n", nrow(m)))

# --- Concordance + replication statistics ------------------------------------

m[, concordant_direction := sign(beta_eur) == sign(beta_kim)]
m[, replication_z := beta_kim / se_kim]
m[, replication_pval := 2 * pnorm(-abs(replication_z))]
m[, h4_replicates := concordant_direction & replication_pval < 0.05]

# --- Per-panel summary ------------------------------------------------------

panels <- c("tier1_apoe_inclusive", "tier1_apoe_narrow_excl",
            "tier1_apoe_wide_excl", "tier2_apoe_narrow_excl",
            "tier2_apoe_wide_excl")

panel_summary <- list()
for (pn in panels) {
  pn_ivs <- m[grepl(pn, panels, fixed = TRUE)]
  panel_summary[[length(panel_summary) + 1]] <- data.table(
    panel = pn,
    n_iv_eur = nrow(pn_ivs),
    n_concordant_direction = sum(pn_ivs$concordant_direction, na.rm = TRUE),
    frac_concordant = if (nrow(pn_ivs) > 0) mean(pn_ivs$concordant_direction, na.rm = TRUE) else NA_real_,
    n_h4_replicates = sum(pn_ivs$h4_replicates, na.rm = TRUE),
    min_repl_pval = if (nrow(pn_ivs) > 0) min(pn_ivs$replication_pval, na.rm = TRUE) else NA_real_
  )
}
panel_summary_dt <- rbindlist(panel_summary)

# --- Write ------------------------------------------------------------------

write_parquet(m, file.path(dir_processed, "replication_kim_korean_per_iv.parquet"))
write_parquet(panel_summary_dt, file.path(dir_processed, "replication_kim_korean_per_panel.parquet"))
cat(sprintf("[%s] Wrote replication_kim_korean_*.parquet\n", Sys.time()))

# --- Report ----------------------------------------------------------------

cat("\n=== Per-IV replication (EUR vs Kim Korean N=3,387) ===\n")
# p_eur may be absent from older lookup parquets; render only what exists
cols_to_print <- c("rsid", "chr", "beta_eur", "beta_kim", "concordant_direction",
                    if ("p_eur" %in% names(m)) "p_eur",
                    "replication_pval", "h4_replicates")
print(m[order(replication_pval), ..cols_to_print])

cat("\n=== Per-panel replication summary ===\n")
print(panel_summary_dt)

cat(sprintf("\n=== H4 outcome (Amendment 4 Component A) ===\n"))
n_total_replicates <- sum(m$h4_replicates, na.rm = TRUE)
cat(sprintf("Total IVs replicating at p<0.05 + same direction: %d / %d\n",
            n_total_replicates, nrow(m)))
non_apoe_replicates <- m[grepl("wide_excl", panels) & h4_replicates == TRUE]
cat(sprintf("Non-APOE IVs (in any wide-excl panel) replicating: %d\n",
            nrow(non_apoe_replicates)))
if (nrow(non_apoe_replicates) > 0) {
  cat("\nH4 HYPOTHESIS SUPPORTED — non-APOE IVs replicate in Kim Korean.\n")
  print(non_apoe_replicates[, .(rsid, beta_eur, beta_kim,
                                  p_kim = replication_pval)])
} else {
  cat("\nH4 hypothesis NOT supported in Kim Korean at p<0.05 — non-APOE IVs underpowered or non-replicating.\n")
}
