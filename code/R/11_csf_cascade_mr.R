# =============================================================================
# Script 11: CSF biomarker cascade MR (amyloid-PET → CSF Aβ42 / t-tau / p-tau181)
# =============================================================================
#
# Tests whether genetically-predicted amyloid-PET burden causally affects three
# downstream CSF biomarkers (Timsina 2026, N=18,948 EUR). Pre-registered as the
# Amendment 2 multi-biomarker cascade deliverable.
#
# Biological expectation (a-priori):
#   amyloid-PET → CSF Aβ42:  NEGATIVE effect (more brain amyloid plaque → less
#                            soluble Aβ42 in CSF — classical AD biomarker pattern)
#   amyloid-PET → CSF t-tau: POSITIVE effect (downstream tau pathology cascade)
#   amyloid-PET → CSF p-tau181: POSITIVE effect (downstream tau pathology cascade)
#
# All-three concordant directions (− + +) strengthen the amyloid-causal-cascade
# interpretation. Discordant directions raise concerns about pleiotropy or
# unmodeled selection.
#
# Per-panel results enable: APOE-inclusive (known biology validation) vs
# non-APOE wide-excl (the question that actually matters for amyloid-direct
# causal claims).
#
# Pre-registration: Amendment 2 (osf/amendment_log.md, pushed to OSF 2026-05-20).
#
# Input:
#   - data/reference/timsina_2026/GCST90726396_tsv (Aβ42)
#   - data/reference/timsina_2026/GCST90726397_tsv (t-tau)
#   - data/reference/timsina_2026/GCST90726398_tsv (p-tau181)
#   - data/interim/instruments/unified_iv_list.tsv (16 SNPs × 5 panel memberships)
#
# Output:
#   - data/processed/csf_cascade_mr_results.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(arrow)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)

TIMSINA_FILES <- list(
  csf_abeta42 = "data/reference/timsina_2026/GCST90726396_tsv",
  csf_ttau    = "data/reference/timsina_2026/GCST90726397_tsv",
  csf_ptau181 = "data/reference/timsina_2026/GCST90726398_tsv"
)
N_timsina <- 18948

SEED <- 81
set.seed(SEED)

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(dir_interim, "instruments", "unified_iv_list.tsv"))
PANELS <- c("tier1_apoe_inclusive",
            "tier1_apoe_narrow_excl",
            "tier1_apoe_wide_excl",
            "tier2_apoe_narrow_excl",
            "tier2_apoe_wide_excl")

cat(sprintf("[%s] Unified IVs: %d; testing %d Timsina biomarkers × %d panels\n",
            Sys.time(), nrow(iv_unified), length(TIMSINA_FILES), length(PANELS)))

# --- Build outcome dat for each Timsina biomarker -----------------------------

# Timsina schema: chromosome, base_pair_location, effect_allele, other_allele,
# beta, standard_error, effect_allele_frequency, p_value
# All GRCh38, EUR ancestry.

build_outcome_dat <- function(timsina_fp, label) {
  cat(sprintf("[%s] Loading %s (Timsina 2026)\n", Sys.time(), label))
  d <- fread(timsina_fp)
  # Build a lookup keyed by chr:pos for efficient join against iv_unified
  setkey(d, chromosome, base_pair_location)
  # Subset to just the 16 IV positions
  ivpos <- iv_unified[, .(chr = as.integer(chr), pos = as.integer(hg38_pos),
                          rsid, effect_allele, other_allele)]
  m <- merge(ivpos,
             d[, .(chr = chromosome, pos = base_pair_location,
                   ea = effect_allele, oa = other_allele,
                   beta_out = beta, se_out = standard_error, p_out = p_value)],
             by = c("chr", "pos"))
  if (nrow(m) == 0) {
    cat(sprintf("  WARN: 0 IVs matched in %s\n", label))
    return(NULL)
  }
  # Harmonise: if alleles flipped, negate beta
  m[, flip := FALSE]
  m[toupper(effect_allele) == toupper(oa) & toupper(other_allele) == toupper(ea), flip := TRUE]
  m[flip == TRUE, beta_out := -beta_out]
  # Drop mismatched alleles (neither orientation matches)
  m <- m[(toupper(effect_allele) == toupper(ea) & toupper(other_allele) == toupper(oa)) |
         (toupper(effect_allele) == toupper(oa) & toupper(other_allele) == toupper(ea))]
  cat(sprintf("  IVs matched (post-harmonise): %d / 16\n", nrow(m)))
  data.frame(
    SNP                    = m$rsid,
    beta.outcome           = m$beta_out,
    se.outcome             = m$se_out,
    pval.outcome           = m$p_out,
    effect_allele.outcome  = m$effect_allele,
    other_allele.outcome   = m$other_allele,
    eaf.outcome            = NA_real_,
    outcome                = label,
    id.outcome             = label,
    stringsAsFactors       = FALSE
  )
}

outcome_dats <- lapply(seq_along(TIMSINA_FILES), function(i) {
  build_outcome_dat(TIMSINA_FILES[[i]], names(TIMSINA_FILES)[i])
})
names(outcome_dats) <- names(TIMSINA_FILES)

# --- Format exposure (per-panel) ---------------------------------------------

format_exposure <- function(iv_subset, panel) {
  data.frame(
    SNP                    = iv_subset$rsid,
    beta.exposure          = iv_subset$beta,
    se.exposure            = iv_subset$se,
    pval.exposure          = iv_subset$p,
    effect_allele.exposure = iv_subset$effect_allele,
    other_allele.exposure  = iv_subset$other_allele,
    eaf.exposure           = iv_subset$eaf,
    exposure               = "amyloid_PET",
    id.exposure            = panel,
    stringsAsFactors       = FALSE
  )
}

# --- Run MR per panel × biomarker -------------------------------------------

results <- list()

for (panel in PANELS) {
  iv_sub <- iv_unified[grepl(panel, panels, fixed = TRUE)]
  if (nrow(iv_sub) == 0) next
  exposure_dat <- format_exposure(iv_sub, panel)

  for (biomarker in names(TIMSINA_FILES)) {
    outcome_dat <- outcome_dats[[biomarker]]
    if (is.null(outcome_dat) || nrow(outcome_dat) == 0) next

    harm <- tryCatch(
      TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2),
      error = function(e) NULL
    )
    if (is.null(harm) || sum(harm$mr_keep, na.rm = TRUE) == 0) next
    harm <- harm[harm$mr_keep == TRUE, ]
    n_iv <- nrow(harm)

    methods_to_run <- if (n_iv == 1) "mr_wald_ratio" else
                     if (n_iv == 2) "mr_ivw" else
                     c("mr_ivw", "mr_weighted_median", "mr_egger_regression",
                       "mr_weighted_mode")

    for (method_name in methods_to_run) {
      mr_res <- tryCatch(
        get(method_name)(harm$beta.exposure, harm$beta.outcome,
                          harm$se.exposure, harm$se.outcome),
        error = function(e) NULL
      )
      if (is.null(mr_res)) next
      results[[length(results) + 1]] <- data.table(
        panel = panel, biomarker = biomarker, method = method_name,
        n_iv = n_iv,
        b = mr_res$b %||% NA_real_,
        se = mr_res$se %||% NA_real_,
        pval = mr_res$pval %||% NA_real_
      )
    }

    # Heterogeneity + Egger intercept
    if (n_iv >= 3) {
      het <- tryCatch(
        mr_heterogeneity(harm, method_list = "mr_ivw"),
        error = function(e) NULL
      )
      pleio <- tryCatch(mr_pleiotropy_test(harm), error = function(e) NULL)
      if (!is.null(het) && nrow(het) > 0) {
        results[[length(results) + 1]] <- data.table(
          panel = panel, biomarker = biomarker, method = "heterogeneity_Q",
          n_iv = n_iv, b = het$Q[1], se = NA_real_, pval = het$Q_pval[1]
        )
      }
      if (!is.null(pleio) && nrow(pleio) > 0) {
        results[[length(results) + 1]] <- data.table(
          panel = panel, biomarker = biomarker, method = "egger_intercept",
          n_iv = n_iv,
          b = pleio$egger_intercept[1], se = pleio$se[1], pval = pleio$pval[1]
        )
      }
    }
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (length(results) == 0) {
  stop("No MR results produced. Check unified_iv_list + Timsina file paths.")
}

results_dt <- rbindlist(results, fill = TRUE)

out_fp <- file.path(dir_processed,
                    sprintf("csf_cascade_mr_results_%s.parquet",
                             format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(results_dt, out_fp)
write_parquet(results_dt, file.path(dir_processed, "csf_cascade_mr_results.parquet"))
cat(sprintf("[%s] Wrote %s (%d rows)\n", Sys.time(), basename(out_fp), nrow(results_dt)))

# --- Summary table for manuscript Table -------------------------------------

cat("\n=== CSF cascade MR summary (IVW estimates) ===\n")
ivw_summary <- dcast(results_dt[method == "mr_ivw"],
                     panel ~ biomarker, value.var = c("b", "se", "pval"))
print(ivw_summary)

cat("\n=== Per-panel concordance with biological expectation ===\n")
cat("  Expected directions: amyloid → CSF Aβ42 NEGATIVE; → t-tau POSITIVE; → p-tau181 POSITIVE\n\n")
for (pn in PANELS) {
  r <- results_dt[panel == pn & method == "ivw"]
  if (nrow(r) == 0) next
  cat(sprintf("Panel: %s\n", pn))
  for (bm in names(TIMSINA_FILES)) {
    row <- r[biomarker == bm]
    if (nrow(row) == 0) next
    expected_sign <- if (bm == "csf_abeta42") "-" else "+"
    actual_sign <- if (is.na(row$b)) "?" else if (row$b > 0) "+" else if (row$b < 0) "-" else "0"
    concord <- if (actual_sign == expected_sign) "✓" else if (actual_sign == "0") "—" else "✗"
    cat(sprintf("  %-15s β=%+.4f, p=%.2e  [expected %s, observed %s] %s\n",
                bm, row$b, row$pval, expected_sign, actual_sign, concord))
  }
  cat("\n")
}
