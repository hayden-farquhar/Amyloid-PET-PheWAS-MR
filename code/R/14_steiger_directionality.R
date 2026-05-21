# =============================================================================
# Script 14: Steiger directionality test (reverse-causation guard)
# =============================================================================
#
# Per pre-reg §10.3: tests whether each instrument-outcome pair satisfies the
# Steiger directionality assumption that the SNP explains MORE variance in the
# exposure (amyloid-PET) than in the outcome (each FinnGen phecode). When a
# SNP explains more variance in the outcome than the exposure, the MR causal
# direction may be inverted (reverse causation or pleiotropy).
#
# Methods:
#   - Exposure R² (continuous, Z-scored amyloid-PET): Burgess & Thompson 2017
#       R² = β² / (β² + n × SE²)
#   - Outcome R² (binary phecode): Lee et al. 2012 prevalence-adjusted formula
#       (handled internally by TwoSampleMR::directionality_test, which also
#       computes the Steiger Z and p-value).
#
# Output:
#   - data/processed/steiger_results.parquet — one row per (phecode × panel × IV)
#   - data/processed/steiger_summary.parquet — per (phecode × panel) summary:
#       fraction of IVs passing Steiger, mean Steiger p-value, count of
#       directionally-incorrect IVs.
#
# Integrates into the Causal Atlas via R/09: any phecode×panel with >50% of IVs
# failing Steiger gets flagged `steiger_concern = TRUE`.
#
# Run as:
#   Rscript R/14_steiger_directionality.R [p_threshold]
# where p_threshold (default 1e-4) restricts to Phase 2 hits worth testing.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(arrow)
  library(R.utils)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
DRIVE_MIRROR <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12/sliced")
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")

# --- Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
p_threshold <- if (length(args) >= 1) as.numeric(args[1]) else 1e-4

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(dir_interim, "instruments/unified_iv_list.tsv"))
phecode_meta <- fread(file.path(proj_root, "osf/phecode_filter_prespec.csv"))
N_AMYLOID <- 13409  # Ali 2023 multi-ethnic N

phase2_fp <- file.path(dir_processed, "locus_mr_results.parquet")
phase2_ck <- file.path(dir_processed, "locus_mr_checkpoint.parquet")
src <- if (file.exists(phase2_fp) && file.size(phase2_fp) > file.size(phase2_ck) %||% 0)
       phase2_fp else phase2_ck
phase2 <- as.data.table(read_parquet(src))
cat(sprintf("[%s] Phase 2 source: %s (%d rows)\n",
            Sys.time(), basename(src), nrow(phase2)))

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# --- Identify phecodes to test -----------------------------------------------

target_panels <- c("tier1_apoe_wide_excl", "tier2_apoe_narrow_excl",
                    "tier2_apoe_wide_excl", "tier1_apoe_narrow_excl",
                    "tier1_apoe_inclusive")
candidates <- phase2[!is.na(pval) & pval < p_threshold &
                       panel %in% target_panels & method == "ivw",
                     unique(phecode)]
cat(sprintf("[%s] Phecodes with IVW p < %.0e in any panel: %d\n",
            Sys.time(), p_threshold, length(candidates)))
if (length(candidates) == 0) {
  cat("No candidate phecodes — exiting. Re-run after Phase 2 completes or relax threshold.\n")
  quit(status = 0)
}

# --- Slice reader (same as Phase 2) ------------------------------------------

read_slice <- function(phecode_id) {
  fp <- file.path(DRIVE_MIRROR, paste0(phecode_id, ".tsv.gz"))
  if (!file.exists(fp)) return(NULL)
  sl <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(sl) || nrow(sl) == 0) return(NULL)
  setnames(sl, "#chrom", "chrom")
  sl[, row_idx := .I]
  rsid_long <- sl[, .(rsid = unlist(strsplit(rsids, ","))), by = row_idx]
  rsid_long <- merge(rsid_long, sl[, .(row_idx, chrom, pos, ref, alt, beta, sebeta,
                                       pval, af_alt)], by = "row_idx")
  rsid_long <- rsid_long[!is.na(beta) & !is.na(sebeta) & sebeta > 0]
  if (nrow(rsid_long) == 0) return(NULL)
  data.frame(
    SNP                    = rsid_long$rsid,
    beta.outcome           = rsid_long$beta,
    se.outcome             = rsid_long$sebeta,
    pval.outcome           = rsid_long$pval,
    effect_allele.outcome  = rsid_long$alt,
    other_allele.outcome   = rsid_long$ref,
    eaf.outcome            = rsid_long$af_alt,
    outcome                = phecode_id,
    id.outcome             = phecode_id,
    stringsAsFactors       = FALSE
  )
}

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

# --- Driver per (phecode × panel) -------------------------------------------

cat(sprintf("[%s] Steiger sweep: %d phecodes × %d panels\n",
            Sys.time(), length(candidates), length(target_panels)))

all_iv_rows <- list()
all_summary_rows <- list()
t0 <- Sys.time()

for (i in seq_along(candidates)) {
  phecode_id <- candidates[i]
  outcome_dat <- read_slice(phecode_id)
  if (is.null(outcome_dat)) next
  outcome_meta <- phecode_meta[phenocode == phecode_id]
  n_cases <- outcome_meta$num_cases[1]
  n_controls <- outcome_meta$num_controls[1]
  n_outcome <- n_cases + n_controls
  prevalence <- n_cases / n_outcome

  for (panel in target_panels) {
    iv_sub <- iv_unified[grepl(panel, panels, fixed = TRUE)]
    if (nrow(iv_sub) == 0) next
    exposure_dat <- format_exposure(iv_sub, panel)
    harm <- tryCatch(
      TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2),
      error = function(e) NULL
    )
    if (is.null(harm) || nrow(harm) == 0) next
    harm <- harm[harm$mr_keep == TRUE, ]
    if (nrow(harm) == 0) next

    # TwoSampleMR's directionality_test expects per-SNP units/eaf
    # samplesize columns. Add them.
    harm$samplesize.exposure <- N_AMYLOID
    harm$samplesize.outcome  <- n_outcome
    harm$units.exposure <- "SD"  # amyloid-PET is Z-scored in Ali 2023
    harm$units.outcome  <- "log_odds"
    harm$ncase.outcome <- n_cases
    harm$ncontrol.outcome <- n_controls
    harm$prevalence.outcome <- prevalence

    steiger <- tryCatch(
      TwoSampleMR::directionality_test(harm),
      error = function(e) NULL
    )

    if (is.null(steiger) || nrow(steiger) == 0) next

    # Per-IV Steiger via TwoSampleMR::mr_steiger_filter (or compute manually)
    # The directionality_test gives aggregate (per phecode×panel).
    # For per-IV, compute manually:
    iv_steiger <- copy(harm)
    setDT(iv_steiger)
    iv_steiger[, r2_exposure := beta.exposure^2 /
                 (beta.exposure^2 + N_AMYLOID * se.exposure^2)]
    # Lee 2012 liability-scale R² for binary outcome
    # R²_liab = R²_observed × K(1-K) / (z² × (1-K-z²×K(1-K)))
    # Simplified: per-SNP outcome R² on observed scale, then transform
    z2 <- qnorm(prevalence)^2
    iv_steiger[, r2_outcome_obs := beta.outcome^2 /
                 (beta.outcome^2 + n_outcome * se.outcome^2)]
    iv_steiger[, r2_outcome := r2_outcome_obs *
                  prevalence * (1 - prevalence) /
                  (dnorm(qnorm(prevalence))^2)]
    iv_steiger[, steiger_pass := r2_exposure > r2_outcome]
    iv_steiger[, phecode := phecode_id]
    iv_steiger[, panel := panel]

    all_iv_rows[[length(all_iv_rows) + 1]] <- iv_steiger[, .(
      phecode, panel, SNP, beta.exposure, se.exposure, beta.outcome, se.outcome,
      r2_exposure, r2_outcome, steiger_pass
    )]

    # Per phecode×panel summary
    all_summary_rows[[length(all_summary_rows) + 1]] <- data.table(
      phecode = phecode_id, panel = panel,
      n_iv = nrow(iv_steiger),
      n_pass = sum(iv_steiger$steiger_pass, na.rm = TRUE),
      frac_pass = mean(iv_steiger$steiger_pass, na.rm = TRUE),
      aggregate_correct_direction = steiger$correct_causal_direction[1],
      aggregate_steiger_pval = steiger$steiger_pval[1]
    )
  }
  if (i %% 10 == 0 || i == length(candidates)) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    eta_min <- (length(candidates) - i) / (i/elapsed) / 60
    cat(sprintf("[%s] %d/%d (%.1f%%) — ETA %.1f min\n",
                format(Sys.time(), "%H:%M:%S"), i, length(candidates),
                100*i/length(candidates), eta_min))
  }
}

iv_results <- rbindlist(all_iv_rows, fill = TRUE)
summary_results <- rbindlist(all_summary_rows, fill = TRUE)

if (nrow(iv_results) == 0) {
  cat("No Steiger results produced.\n"); quit(status = 0)
}

# --- Write outputs ---------------------------------------------------------

out_iv <- file.path(dir_processed, "steiger_results.parquet")
out_sum <- file.path(dir_processed, "steiger_summary.parquet")
write_parquet(iv_results, out_iv)
write_parquet(summary_results, out_sum)
cat(sprintf("[%s] Wrote %s (%d rows) + %s (%d rows)\n",
            Sys.time(), basename(out_iv), nrow(iv_results),
            basename(out_sum), nrow(summary_results)))

# --- Summary report --------------------------------------------------------

cat("\n=== Steiger summary by panel ===\n")
print(summary_results[, .(
  n_phecodes_tested = .N,
  median_frac_pass = round(median(frac_pass, na.rm = TRUE), 3),
  n_concern_phecodes = sum(frac_pass < 0.5, na.rm = TRUE)
), by = panel])

cat("\n=== Phecodes with Steiger concern (>50% IVs failing) ===\n")
concern <- summary_results[frac_pass < 0.5][order(frac_pass)]
if (nrow(concern) > 0) {
  print(concern[1:min(20, .N), .(phecode, panel, n_iv, n_pass, frac_pass,
                                  aggregate_steiger_pval = signif(aggregate_steiger_pval, 3))])
} else {
  cat("  (no phecodes with majority IVs failing Steiger — clean directionality)\n")
}
