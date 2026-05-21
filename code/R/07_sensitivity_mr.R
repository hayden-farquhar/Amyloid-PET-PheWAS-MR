# =============================================================================
# Script 07: Phase 3 — Sensitivity MR (MR-PRESSO + MR-RAPS)
# =============================================================================
#
# Runs bootstrap-heavy sensitivity MR methods on the subset of phecodes that
# look interesting in the Phase 2 locus-level sweep. Specifically:
#
#   MR-PRESSO (Verbanck 2018) — global pleiotropy test + outlier-corrected IVW.
#                                Requires n_IV >= 4. NbDistribution = 1000.
#   MR-RAPS   (Zhao 2020)     — Robust Adjusted Profile Score; primary
#                                estimator under Tier 2 relaxed panels.
#
# Pre-registration: §6.1 (MR-RAPS as Tier 2 primary), §13.3 (MR-PRESSO global +
# outlier test as primary pleiotropy diagnostic).
#
# Filter criteria for "interesting" — pre-registered §13.4 retraction triggers
# downstream + standard practice screen:
#   - Any panel IVW p < 1e-4 (~ Bonferroni 0.05 / 5*1574 ≈ 6e-6 would be too
#     strict for a sensitivity-suite screen; 1e-4 catches the tail liberally),
#   AND
#   - n_IV >= 4 in at least one Tier 2 panel (MR-PRESSO requires this).
#
# Input:
#   - data/processed/locus_mr_results.parquet (or its checkpoint while Phase 2
#     is still running)
#   - data/interim/instruments/unified_iv_list.tsv
#   - <DRIVE>/_shared_mirror/finngen_R12/sliced/<phecode>.tsv.gz
#
# Output:
#   - data/processed/sensitivity_mr_results.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(arrow)
  library(R.utils)
})

# --- Project paths -----------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
DRIVE_MIRROR <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12/sliced")
dir_processed <- file.path(proj_root, "data", "processed")
SEED <- 81
set.seed(SEED)

# --- Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
p_threshold <- if (length(args) >= 1) as.numeric(args[1]) else 1e-4
panels_to_run <- if (length(args) >= 2) {
  strsplit(args[2], ",")[[1]]
} else {
  c("tier2_apoe_narrow_excl", "tier2_apoe_wide_excl", "tier1_apoe_wide_excl")
}
n_bootstrap <- if (length(args) >= 3) as.integer(args[3]) else 1000L

cat(sprintf("[%s] Sensitivity MR: filter p < %.0e, panels = [%s], n_boot = %d\n",
            Sys.time(), p_threshold, paste(panels_to_run, collapse=", "),
            n_bootstrap))

# --- Identify interesting phecodes from Phase 2 ------------------------------

# Try main parquet first; fall back to checkpoint if Phase 2 still running
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

phase2_fp <- file.path(dir_processed, "locus_mr_results.parquet")
phase2_ck <- file.path(dir_processed, "locus_mr_checkpoint.parquet")
ck_size <- if (file.exists(phase2_ck)) file.size(phase2_ck) else 0
if (file.exists(phase2_fp) && file.size(phase2_fp) > ck_size) {
  cat(sprintf("[%s] Reading Phase 2 results from %s\n", Sys.time(), phase2_fp))
  phase2 <- as.data.table(read_parquet(phase2_fp))
} else if (file.exists(phase2_ck)) {
  cat(sprintf("[%s] Reading Phase 2 CHECKPOINT (Phase 2 still running)\n",
              Sys.time()))
  phase2 <- as.data.table(read_parquet(phase2_ck))
} else {
  stop("Neither Phase 2 results nor checkpoint exist. Run R/05_phase2_locus_mr.R first.")
}
cat(sprintf("  Phase 2 rows: %d, phecodes: %d\n",
            nrow(phase2), uniqueN(phase2$phecode)))

# Pull the interesting subset: any panel-method with p < threshold AND n_iv >= 4
interesting <- phase2[!is.na(pval) & pval < p_threshold &
                       panel %in% panels_to_run &
                       method == "ivw" & n_iv >= 4,
                      unique(phecode)]
cat(sprintf("[%s] Selected %d interesting phecodes for sensitivity sweep\n",
            Sys.time(), length(interesting)))

if (length(interesting) == 0) {
  cat("No interesting phecodes — exiting. Re-run with --p_threshold relaxed if needed.\n")
  quit(status = 0)
}

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(proj_root, "data/interim/instruments/unified_iv_list.tsv"))
phecode_meta <- fread(file.path(proj_root, "osf/phecode_filter_prespec.csv"))

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

read_outcome_slice <- function(phecode_id, label) {
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
    SNP                   = rsid_long$rsid,
    beta.outcome          = rsid_long$beta,
    se.outcome            = rsid_long$sebeta,
    pval.outcome          = rsid_long$pval,
    effect_allele.outcome = rsid_long$alt,
    other_allele.outcome  = rsid_long$ref,
    eaf.outcome           = rsid_long$af_alt,
    outcome               = label,
    id.outcome            = phecode_id,
    stringsAsFactors      = FALSE
  )
}

# --- Per-phecode × per-panel sensitivity runner -----------------------------

run_sensitivity_one <- function(phecode_id, label, panel) {
  iv_sub <- iv_unified[grepl(panel, panels, fixed = TRUE)]
  if (nrow(iv_sub) < 4) {
    return(data.table(phecode = phecode_id, panel = panel, method = "skip",
                      n_iv = nrow(iv_sub),
                      reason = "n_iv < 4 (MR-PRESSO requires >= 4)"))
  }
  outcome_dat <- read_outcome_slice(phecode_id, label)
  if (is.null(outcome_dat)) {
    return(data.table(phecode = phecode_id, panel = panel, method = "skip",
                      n_iv = NA_integer_, reason = "missing_or_empty_slice"))
  }
  harm <- tryCatch(
    TwoSampleMR::harmonise_data(format_exposure(iv_sub, panel), outcome_dat,
                                action = 2),
    error = function(e) NULL
  )
  if (is.null(harm) || sum(harm$mr_keep) < 4) {
    return(data.table(phecode = phecode_id, panel = panel, method = "skip",
                      n_iv = sum(harm$mr_keep %||% 0),
                      reason = "post-harm n_iv < 4"))
  }
  harm <- harm[harm$mr_keep == TRUE, ]
  n_iv <- nrow(harm)

  rows <- list()

  # --- MR-PRESSO -----------------------------------------------------------
  presso <- tryCatch({
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = harm, NbDistribution = n_bootstrap, SignifThreshold = 0.05
    )
  }, error = function(e) NULL)

  if (!is.null(presso)) {
    raw <- presso$`Main MR results`
    raw_idx <- which(raw$`MR Analysis` == "Raw")
    out_idx <- which(raw$`MR Analysis` == "Outlier-corrected")
    rows[[length(rows) + 1]] <- data.table(
      phecode = phecode_id, panel = panel, method = "presso_raw",
      n_iv = n_iv, b = raw$`Causal Estimate`[raw_idx],
      se = raw$Sd[raw_idx], pval = raw$`P-value`[raw_idx],
      global_p = presso$`MR-PRESSO results`$`Global Test`$Pvalue,
      n_outliers = if (length(out_idx) > 0)
        length(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
        else 0L,
      reason = NA_character_
    )
    if (length(out_idx) > 0 && !is.na(raw$`Causal Estimate`[out_idx])) {
      rows[[length(rows) + 1]] <- data.table(
        phecode = phecode_id, panel = panel, method = "presso_corrected",
        n_iv = n_iv - length(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
        b = raw$`Causal Estimate`[out_idx], se = raw$Sd[out_idx],
        pval = raw$`P-value`[out_idx],
        global_p = presso$`MR-PRESSO results`$`Global Test`$Pvalue,
        n_outliers = length(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
        reason = NA_character_
      )
    }
  }

  # --- MR-RAPS -------------------------------------------------------------
  raps <- tryCatch({
    TwoSampleMR::mr_raps(harm$beta.exposure, harm$beta.outcome,
                         harm$se.exposure, harm$se.outcome)
  }, error = function(e) NULL)
  if (!is.null(raps)) {
    rows[[length(rows) + 1]] <- data.table(
      phecode = phecode_id, panel = panel, method = "mr_raps",
      n_iv = n_iv, b = raps$b, se = raps$se, pval = raps$pval,
      global_p = NA_real_, n_outliers = NA_integer_,
      reason = NA_character_
    )
  }

  rbindlist(rows, fill = TRUE)
}

# --- Driver ------------------------------------------------------------------

n_jobs <- length(interesting) * length(panels_to_run)
cat(sprintf("[%s] Running %d phecodes × %d panels = %d sensitivity jobs\n",
            Sys.time(), length(interesting), length(panels_to_run), n_jobs))

t0 <- Sys.time()
all_rows <- list()
job_i <- 0
for (phecode_id in interesting) {
  label <- phecode_meta[phenocode == phecode_id, phenotype][1]
  for (panel in panels_to_run) {
    job_i <- job_i + 1
    rows <- run_sensitivity_one(phecode_id, label, panel)
    all_rows[[length(all_rows) + 1]] <- rows
    if (job_i %% 5 == 0 || job_i == n_jobs) {
      elapsed <- as.numeric(Sys.time() - t0, units = "secs")
      rate <- job_i / elapsed
      eta_min <- (n_jobs - job_i) / rate / 60
      cat(sprintf("[%s] job %d/%d (%.1f%%) — %.1f sec/job, ETA %.1f min\n",
                  format(Sys.time(), "%H:%M:%S"),
                  job_i, n_jobs, 100 * job_i / n_jobs,
                  1/rate, eta_min))
    }
  }
}

results_all <- rbindlist(all_rows, fill = TRUE)
results_all <- merge(results_all,
                     phecode_meta[, .(phecode = phenocode, phenotype,
                                       category, num_cases, num_controls)],
                     by = "phecode", all.x = TRUE, sort = FALSE)

out_fp <- file.path(dir_processed,
                    sprintf("sensitivity_mr_results_%s.parquet",
                             format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(results_all, out_fp)
write_parquet(results_all,
              file.path(dir_processed, "sensitivity_mr_results.parquet"))
cat(sprintf("[%s] Wrote %s (%d rows)\n",
            Sys.time(), out_fp, nrow(results_all)))

# Summary
cat("\nResults by panel × method:\n")
print(results_all[!is.na(b), .(.N, n_iv_med = median(n_iv, na.rm = TRUE),
                                global_p_min = min(global_p, na.rm = TRUE)),
                  by = .(panel, method)])
