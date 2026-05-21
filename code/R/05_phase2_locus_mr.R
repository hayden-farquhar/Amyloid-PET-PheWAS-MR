# =============================================================================
# Script 05: Phase 2 — Per-phecode locus-level MR sweep
# =============================================================================
#
# For each of the 1,574 FinnGen R12 phecodes and each of the 5 instrument
# panels (Tier 1 inclusive/narrow/wide-excl + Tier 2 narrow/wide-excl), runs:
#   - IVW (random-effects)
#   - Weighted median
#   - MR-Egger (with intercept p-value)
#   - Weighted mode
#   - Cochran's Q heterogeneity test
#
# Sensitivity methods (MR-PRESSO, MR-RAPS) are deferred to script 06.
#
# Per pre-reg §6.1: Tier 1 wide-excluded is the conservative non-APOE primary
# under Pivot A activation. Tier 2 panels run with MR-RAPS as primary
# locus-level estimator. Tier 1 APOE-inclusive runs descriptively (APOE-
# dominated; cannot disentangle amyloid effect from APOE per se).
#
# Input:
#   - data/interim/instruments/unified_iv_list.tsv (16 rsIDs × panel membership)
#   - <DRIVE>/_shared_mirror/finngen_R12/sliced/<phecode>.tsv.gz
#   - osf/phecode_filter_prespec.csv (1,574 phecode metadata)
#
# Output:
#   - data/processed/locus_mr_results.parquet (one row per phecode × panel × method)
#
# Pre-registration: Sections 9, 12, 13 (primary MR); Section 27 (Pivot A)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(arrow)
  library(R.utils)
})

# --- Project paths -----------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)

DRIVE_MIRROR <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12/sliced")

dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)

SEED <- 81
set.seed(SEED)

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(dir_interim, "instruments", "unified_iv_list.tsv"))
phecode_meta <- fread(file.path(proj_root, "osf", "phecode_filter_prespec.csv"))

PANELS <- c("tier1_apoe_inclusive",
            "tier1_apoe_narrow_excl",
            "tier1_apoe_wide_excl",
            "tier2_apoe_narrow_excl",
            "tier2_apoe_wide_excl")

cat(sprintf("[%s] Loaded %d unified IVs across %d phecodes; %d panels\n",
            Sys.time(), nrow(iv_unified), nrow(phecode_meta), length(PANELS)))

# --- Helpers -----------------------------------------------------------------

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

# Read a single FinnGen slice and reformat as TwoSampleMR outcome.
# Returns NULL if slice missing / empty.
read_outcome_slice <- function(phecode_id, label, category) {
  fp <- file.path(DRIVE_MIRROR, paste0(phecode_id, ".tsv.gz"))
  if (!file.exists(fp)) return(NULL)
  sl <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(sl) || nrow(sl) == 0) return(NULL)
  setnames(sl, "#chrom", "chrom")
  # FinnGen multi-rsID rows have comma-separated rsids; we accept on any match.
  # Expand: one row per individual rsID for join against iv_unified$rsid.
  sl[, row_idx := .I]
  rsid_long <- sl[, .(rsid = unlist(strsplit(rsids, ","))), by = row_idx]
  rsid_long <- merge(rsid_long, sl[, .(row_idx, chrom, pos, ref, alt, beta, sebeta,
                                       pval, af_alt)],
                     by = "row_idx")
  # If FinnGen reports beta as NA, drop
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

run_mr_for_pair <- function(harm, n_iv) {
  # Return a list of method × stat rows. NULL if MR can't run.
  res <- list()
  if (n_iv == 0) return(NULL)
  if (n_iv == 1) {
    wr <- tryCatch(mr_wald_ratio(harm$beta.exposure, harm$beta.outcome,
                                  harm$se.exposure, harm$se.outcome),
                   error = function(e) NULL)
    if (!is.null(wr)) res[["wald_ratio"]] <- wr
    return(res)
  }
  if (n_iv >= 2) {
    ivw <- tryCatch(mr_ivw(harm$beta.exposure, harm$beta.outcome,
                            harm$se.exposure, harm$se.outcome),
                    error = function(e) NULL)
    if (!is.null(ivw)) res[["ivw"]] <- ivw
  }
  if (n_iv >= 3) {
    wm <- tryCatch(mr_weighted_median(harm$beta.exposure, harm$beta.outcome,
                                       harm$se.exposure, harm$se.outcome),
                   error = function(e) NULL)
    if (!is.null(wm)) res[["weighted_median"]] <- wm

    egger <- tryCatch(mr_egger_regression(harm$beta.exposure, harm$beta.outcome,
                                            harm$se.exposure, harm$se.outcome),
                       error = function(e) NULL)
    if (!is.null(egger)) res[["egger"]] <- egger

    wmode <- tryCatch(mr_weighted_mode(harm$beta.exposure, harm$beta.outcome,
                                         harm$se.exposure, harm$se.outcome),
                       error = function(e) NULL)
    if (!is.null(wmode)) res[["weighted_mode"]] <- wmode
  }
  res
}

# Run all 5 panels for one phecode. Returns a data.table of result rows.
run_panel_sweep <- function(phecode_id, label, category) {
  outcome_dat <- read_outcome_slice(phecode_id, label, category)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    return(data.table(phecode = phecode_id, panel = NA_character_,
                      method = NA_character_, n_iv = 0,
                      b = NA_real_, se = NA_real_, pval = NA_real_,
                      q_stat = NA_real_, q_pval = NA_real_,
                      egger_intercept = NA_real_, egger_intercept_pval = NA_real_,
                      skipped_reason = "missing_or_empty_slice"))
  }

  rows <- list()
  for (panel in PANELS) {
    iv_subset <- iv_unified[grepl(panel, panels, fixed = TRUE)]
    if (nrow(iv_subset) == 0) next
    exposure_dat <- format_exposure(iv_subset, panel)
    harm <- tryCatch(
      TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2),
      error = function(e) NULL
    )
    if (is.null(harm) || nrow(harm) == 0) next
    harm <- harm[harm$mr_keep == TRUE, ]
    n_iv <- nrow(harm)
    if (n_iv == 0) next

    # Heterogeneity (Q-stat from IVW) — only meaningful for n>=2
    het_q <- het_p <- NA_real_
    if (n_iv >= 2) {
      het <- tryCatch(
        mr_heterogeneity(data.frame(
          id.exposure = panel, id.outcome = phecode_id,
          exposure = "amyloid_PET", outcome = label,
          mr_keep = TRUE,
          beta.exposure = harm$beta.exposure, beta.outcome = harm$beta.outcome,
          se.exposure = harm$se.exposure, se.outcome = harm$se.outcome
        ), method_list = "mr_ivw"),
        error = function(e) NULL
      )
      if (!is.null(het) && nrow(het) > 0) {
        het_q <- het$Q[1]; het_p <- het$Q_pval[1]
      }
    }

    # Egger intercept p
    egger_int <- egger_intp <- NA_real_
    if (n_iv >= 3) {
      pi <- tryCatch(
        mr_pleiotropy_test(data.frame(
          id.exposure = panel, id.outcome = phecode_id,
          exposure = "amyloid_PET", outcome = label,
          mr_keep = TRUE,
          beta.exposure = harm$beta.exposure, beta.outcome = harm$beta.outcome,
          se.exposure = harm$se.exposure, se.outcome = harm$se.outcome
        )),
        error = function(e) NULL
      )
      if (!is.null(pi) && nrow(pi) > 0) {
        egger_int <- pi$egger_intercept[1]
        egger_intp <- pi$pval[1]
      }
    }

    mr_res <- run_mr_for_pair(harm, n_iv)
    if (is.null(mr_res)) next
    for (method_name in names(mr_res)) {
      r <- mr_res[[method_name]]
      rows[[length(rows) + 1]] <- data.table(
        phecode = phecode_id,
        panel = panel,
        method = method_name,
        n_iv = n_iv,
        b = if (!is.null(r$b)) r$b else NA_real_,
        se = if (!is.null(r$se)) r$se else NA_real_,
        pval = if (!is.null(r$pval)) r$pval else NA_real_,
        q_stat = het_q,
        q_pval = het_p,
        egger_intercept = egger_int,
        egger_intercept_pval = egger_intp,
        skipped_reason = NA_character_
      )
    }
  }
  if (length(rows) == 0) {
    return(data.table(phecode = phecode_id, panel = NA_character_,
                      method = NA_character_, n_iv = 0,
                      b = NA_real_, se = NA_real_, pval = NA_real_,
                      q_stat = NA_real_, q_pval = NA_real_,
                      egger_intercept = NA_real_, egger_intercept_pval = NA_real_,
                      skipped_reason = "no_harm_overlap"))
  }
  rbindlist(rows, fill = TRUE)
}

# --- Driver ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "dry-run"
workers <- if (length(args) >= 2) as.integer(args[2]) else 4

# Build the work table. If a slice file is missing the phecode is still listed
# (will return skipped_reason).
work_tbl <- phecode_meta[, .(phenocode, phenotype, category)]

if (mode == "dry-run") work_tbl <- head(work_tbl, 10)

n_total <- nrow(work_tbl)
cat(sprintf("[%s] %s mode: %d phecodes, sequential with progress\n",
            Sys.time(), mode, n_total))

# Resume support: skip phecodes that already exist in a checkpoint
checkpoint_fp <- file.path(dir_processed, "locus_mr_checkpoint.parquet")
already_done <- character(0)
if (file.exists(checkpoint_fp)) {
  ck <- as.data.table(read_parquet(checkpoint_fp))
  already_done <- unique(ck$phecode)
  cat(sprintf("[%s] Resuming: %d phecodes already in checkpoint\n",
              Sys.time(), length(already_done)))
}

results <- vector("list", n_total)
checkpoint_interval <- 100
t0 <- Sys.time()
for (i in seq_len(n_total)) {
  pid <- work_tbl$phenocode[i]
  if (pid %in% already_done) {
    results[[i]] <- ck[phecode == pid]
    next
  }
  results[[i]] <- run_panel_sweep(pid, work_tbl$phenotype[i], work_tbl$category[i])

  if (i %% 50 == 0 || i == n_total) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    rate <- i / elapsed
    eta_sec <- (n_total - i) / rate
    cat(sprintf("[%s] %d/%d (%.1f%%) — %.2f sec/phecode, ETA %.0f min\n",
                format(Sys.time(), "%H:%M:%S"),
                i, n_total, 100 * i / n_total, 1/rate, eta_sec / 60))
  }

  if (i %% checkpoint_interval == 0 || i == n_total) {
    partial <- rbindlist(results[seq_len(i)], fill = TRUE)
    write_parquet(partial, checkpoint_fp)
  }
}

results_all <- rbindlist(results, fill = TRUE)
# Attach metadata
results_all <- merge(results_all, phecode_meta[, .(phecode = phenocode,
                                                    phenotype, category,
                                                    num_cases, num_controls)],
                     by = "phecode", all.x = TRUE, sort = FALSE)

out_fp <- file.path(dir_processed,
                    paste0("locus_mr_results_", format(Sys.time(), "%Y%m%d-%H%M%S"),
                           ".parquet"))
write_parquet(results_all, out_fp)

# Also write a "current" alias
write_parquet(results_all, file.path(dir_processed, "locus_mr_results.parquet"))

cat(sprintf("[%s] Wrote %s (%d rows)\n",
            Sys.time(), out_fp, nrow(results_all)))
cat("Summary by panel × method:\n")
print(results_all[!is.na(method),
                  .(.N, n_iv_med = median(n_iv, na.rm = TRUE)),
                  by = .(panel, method)])
