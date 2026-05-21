# =============================================================================
# Script 15: Amendment 4 Component B — Kim Korean as secondary primary discovery
# =============================================================================
#
# Re-runs the per-phecode locus-MR pipeline using Kim 2025 Korean (GCST90483382,
# N=3,387 EAS) exposure effects on the same 16 unified instrument variants,
# against the same FinnGen R12 outcome slices.
#
# This is the "wide-net" Component B implementation: we use the same instrument
# set the EUR primary uses, but plug in Kim Korean effect estimates as the
# exposure side. This is conceptually the cleanest test of:
#   "If we instrument amyloid PET with the same SNPs but read their effects
#    off the East Asian cohort instead of the European one, do we still see
#    the same downstream phecode effects?"
# Concordance between Ali EUR locus-MR (Phase 2, locked primary) and this Kim
# Korean re-run is methodologically informative triangulation.
#
# Note: a stricter Component B implementation would extract new instruments
# from Kim Korean directly (clump, threshold, etc.) and run a fresh FinnGen
# tabix sweep against those. That requires substantial additional infrastructure
# and is deferred to a future amendment if the wide-net result is interesting.
#
# Pre-reg: Amendment 4 Component B (osf/amendment_log.md, pushed to OSF 2026-05-20)
#
# Input:
#   - data/reference/kim_2025/GCST90483382.tsv (Kim Korean, N=3,387, GRCh38)
#   - data/interim/instruments/unified_iv_list.tsv (16 unified IVs)
#   - <DRIVE>/_shared_mirror/finngen_R12/sliced/*.tsv.gz (existing slices)
#
# Output:
#   - data/processed/kim_korean_locus_mr_results.parquet
#   - data/processed/kim_korean_iv_lookup.parquet (per-IV Kim Korean effects)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(arrow)
  library(R.utils)
})

proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)

KIM_FP <- file.path(proj_root, "data/reference/kim_2025/GCST90483382.tsv")
DRIVE_MIRROR <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12/sliced")
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
N_KIM <- 3387

SEED <- 81
set.seed(SEED)

# --- Inputs ------------------------------------------------------------------

iv_unified <- fread(file.path(dir_interim, "instruments/unified_iv_list.tsv"))
phecode_meta <- fread(file.path(proj_root, "osf/phecode_filter_prespec.csv"))
PANELS <- c("tier1_apoe_inclusive", "tier1_apoe_narrow_excl",
            "tier1_apoe_wide_excl", "tier2_apoe_narrow_excl",
            "tier2_apoe_wide_excl")

cat(sprintf("[%s] Loaded %d unified IVs; targeting %d phecodes × %d panels\n",
            Sys.time(), nrow(iv_unified), nrow(phecode_meta), length(PANELS)))

# --- Look up our 16 IVs in Kim Korean ----------------------------------------

cat(sprintf("[%s] Loading Kim 2025 Korean sumstats (~315 MB)\n", Sys.time()))
kim <- fread(KIM_FP)
setnames(kim, "rs_id", "rsid")
cat(sprintf("  Loaded: %d SNPs × %d cols\n", nrow(kim), ncol(kim)))

# Match by rsID
kim_ivs <- merge(iv_unified[, .(rsid, chr, hg38_pos, ea_eur = effect_allele,
                                 oa_eur = other_allele, beta_eur = beta,
                                 se_eur = se, panels)],
                 kim[, .(rsid, chr_kim = chromosome, pos_kim = base_pair_location,
                         ea_kim = effect_allele, oa_kim = other_allele,
                         beta_kim = beta, se_kim = standard_error,
                         eaf_kim = effect_allele_frequency, p_kim = p_value)],
                 by = "rsid")
cat(sprintf("[%s] IV lookup in Kim Korean: %d / %d matched\n",
            Sys.time(), nrow(kim_ivs), nrow(iv_unified)))

# Harmonise: flip Kim beta if EA differs from EUR's EA
kim_ivs[, flip := FALSE]
kim_ivs[toupper(ea_eur) == toupper(oa_kim) & toupper(oa_eur) == toupper(ea_kim), flip := TRUE]
kim_ivs[flip == TRUE, beta_kim := -beta_kim]
# Drop full mismatches
kim_ivs <- kim_ivs[(toupper(ea_eur) == toupper(ea_kim) & toupper(oa_eur) == toupper(oa_kim)) |
                    (toupper(ea_eur) == toupper(oa_kim) & toupper(oa_eur) == toupper(ea_kim))]
cat(sprintf("[%s] Post-harmonise: %d / 16 IVs usable (dropped %d for allele mismatch)\n",
            Sys.time(), nrow(kim_ivs), nrow(iv_unified) - nrow(kim_ivs)))

write_parquet(kim_ivs, file.path(dir_processed, "kim_korean_iv_lookup.parquet"))
cat(sprintf("[%s] Wrote per-IV Kim Korean lookup parquet\n", Sys.time()))

# Per-IV concordance check vs EUR (informative diagnostic)
cat("\n=== Per-IV concordance (EUR β vs Kim Korean β, harmonised on EUR EA) ===\n")
print(kim_ivs[order(chr, hg38_pos),
              .(rsid, chr, beta_eur = round(beta_eur, 3),
                beta_kim = round(beta_kim, 3),
                same_direction = sign(beta_eur) == sign(beta_kim),
                p_kim = signif(p_kim, 3))])
n_concord <- sum(sign(kim_ivs$beta_eur) == sign(kim_ivs$beta_kim))
cat(sprintf("\nConcordant direction: %d / %d IVs\n", n_concord, nrow(kim_ivs)))

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

# Format Kim Korean as exposure data per panel
format_exposure_kim <- function(iv_subset_with_kim, panel) {
  data.frame(
    SNP                    = iv_subset_with_kim$rsid,
    beta.exposure          = iv_subset_with_kim$beta_kim,
    se.exposure            = iv_subset_with_kim$se_kim,
    pval.exposure          = iv_subset_with_kim$p_kim,
    effect_allele.exposure = iv_subset_with_kim$ea_eur,  # already harmonised to EUR EA
    other_allele.exposure  = iv_subset_with_kim$oa_eur,
    eaf.exposure           = iv_subset_with_kim$eaf_kim,
    exposure               = "amyloid_PET_Kim_EAS",
    id.exposure            = panel,
    stringsAsFactors       = FALSE
  )
}

run_mr_for_pair <- function(harm, n_iv) {
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

# Per-phecode runner using Kim Korean exposure
run_panel_sweep_kim <- function(phecode_id, label, category) {
  outcome_dat <- read_slice(phecode_id)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    return(data.table(phecode = phecode_id, panel = NA_character_, method = NA_character_,
                      n_iv = 0L, b = NA_real_, se = NA_real_, pval = NA_real_,
                      skipped_reason = "missing_or_empty_slice"))
  }
  rows <- list()
  for (panel in PANELS) {
    iv_panel <- kim_ivs[grepl(panel, panels, fixed = TRUE)]
    if (nrow(iv_panel) == 0) next
    exposure_dat <- format_exposure_kim(iv_panel, panel)
    harm <- tryCatch(TwoSampleMR::harmonise_data(exposure_dat, outcome_dat, action = 2),
                     error = function(e) NULL)
    if (is.null(harm) || nrow(harm) == 0) next
    harm <- harm[harm$mr_keep == TRUE, ]
    n_iv <- nrow(harm)
    if (n_iv == 0) next
    mr_res <- run_mr_for_pair(harm, n_iv)
    if (is.null(mr_res)) next
    for (mname in names(mr_res)) {
      r <- mr_res[[mname]]
      rows[[length(rows) + 1]] <- data.table(
        phecode = phecode_id, panel = panel, method = mname, n_iv = n_iv,
        b = if (!is.null(r$b)) r$b else NA_real_,
        se = if (!is.null(r$se)) r$se else NA_real_,
        pval = if (!is.null(r$pval)) r$pval else NA_real_,
        skipped_reason = NA_character_
      )
    }
  }
  if (length(rows) == 0) {
    return(data.table(phecode = phecode_id, panel = NA_character_, method = NA_character_,
                      n_iv = 0L, b = NA_real_, se = NA_real_, pval = NA_real_,
                      skipped_reason = "no_harm_overlap"))
  }
  rbindlist(rows, fill = TRUE)
}

# --- Driver ------------------------------------------------------------------

n_total <- nrow(phecode_meta)
cat(sprintf("\n[%s] Running Kim Korean locus-MR sweep over %d phecodes...\n",
            Sys.time(), n_total))

t0 <- Sys.time()
results <- vector("list", n_total)
for (i in seq_len(n_total)) {
  results[[i]] <- run_panel_sweep_kim(
    phecode_meta$phenocode[i], phecode_meta$phenotype[i], phecode_meta$category[i])
  if (i %% 100 == 0 || i == n_total) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    rate <- i / elapsed
    eta_min <- (n_total - i) / rate / 60
    cat(sprintf("[%s] %d/%d (%.1f%%) — %.2f sec/phecode, ETA %.0f min\n",
                format(Sys.time(), "%H:%M:%S"),
                i, n_total, 100*i/n_total, 1/rate, eta_min))
  }
}

results_all <- rbindlist(results, fill = TRUE)
results_all <- merge(results_all,
                     phecode_meta[, .(phecode = phenocode, phenotype, category,
                                       num_cases, num_controls)],
                     by = "phecode", all.x = TRUE, sort = FALSE)

out_fp <- file.path(dir_processed,
                    sprintf("kim_korean_locus_mr_results_%s.parquet",
                             format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(results_all, out_fp)
write_parquet(results_all,
              file.path(dir_processed, "kim_korean_locus_mr_results.parquet"))
cat(sprintf("[%s] Wrote %s (%d rows)\n", Sys.time(), basename(out_fp), nrow(results_all)))

# Summary by panel × method
cat("\n=== Kim Korean MR summary (n with method assigned) ===\n")
print(results_all[!is.na(method), .N, by = .(panel, method)])
