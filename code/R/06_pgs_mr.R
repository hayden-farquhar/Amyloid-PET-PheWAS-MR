# =============================================================================
# Script 06: PGS-MR — Bayesian-shrunk polygenic-score Mendelian randomization
# =============================================================================
#
# Under Pivot A activation (Tier 1 wide-excluded n_IV = 3 < 5 threshold), the
# PGS-MR arm becomes the primary inferential layer at the locus/genome level.
#
# Two scopes are supported:
#
#   1. NARROW (--scope narrow): use the existing 16 unified IVs with PRS-CS
#      posterior weights replacing the standard 1/SE² IVW weights. Conservative
#      sensitivity check on the locus-level IVW; runs immediately on existing
#      FinnGen slices. ~5 min wall-clock.
#
#   2. WIDE (--scope wide): use the top-N PRS-CS-weighted SNPs (default 500)
#      via a new FinnGen tabix sweep. Implements the textbook Burgess aggregate
#      PGS-MR; ~30-60 min additional sweep + ~5 min computation.
#
# Pre-registration: §6.10, §27 (Pivot A)
#
# Input:
#   - data/interim/prscs/amyloid_pgs_pst_eff_a1_b0.5_phi1e-02_chr{1..22}.txt
#   - data/interim/instruments/unified_iv_list.tsv (NARROW)
#   - <DRIVE>/_shared_mirror/finngen_R12/sliced/<phecode>.tsv.gz (NARROW)
#   - <DRIVE>/_shared_mirror/finngen_R12_wide/sliced/<phecode>.tsv.gz (WIDE; needs prep step)
#   - osf/phecode_filter_prespec.csv
#
# Output:
#   - data/processed/pgs_mr_results_<scope>.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(R.utils)
})

# --- Project paths -----------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)

DRIVE_MIRROR_NARROW <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12/sliced")
DRIVE_MIRROR_WIDE <- file.path(Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-hayden.l.farquhar@gmail.com",
  "My Drive/Research Projects/_shared_mirror/finngen_R12_wide/sliced")

dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)

SEED <- 81
set.seed(SEED)

# --- Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
scope <- if (length(args) >= 1) args[1] else "narrow"
top_n_wide <- if (length(args) >= 2) as.integer(args[2]) else 500L

stopifnot(scope %in% c("narrow", "wide"))

# --- Load PRS-CS posterior weights -------------------------------------------

prscs_fp_pattern <- file.path(dir_interim, "prscs",
                              "amyloid_pgs_pst_eff_a1_b0.5_phi1e-02_chr%d.txt")
chr_files <- sapply(1:22, function(chr) sprintf(prscs_fp_pattern, chr))
chr_avail <- file.exists(chr_files) & file.size(chr_files) > 0
cat(sprintf("[%s] PRS-CS chr files: %d/22 available\n",
            Sys.time(), sum(chr_avail)))
if (sum(chr_avail) < 22) {
  cat("  Missing chrs:", paste(which(!chr_avail), collapse = ", "), "\n")
  stop("PGS-MR requires all 22 PRS-CS chromosome outputs.")
}

prscs <- rbindlist(lapply(chr_files, function(fp) {
  fread(fp, col.names = c("chr", "rsid", "bp", "A1", "A2", "beta_pgs"),
        colClasses = c("integer", "character", "integer",
                       "character", "character", "numeric"))
}))
cat(sprintf("[%s] Loaded %d PRS-CS-weighted SNPs across 22 chromosomes\n",
            Sys.time(), nrow(prscs)))
cat(sprintf("  |beta_pgs| range: %.2e to %.2e (median %.2e)\n",
            min(abs(prscs$beta_pgs)), max(abs(prscs$beta_pgs)),
            median(abs(prscs$beta_pgs))))

# --- Compute aggregate PGS-MR helper -----------------------------------------
#
# Burgess aggregate MR with PRS-CS shrunken weights:
#
#   β_PGS-MR = Σ (w_j × β_outcome_j) / Σ (w_j × β_exposure_j)
#
# Variance via delta method on per-SNP IVs:
#   var(β_outcome / β_exposure) ≈ var(β_outcome) / β_exposure² +
#                                  β_outcome² × var(β_exposure) / β_exposure⁴
# Then weighted average with weights w_j² × β_exposure_j².
#
# We use β_pgs as both the PGS weight AND the exposure effect (it IS the
# shrunken exposure estimate). The aggregate then becomes:
#   β_PGS-MR = Σ (β_pgs_j × β_outcome_j) / Σ (β_pgs_j²)
#   var(β_PGS-MR) ≈ Σ (β_pgs_j² × se_outcome_j²) / (Σ β_pgs_j²)²
#
# This is the "summary-statistic projection" of Vilhjálmsson et al. 2015
# applied to MR-IVW.
pgs_mr_estimate <- function(beta_pgs, beta_out, se_out) {
  denom <- sum(beta_pgs^2)
  if (denom <= 0) return(list(b = NA_real_, se = NA_real_, n = 0L))
  numer <- sum(beta_pgs * beta_out)
  var_b <- sum((beta_pgs^2) * (se_out^2)) / denom^2
  list(b = numer / denom, se = sqrt(var_b), n = length(beta_pgs))
}

# --- Outcome slice reader ----------------------------------------------------

read_slice <- function(phecode_id, slice_dir) {
  fp <- file.path(slice_dir, paste0(phecode_id, ".tsv.gz"))
  if (!file.exists(fp)) return(NULL)
  sl <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(sl) || nrow(sl) == 0) return(NULL)
  setnames(sl, "#chrom", "chrom")
  sl[, row_idx := .I]
  rsid_long <- sl[, .(rsid = unlist(strsplit(rsids, ","))), by = row_idx]
  rsid_long <- merge(rsid_long, sl[, .(row_idx, chrom, pos, ref, alt, beta, sebeta,
                                       pval, af_alt)], by = "row_idx")
  rsid_long[!is.na(beta) & !is.na(sebeta) & sebeta > 0]
}

# --- Allele harmonisation: flip outcome beta if alleles swapped -------------
#
# PRS-CS A1 = effect allele, A2 = other allele.
# FinnGen alt = effect allele, ref = other allele.
# If (A1==alt, A2==ref): same direction.
# If (A1==ref, A2==alt): flip sign on beta_out.
# If alleles mismatch entirely: drop SNP.
harmonise_outcome <- function(prscs_row, outcome_row) {
  pa1 <- toupper(prscs_row$A1); pa2 <- toupper(prscs_row$A2)
  oea <- toupper(outcome_row$alt); ooa <- toupper(outcome_row$ref)
  if (pa1 == oea && pa2 == ooa) return(list(beta = outcome_row$beta, se = outcome_row$sebeta, keep = TRUE))
  if (pa1 == ooa && pa2 == oea) return(list(beta = -outcome_row$beta, se = outcome_row$sebeta, keep = TRUE))
  list(beta = NA_real_, se = NA_real_, keep = FALSE)
}

# --- Driver per phecode -----------------------------------------------------

run_pgs_mr_per_phecode <- function(phecode_id, prscs_sub, slice_dir) {
  outcome <- read_slice(phecode_id, slice_dir)
  if (is.null(outcome) || nrow(outcome) == 0) {
    return(data.table(phecode = phecode_id, scope = scope,
                      n_snps_overlap = 0L,
                      b = NA_real_, se = NA_real_, pval = NA_real_,
                      skipped_reason = "missing_or_empty_slice"))
  }
  m <- merge(prscs_sub[, .(rsid, A1, A2, beta_pgs)],
             outcome[, .(rsid, alt, ref, beta, sebeta)],
             by = "rsid")
  if (nrow(m) == 0) {
    return(data.table(phecode = phecode_id, scope = scope, n_snps_overlap = 0L,
                      b = NA_real_, se = NA_real_, pval = NA_real_,
                      skipped_reason = "no_snp_overlap_prscs_outcome"))
  }
  harm <- m[, {
    h <- harmonise_outcome(.SD, .SD)
    list(beta_out = h$beta, se_out = h$se, keep = h$keep)
  }, by = .(rsid, beta_pgs)]
  harm <- harm[keep == TRUE & !is.na(beta_out) & se_out > 0]
  if (nrow(harm) == 0) {
    return(data.table(phecode = phecode_id, scope = scope, n_snps_overlap = 0L,
                      b = NA_real_, se = NA_real_, pval = NA_real_,
                      skipped_reason = "all_dropped_allele_mismatch"))
  }
  est <- pgs_mr_estimate(harm$beta_pgs, harm$beta_out, harm$se_out)
  z <- est$b / est$se
  data.table(phecode = phecode_id, scope = scope,
             n_snps_overlap = est$n,
             b = est$b, se = est$se,
             pval = 2 * pnorm(-abs(z)),
             skipped_reason = NA_character_)
}

# --- Determine SNP subset and slice directory per scope ----------------------

phecode_meta <- fread(file.path(proj_root, "osf", "phecode_filter_prespec.csv"))

if (scope == "narrow") {
  iv_unified <- fread(file.path(dir_interim, "instruments", "unified_iv_list.tsv"))
  cat(sprintf("[%s] NARROW: 16 unified IVs intersected with PRS-CS weights\n",
              Sys.time()))
  prscs_sub <- prscs[rsid %in% iv_unified$rsid]
  cat(sprintf("  PRS-CS-weighted overlap: %d / 16 IVs\n", nrow(prscs_sub)))
  slice_dir <- DRIVE_MIRROR_NARROW
} else {
  # WIDE: use the pre-built non-APOE GRCh38-lifted IV list
  # (built by scripts/build_wide_pgs_noapoe.py to address the PRS-CS GRCh37 →
  # FinnGen GRCh38 build mismatch and the wide-APOE exclusion requirement).
  wide_iv_fp <- file.path(dir_interim, "instruments", "pgs_wide_iv_list_noAPOE.tsv")
  if (!file.exists(wide_iv_fp)) {
    cat("\n>>> WIDE scope requires the lifted non-APOE IV list at:\n>>> ", wide_iv_fp, "\n\n")
    cat("Build it with:\n")
    cat("  .venv/bin/python3 scripts/build_wide_pgs_noapoe.py 50\n\n")
    stop("WIDE IV list missing")
  }
  iv_wide <- fread(wide_iv_fp)
  cat(sprintf("[%s] WIDE (non-APOE, GRCh38-lifted): %d SNPs from %s\n",
              Sys.time(), nrow(iv_wide), basename(wide_iv_fp)))
  cat(sprintf("  Chromosomes: %s\n", paste(sort(unique(iv_wide$chr)), collapse = ",")))
  cat(sprintf("  Weight range: %.2e to %.2e\n",
              min(abs(iv_wide$beta_pgs)), max(abs(iv_wide$beta_pgs))))

  # Match the lifted IV list with the full PRS-CS posteriors by rsID to get
  # the canonical beta_pgs (this preserves the A1/A2 alignment from PRS-CS).
  prscs_sub <- merge(iv_wide[, .(rsid, chr_lifted = chr, hg38_pos)],
                     prscs[, .(rsid, A1, A2, beta_pgs)], by = "rsid")
  cat(sprintf("  PRS-CS×IV-list join: %d / %d\n", nrow(prscs_sub), nrow(iv_wide)))

  if (!dir.exists(DRIVE_MIRROR_WIDE)) {
    cat("\n>>> WIDE scope requires the FinnGen sweep to be run:\n")
    cat("  bash scripts/sweep_finngen_r12.sh \\\n")
    cat("    \"$DRIVE/Research Projects/_shared_mirror/finngen_R12_wide/\" \\\n")
    cat("    ", wide_iv_fp, " 4\n", sep = "")
    stop("WIDE scope sweep directory missing")
  }
  slice_dir <- DRIVE_MIRROR_WIDE
}

# --- Run per phecode ---------------------------------------------------------

cat(sprintf("[%s] Running PGS-MR (%s scope) over %d phecodes...\n",
            Sys.time(), scope, nrow(phecode_meta)))

t0 <- Sys.time()
results <- vector("list", nrow(phecode_meta))
for (i in seq_len(nrow(phecode_meta))) {
  results[[i]] <- run_pgs_mr_per_phecode(phecode_meta$phenocode[i], prscs_sub, slice_dir)
  if (i %% 100 == 0 || i == nrow(phecode_meta)) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    rate <- i / elapsed
    eta_min <- (nrow(phecode_meta) - i) / rate / 60
    cat(sprintf("[%s] %d/%d (%.1f%%) — %.2f sec/phecode, ETA %.1f min\n",
                format(Sys.time(), "%H:%M:%S"),
                i, nrow(phecode_meta), 100 * i / nrow(phecode_meta),
                1/rate, eta_min))
  }
}

results_all <- rbindlist(results, fill = TRUE)
results_all <- merge(results_all,
                     phecode_meta[, .(phecode = phenocode, phenotype,
                                       category, num_cases, num_controls)],
                     by = "phecode", all.x = TRUE, sort = FALSE)

out_fp <- file.path(dir_processed,
                    sprintf("pgs_mr_results_%s_%s.parquet",
                             scope, format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(results_all, out_fp)
write_parquet(results_all,
              file.path(dir_processed, sprintf("pgs_mr_results_%s.parquet", scope)))

cat(sprintf("[%s] Wrote %s (%d rows)\n",
            Sys.time(), out_fp, nrow(results_all)))

# Quick summary of hits (p < 1e-4, n_snps_overlap > 1)
hits <- results_all[!is.na(b) & pval < 1e-4 & n_snps_overlap >= 2]
cat(sprintf("\nPotential hits (p < 1e-4, n_snps >= 2): %d phecodes\n",
            nrow(hits)))
if (nrow(hits) > 0) {
  print(hits[order(pval),
             .(phecode, phenotype, category, n_snps_overlap,
               b = round(b, 4), se = round(se, 4),
               pval = signif(pval, 3))][1:min(20, .N)])
}
