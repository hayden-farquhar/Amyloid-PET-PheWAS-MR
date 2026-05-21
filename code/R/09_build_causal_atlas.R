# =============================================================================
# Script 09: Build the Amyloid-PET Causal Atlas (primary deliverable)
# =============================================================================
#
# Under Pivot A activation (per Phase 1), the Amyloid-PET Causal Atlas is the
# manuscript's primary deliverable (per pre-reg §27 and §38.5). This script
# assembles all evidence layers into one wide parquet:
#
#   - Locus-level MR (Phase 2): IVW, Weighted Median, MR-Egger, Weighted Mode
#                               + heterogeneity Q + Egger intercept (per panel)
#   - Sensitivity MR (Phase 3): MR-PRESSO raw/corrected, MR-RAPS (per panel)
#   - PGS-MR (Phase 6 narrow scope): PRS-CS-weighted aggregate ratio
#   - Colocalisation (Phase 4): best PP.H4 across brain tissues + genes
#
# Composite "robust hit" flag is computed per pre-reg §28: a hit must pass
# IVW + Weighted Median + Egger directional consistency + non-significant
# Egger intercept + non-significant Cochran's Q at the panel's primary level.
#
# Input parquets (any missing layer is tolerated — atlas is built incrementally):
#   - data/processed/locus_mr_results.parquet (Phase 2 — required)
#   - data/processed/sensitivity_mr_results.parquet (Phase 3 — optional)
#   - data/processed/pgs_mr_results_narrow.parquet (PGS-MR narrow — optional)
#   - data/processed/pgs_mr_results_wide.parquet (PGS-MR wide — optional)
#   - data/processed/coloc_results.parquet (Phase 4 — optional)
#
# Output:
#   - data/processed/amyloid_pet_causal_atlas.parquet (one row per phecode × panel,
#     ~7,870 rows across 5 panels × 1,574 phecodes)
#   - data/processed/causal_atlas_top_hits.parquet (filtered to "robust" + p signal)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_processed <- file.path(proj_root, "data", "processed")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)

# --- Args --------------------------------------------------------------------
# Allow override of robust-hit thresholds via CLI
args <- commandArgs(trailingOnly = TRUE)
robust_p_threshold <- if (length(args) >= 1) as.numeric(args[1]) else 1e-4
egger_int_threshold <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
het_q_threshold     <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
coloc_pph4_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.5

# --- Read each evidence layer (graceful missing) -----------------------------

load_layer <- function(fp_pattern, label) {
  candidates <- c(file.path(dir_processed, paste0(fp_pattern, ".parquet")),
                  list.files(dir_processed,
                             pattern = paste0("^", fp_pattern, ".*\\.parquet$"),
                             full.names = TRUE))
  candidates <- unique(candidates[file.exists(candidates)])
  # Prefer the alias (no timestamp suffix), then most recent timestamped
  if (length(candidates) == 0) {
    cat(sprintf("  [%s] missing — skipping\n", label)); return(NULL)
  }
  fp <- candidates[1]
  d <- as.data.table(read_parquet(fp))
  cat(sprintf("  [%s] loaded: %s (%d rows)\n", label, basename(fp), nrow(d)))
  d
}

cat(sprintf("[%s] Loading evidence layers...\n", Sys.time()))
phase2 <- load_layer("locus_mr_results", "Phase 2 locus MR")
sens   <- load_layer("sensitivity_mr_results", "Phase 3 sensitivity")
pgs_n  <- load_layer("pgs_mr_results_narrow", "PGS-MR narrow")
pgs_w  <- load_layer("pgs_mr_results_wide",   "PGS-MR wide")
coloc  <- load_layer("coloc_results", "Phase 4 coloc")
kim    <- load_layer("kim_korean_locus_mr_results", "Kim Korean secondary primary (Amendment 4 B)")
kim_rep <- load_layer("replication_kim_korean_per_panel", "Kim Korean H4 replication (Amendment 4 A)")

if (is.null(phase2)) {
  stop("Phase 2 locus MR results are required to build the atlas. Run R/05 first.")
}

# Metadata may not be present on the checkpoint (only added at Phase 2 final step);
# merge it in here regardless so the script works on intermediate + final inputs.
if (!all(c("phenotype", "category", "num_cases", "num_controls") %in% names(phase2))) {
  cat(sprintf("  [Phase 2] metadata columns missing — merging from osf/phecode_filter_prespec.csv\n"))
  phecode_meta <- fread(file.path(proj_root, "osf/phecode_filter_prespec.csv"))
  meta_sub <- phecode_meta[, .(phecode = phenocode, phenotype, category,
                                num_cases, num_controls)]
  phase2 <- merge(phase2, meta_sub, by = "phecode", all.x = TRUE, sort = FALSE)
}

# --- Build base atlas from Phase 2 -------------------------------------------
# Reshape locus MR from long (method per row) to wide (one row per phecode × panel)

base <- dcast(phase2[!is.na(method)],
              phecode + panel + n_iv + phenotype + category + num_cases + num_controls ~ method,
              value.var = c("b", "se", "pval"))
# Q + Egger intercept are stored on every method-row in Phase 2; take the per-pair maximum-coverage value
het_egger <- phase2[!is.na(method),
                    .(q_stat = max(q_stat, na.rm = TRUE),
                      q_pval = max(q_pval, na.rm = TRUE),
                      egger_intercept = max(egger_intercept, na.rm = TRUE),
                      egger_intercept_pval = max(egger_intercept_pval, na.rm = TRUE)),
                    by = .(phecode, panel)]
# max-on-all-NA returns -Inf; convert
for (col in c("q_stat", "q_pval", "egger_intercept", "egger_intercept_pval")) {
  het_egger[is.infinite(get(col)), (col) := NA_real_]
}
atlas <- merge(base, het_egger, by = c("phecode", "panel"), all.x = TRUE)

cat(sprintf("[%s] Base atlas: %d rows × %d cols (after Phase 2 reshape)\n",
            Sys.time(), nrow(atlas), ncol(atlas)))

# --- Join Phase 3 sensitivity ------------------------------------------------

if (!is.null(sens) && nrow(sens) > 0) {
  sens_wide <- dcast(sens[!is.na(b)],
                      phecode + panel ~ method,
                      value.var = c("b", "se", "pval"),
                      fun.aggregate = function(x) x[1])
  # Rename to clarify these are sensitivity methods
  setnames(sens_wide,
           old = grep("^(b|se|pval)_", names(sens_wide), value = TRUE),
           new = paste0("sens_", grep("^(b|se|pval)_", names(sens_wide), value = TRUE)))
  atlas <- merge(atlas, sens_wide, by = c("phecode", "panel"), all.x = TRUE)
  cat(sprintf("[%s] Sensitivity merged (Phase 3): atlas now %d cols\n",
              Sys.time(), ncol(atlas)))
}

# --- Join PGS-MR (narrow + wide) --------------------------------------------

if (!is.null(pgs_n) && nrow(pgs_n) > 0) {
  pgs_n_sub <- pgs_n[, .(phecode,
                          pgs_narrow_n_snps = n_snps_overlap,
                          pgs_narrow_b = b,
                          pgs_narrow_se = se,
                          pgs_narrow_pval = pval)]
  atlas <- merge(atlas, pgs_n_sub, by = "phecode", all.x = TRUE)
  cat(sprintf("[%s] PGS-MR narrow merged: atlas now %d cols\n", Sys.time(), ncol(atlas)))
}
if (!is.null(pgs_w) && nrow(pgs_w) > 0) {
  pgs_w_sub <- pgs_w[, .(phecode,
                          pgs_wide_n_snps = n_snps_overlap,
                          pgs_wide_b = b,
                          pgs_wide_se = se,
                          pgs_wide_pval = pval)]
  atlas <- merge(atlas, pgs_w_sub, by = "phecode", all.x = TRUE)
  cat(sprintf("[%s] PGS-MR wide merged: atlas now %d cols\n", Sys.time(), ncol(atlas)))
}

# --- Join Kim Korean secondary primary discovery (Amendment 4 Component B) --

if (!is.null(kim) && nrow(kim) > 0 && "method" %in% names(kim)) {
  # Reshape long → wide by method per (phecode × panel), keep only IVW for Atlas
  kim_ivw <- kim[!is.na(method) & method == "ivw",
                  .(phecode, panel,
                    kim_n_iv = n_iv, kim_b_ivw = b, kim_se_ivw = se,
                    kim_pval_ivw = pval)]
  atlas <- merge(atlas, kim_ivw, by = c("phecode", "panel"), all.x = TRUE)
  # Concordance flag: Kim IVW direction matches Ali EUR IVW direction at same panel
  atlas[, kim_concordant_ivw := !is.na(kim_b_ivw) & !is.na(b_ivw) &
                                   sign(kim_b_ivw) == sign(b_ivw)]
  cat(sprintf("[%s] Kim Korean secondary primary merged: atlas now %d cols\n",
              Sys.time(), ncol(atlas)))
}

# --- Join Phase 4 colocalisation --------------------------------------------

if (!is.null(coloc) && nrow(coloc) > 0) {
  # Per phecode: best PP.H4 across all (tissue, gene), plus count of high-PP.H4 evidence
  coloc_summary <- coloc[, .(
    best_pp_h4 = max(PP.H4, na.rm = TRUE),
    n_coloc_tested = .N,
    n_coloc_pph4_gt0p5 = sum(PP.H4 > 0.5, na.rm = TRUE),
    n_coloc_pph4_gt0p8 = sum(PP.H4 > 0.8, na.rm = TRUE),
    top_coloc_gene = gene[which.max(PP.H4)],
    top_coloc_tissue = tissue[which.max(PP.H4)]
  ), by = phecode]
  atlas <- merge(atlas, coloc_summary, by = "phecode", all.x = TRUE)
  cat(sprintf("[%s] Coloc merged: atlas now %d cols\n", Sys.time(), ncol(atlas)))
}

# --- Compute robust hit flag per pre-reg §28 --------------------------------
# A "robust hit" requires:
#   1. IVW p < threshold (default 1e-4)
#   2. Weighted median direction matches IVW (sign(b_wm) == sign(b_ivw))
#   3. Egger intercept non-significant (intercept p > 0.05)
#   4. Cochran's Q non-significant (Q_p > 0.05) — no excess heterogeneity
#   5. (Optional) coloc evidence: best_pp_h4 > 0.5 if available

atlas[, robust_hit := FALSE]
has_ivw  <- !is.na(atlas$pval_ivw) & atlas$pval_ivw < robust_p_threshold
has_wm   <- !is.na(atlas$b_weighted_median) & !is.na(atlas$b_ivw) &
            sign(atlas$b_weighted_median) == sign(atlas$b_ivw)
has_egger_clean <- is.na(atlas$egger_intercept_pval) |
                   atlas$egger_intercept_pval > egger_int_threshold
has_het_clean   <- is.na(atlas$q_pval) |
                   atlas$q_pval > het_q_threshold
atlas[has_ivw & has_wm & has_egger_clean & has_het_clean, robust_hit := TRUE]

# Coloc-strengthened hit: robust_hit AND coloc PP.H4 > threshold
atlas[, coloc_strengthened := FALSE]
if ("best_pp_h4" %in% names(atlas)) {
  atlas[robust_hit & !is.na(best_pp_h4) & best_pp_h4 > coloc_pph4_threshold,
        coloc_strengthened := TRUE]
}

# --- FDR correction within each panel (BH, across phecodes per panel) -------

# Per pre-reg §12: FDR-BH per panel-method (not pooled across panels)
atlas[, ivw_fdr_q_within_panel := NA_real_]
for (pn in unique(atlas$panel)) {
  idx <- atlas$panel == pn & !is.na(atlas$pval_ivw)
  if (sum(idx) >= 2) {
    atlas[idx, ivw_fdr_q_within_panel := p.adjust(pval_ivw, method = "BH")]
  }
}

# --- Write ------------------------------------------------------------------

out_main <- file.path(dir_processed, "amyloid_pet_causal_atlas.parquet")
out_dated <- file.path(dir_processed,
                       sprintf("amyloid_pet_causal_atlas_%s.parquet",
                                format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(atlas, out_main)
write_parquet(atlas, out_dated)
cat(sprintf("[%s] Wrote %s (%d rows × %d cols)\n",
            Sys.time(), basename(out_main), nrow(atlas), ncol(atlas)))

# Top-hits subset
top_hits <- atlas[robust_hit == TRUE | (!is.na(pval_ivw) & pval_ivw < 1e-6)]
out_top <- file.path(dir_processed, "causal_atlas_top_hits.parquet")
write_parquet(top_hits, out_top)
cat(sprintf("[%s] Wrote %s (%d top-hit rows)\n",
            Sys.time(), basename(out_top), nrow(top_hits)))

# --- Summary report ---------------------------------------------------------

cat("\n=== Causal Atlas summary ===\n")
cat(sprintf("Total atlas rows: %d (panels × phecodes)\n", nrow(atlas)))
cat(sprintf("Robust hits (IVW p<%.0e + WM concordance + clean intercept + clean Q): %d\n",
            robust_p_threshold, sum(atlas$robust_hit, na.rm = TRUE)))
if ("best_pp_h4" %in% names(atlas)) {
  cat(sprintf("Coloc-strengthened robust hits (PP.H4 > %.2f): %d\n",
              coloc_pph4_threshold, sum(atlas$coloc_strengthened, na.rm = TRUE)))
}
cat("\nRobust hits by panel:\n")
print(atlas[robust_hit == TRUE, .N, by = panel])

cat("\nTop 20 robust hits by IVW p-value:\n")
top20 <- atlas[robust_hit == TRUE][order(pval_ivw)][1:min(20, .N),
                .(phecode, panel, phenotype, category,
                  n_iv,
                  b = round(b_ivw, 4),
                  pval = signif(pval_ivw, 3),
                  fdr_q = signif(ivw_fdr_q_within_panel, 3))]
print(top20)
