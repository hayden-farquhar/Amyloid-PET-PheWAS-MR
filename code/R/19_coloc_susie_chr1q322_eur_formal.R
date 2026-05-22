#!/usr/bin/env Rscript
# =============================================================================
# Script 19: Formal amyloid-side coloc.susie at chr1q32.2 — EUR-only Ali 2023
# =============================================================================
#
# v3 follow-on to R/17 (which ran eQTL-side SuSiE only because Ali multi-ethnic
# vs 1KG EUR LD broke amyloid-side SuSiE convergence — documented Hukku 2021
# LD-cohort-mismatch failure mode).
#
# Strategy: use the Ali 2023 NHW (non-Hispanic White, N=11,556) sub-meta —
# the EUR-only subset of the same Ali 2023 meta-analysis — which IS ancestry-
# matched to the 1KG EUR LD reference. This should let amyloid-side SuSiE
# converge and enable formal coloc.susie() at the locus.
#
# Pipeline:
#   1. Read Ali NHW sumstats for chr1q32.2 region (newly downloaded)
#   2. Reuse eQTL data + LD matrix from R/17
#   3. Run coloc::runsusie on amyloid (NHW) and per-gene eQTL via dataset interface
#   4. Run coloc::coloc.susie() pairwise (amyloid × each gene)
#   5. Report per-credible-set PP.H4
#   6. Save to data/processed/coloc_susie_chr1q322_eur_formal.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(coloc)
  library(susieR)
})

proj_root <- "/Users/haydenfarquhar/Documents/Research Projects/81 Amyloid PET PheWAS MR"
dir_interim <- file.path(proj_root, "data/interim/coloc_susie")
dir_processed <- file.path(proj_root, "data/processed")

# Inputs
ALI_NHW <- file.path(proj_root, "data/raw/ali_2023",
  "AmyloidPET_GWAS_NHW_metaAnalysis_cohorts13_n11556_hg38_WashU.txt")
EQTL_REGION <- file.path(dir_interim, "chr1q32.2_brain_cortex_raw.tsv")
LD_PREFIX <- file.path(dir_interim, "chr1q32.2_ld")

# Region
CHR <- 1L
WIN_LO <- 207018704L
WIN_HI <- 208018704L

GENES_OF_INTEREST <- c(
  CR1   = "ENSG00000203710",
  CR1L  = "ENSG00000203709",
  CD46  = "ENSG00000117335"
)

N_ALI_NHW    <- 11556  # NHW sub-meta sample size
N_GTEX_BRAIN <- 205

set.seed(81)

if (!file.exists(ALI_NHW)) {
  stop(sprintf("Ali NHW sumstats not found at %s\nDownload via WashU Box: https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l", ALI_NHW))
}

# --- 1. Extract Ali NHW chr1q32.2 region ------------------------------------
cat(sprintf("[%s] Extracting Ali NHW chr1q32.2 region\n", Sys.time()))
out_nhw_region <- file.path(dir_interim, "ali_nhw_chr1q322.tsv")
if (!file.exists(out_nhw_region)) {
  # Same format as multi-ethnic: CHR BP SNP EFFECT_ALLELE OTHER_ALLELE BETA SE P ALT_FREQS
  awk_cmd <- sprintf(
    "head -1 %s > %s; awk -v c=%d -v l=%d -v h=%d 'NR>1 && $1==c && $2>=l && $2<=h' %s >> %s",
    shQuote(ALI_NHW), shQuote(out_nhw_region),
    CHR, WIN_LO, WIN_HI,
    shQuote(ALI_NHW), shQuote(out_nhw_region))
  system(awk_cmd)
}
ali <- fread(out_nhw_region)
cat(sprintf("  Ali NHW rows in chr1q32.2: %d\n", nrow(ali)))
cat(sprintf("  Columns: %s\n", paste(names(ali), collapse=", ")))

# Verify NHW data structure matches multi-ethnic format
if (!all(c("BP","EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P") %in% names(ali))) {
  stop("Ali NHW file structure differs from expected; column names: ",
       paste(names(ali), collapse=", "))
}
ali[, z_amyloid := BETA / SE]

# --- 2. Reuse eQTL data + LD from R/17 --------------------------------------
cat(sprintf("[%s] Loading eQTL + LD (from R/17)\n", Sys.time()))
eqtl <- fread(EQTL_REGION, header = FALSE, col.names = c(
  "variant_id", "r2", "pvalue", "molecular_trait_id", "molecular_trait_object_id",
  "maf", "gene_id", "median_tpm", "beta", "se",
  "an", "ac", "chr_e", "pos_e", "ref_e", "alt_e", "type", "rsid"
))
eqtl <- eqtl[rsid != "." & !is.na(rsid) & rsid != ""]
eqtl[, z_eqtl := beta / se]

eqtl_goi <- eqtl[molecular_trait_id %in% GENES_OF_INTEREST]
eqtl_goi[, abs_z := abs(z_eqtl)]
setorder(eqtl_goi, rsid, molecular_trait_id, -abs_z)
eqtl_goi_uniq <- eqtl_goi[!duplicated(eqtl_goi, by = c("rsid", "molecular_trait_id"))]
eqtl_wide <- dcast(eqtl_goi_uniq,
                   rsid + chr_e + pos_e + ref_e + alt_e ~ molecular_trait_id,
                   value.var = c("z_eqtl", "beta", "se"))
cat(sprintf("  eqtl_wide rows: %d, cols: %d\n", nrow(eqtl_wide), ncol(eqtl_wide)))

# --- 3. Join Ali NHW + eQTL on rsID (via chr+pos in Ali) --------------------
ali_key <- ali[, .(BP, EA = EFFECT_ALLELE, OA = OTHER_ALLELE,
                    BETA, SE, z_amyloid)]
eqtl_pos <- eqtl_wide[, .(rsid, BP = pos_e, ref_e, alt_e)]
ali_w_rsid <- merge(ali_key, eqtl_pos, by = "BP", allow.cartesian = TRUE)
ali_w_rsid[, orientation := fcase(
  EA == alt_e & OA == ref_e, "match",
  EA == ref_e & OA == alt_e, "flip",
  default = "drop"
)]
cat("  Ali-eQTL orientation distribution:\n")
print(table(ali_w_rsid$orientation))
ali_w_rsid <- ali_w_rsid[orientation != "drop"]
ali_w_rsid[orientation == "flip", z_amyloid := -z_amyloid]
ali_w_rsid[orientation == "flip", BETA := -BETA]
ali_compact <- unique(ali_w_rsid[, .(rsid, BETA, SE, z_amyloid,
                                      ali_alt = alt_e, ali_ref = ref_e)],
                     by = "rsid")

merged <- merge(eqtl_wide, ali_compact, by = "rsid")
cat(sprintf("  Merged Ali NHW + eQTL rows: %d\n", nrow(merged)))

# --- 4. Bring in LD matrix (rsID-keyed; reuses R/17 LD extraction) ----------
ld_bim <- fread(paste0(LD_PREFIX, ".bim"), header = FALSE,
                col.names = c("chr", "rsid_ld", "cm", "bp_ld", "a1_ld", "a2_ld"))
ld_bim[, ld_row := .I]
merged <- merge(merged, ld_bim[, .(rsid = rsid_ld, ld_row, a1_ld, a2_ld)],
                by = "rsid")
cat(sprintf("  After LD rsID match: %d SNPs\n", nrow(merged)))

merged[, ld_alt_orient := fcase(
  a1_ld == alt_e & a2_ld == ref_e,  1L,
  a1_ld == ref_e & a2_ld == alt_e, -1L,
  default = NA_integer_
)]
n_amb <- merged[is.na(ld_alt_orient), .N]
cat(sprintf("  Dropped SNPs with LD allele ambiguity: %d\n", n_amb))
merged <- merged[!is.na(ld_alt_orient)]
cat(sprintf("  Final SNP count: %d\n", nrow(merged)))

stopifnot(nrow(merged) >= 20)

# --- 5. Extract LD submatrix + sign-align -----------------------------------
cat(sprintf("[%s] Loading LD submatrix\n", Sys.time()))
ld_full <- as.matrix(fread(paste0(LD_PREFIX, ".ld"), header = FALSE))
ld_sub <- ld_full[merged$ld_row, merged$ld_row]
sign_vec <- merged$ld_alt_orient
ld_sub <- ld_sub * outer(sign_vec, sign_vec)
ld_sub <- (ld_sub + t(ld_sub)) / 2
diag(ld_sub) <- 1
rownames(ld_sub) <- merged$rsid
colnames(ld_sub) <- merged$rsid

# Diagnostic
if ("rs6656401" %in% merged$rsid) {
  idx <- which(merged$rsid == "rs6656401")
  cat(sprintf("  Diagnostic: rs6656401 at idx %d; z_amyloid=%.3f (Ali NHW), p=%.2e\n",
              idx, merged$z_amyloid[idx],
              2*pnorm(-abs(merged$z_amyloid[idx]))))
}
cat(sprintf("  Range of |z_amyloid|: %.2f - %.2f\n",
            min(abs(merged$z_amyloid)), max(abs(merged$z_amyloid))))

# --- 6a. SuSiE-RSS diagnostic — identify problematic SNPs --------------------
# Per https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
# When in-sample LD differs from reference LD (Ali NHW is a meta of 13 cohorts;
# 1KG EUR is 503 unrelated individuals), strong-signal SNPs can have z-scores
# inconsistent with the LD pattern, causing SuSiE to fail. The diagnostic
# computes expected z-scores from neighbouring SNPs' LD-weighted z and compares
# to observed z; outliers are problematic SNPs that should be filtered or
# down-weighted.
cat(sprintf("[%s] SuSiE-RSS diagnostic: identifying LD-inconsistent SNPs\n", Sys.time()))
z <- merged$z_amyloid
kriging_res <- susieR::kriging_rss(z = z, R = ld_sub, n = N_ALI_NHW)
# Outliers: SNPs whose observed z differs from kriged-z by > 3 SD
diag_df <- kriging_res$conditional_dist
diag_df$rsid <- merged$rsid
diag_df$logLR <- diag_df$logLR
diag_df$z_obs <- z
diag_df$abs_diff <- abs(diag_df$z_obs - diag_df$condmean)
# Flag SNPs with abs_diff > 3 OR logLR > 2 (per susieR vignette suggestion)
diag_df$outlier <- (diag_df$abs_diff > 3) | (diag_df$logLR > 2)
n_outliers <- sum(diag_df$outlier)
cat(sprintf("  Total SNPs: %d; LD-inconsistent (outliers): %d\n", nrow(diag_df), n_outliers))
if (n_outliers > 0) {
  cat("  Top 10 outliers by logLR:\n")
  outlier_df <- diag_df[order(-diag_df$logLR), ][1:min(10, n_outliers), ]
  print(outlier_df[, c("rsid", "z_obs", "condmean", "abs_diff", "logLR")])
}

# --- 6b. SuSiE on amyloid: ALWAYS drop kriging outliers first ---------------
# Per the diagnostic above, LD-inconsistent SNPs (whose observed z differs
# from the LD-kriged expectation) cause SuSiE to place "phantom" credible
# sets at the outlier positions rather than at the true signal. The right
# procedure is to drop outliers BEFORE the first SuSiE call, not after a
# failed attempt.
if (n_outliers > 0) {
  cat(sprintf("[%s] Dropping %d kriging-outlier SNPs before SuSiE\n",
              Sys.time(), n_outliers))
  keep_idx <- !diag_df$outlier
  merged <- merged[keep_idx]
  ld_sub <- ld_sub[keep_idx, keep_idx]
  rownames(ld_sub) <- merged$rsid
  colnames(ld_sub) <- merged$rsid
  cat(sprintf("  Remaining SNPs: %d\n", nrow(merged)))
  # Re-check that rs6656401 (the AD GWAS lead) is still present
  if ("rs6656401" %in% merged$rsid) {
    idx <- which(merged$rsid == "rs6656401")
    cat(sprintf("  ✓ rs6656401 retained at idx %d (z=%.2f)\n",
                idx, merged$z_amyloid[idx]))
  }
}

# Try a progression of L values from conservative to permissive
ds_amy <- list(
  beta = merged$BETA, varbeta = merged$SE^2,
  N = N_ALI_NHW, sdY = 1, type = "quant",
  LD = ld_sub, snp = merged$rsid, position = merged$pos_e
)
susie_amy <- NULL
for (L_try in c(1L, 2L, 3L)) {
  cat(sprintf("[%s] runsusie on amyloid-PET NHW (L=%d, no resvar estimation, refine=TRUE)\n",
              Sys.time(), L_try))
  attempt <- tryCatch(
    coloc::runsusie(ds_amy, suffix = sprintf("amyloid_nhw_L%d", L_try),
                    L = L_try, maxit = 1000,
                    estimate_residual_variance = FALSE,
                    refine = TRUE),
    error = function(e) {
      cat(sprintf("  runsusie L=%d error: %s\n", L_try, conditionMessage(e))); NULL }
  )
  if (!is.null(attempt) && !is.null(attempt$sets$cs) && length(attempt$sets$cs) > 0) {
    # Check that lead SNP of CS 1 is near rs6656401 (sanity check)
    cs_leads <- sapply(attempt$sets$cs, function(s) {
      lead_idx <- s[which.max(attempt$pip[s])]
      merged$pos_e[lead_idx]
    })
    dist_to_amyloid_lead_kb <- abs(cs_leads - 207518704) / 1000
    cat(sprintf("  L=%d gave %d credible set(s); lead positions: %s; dist from rs6656401: %s kb\n",
                L_try, length(attempt$sets$cs),
                paste(cs_leads, collapse=", "),
                paste(sprintf("%.1f", dist_to_amyloid_lead_kb), collapse=", ")))
    # If any CS lead is within 100 kb of rs6656401, accept this result
    if (any(dist_to_amyloid_lead_kb < 100)) {
      susie_amy <- attempt
      cat(sprintf("  ✓ At least one CS lead is within 100 kb of rs6656401 — accepting L=%d result\n", L_try))
      break
    } else {
      cat(sprintf("  ✗ All CS leads are >100 kb from rs6656401 — LD pathology; trying next L\n"))
    }
  }
}
if (is.null(susie_amy)) {
  cat("FORMAL coloc.susie failed at amyloid side: no L produced a CS lead near rs6656401.\n")
  cat("This is the documented LD-cohort-mismatch failure mode (Hukku 2021) extending to\n")
  cat("single-ancestry-matched meta-analytical GWAS against single-cohort reference LD.\n")
  diag_out_fp <- file.path(dir_processed, "coloc_susie_diagnostic_nhw.tsv")
  fwrite(diag_df, diag_out_fp, sep = "\t")
  cat(sprintf("  Diagnostic saved to %s\n", diag_out_fp))
  quit(status = 2)
}
cs_amy <- susie_amy$sets$cs
cat(sprintf("  Amyloid NHW credible sets: %d\n",
            if (is.null(cs_amy)) 0 else length(cs_amy)))
if (!is.null(cs_amy)) {
  for (i in seq_along(cs_amy)) {
    snps_i <- cs_amy[[i]]
    lead_idx <- snps_i[which.max(susie_amy$pip[snps_i])]
    cat(sprintf("    CS %d: %d SNPs; lead %s (chr%d:%d) PIP %.3f\n",
                i, length(snps_i),
                merged$rsid[lead_idx],
                merged$chr_e[lead_idx],
                merged$pos_e[lead_idx],
                max(susie_amy$pip[snps_i])))
  }
}

# --- 7. SuSiE per gene + coloc.susie -----------------------------------------
cat(sprintf("[%s] coloc.susie per gene\n", Sys.time()))
results <- list()
for (g_name in names(GENES_OF_INTEREST)) {
  g_id <- GENES_OF_INTEREST[[g_name]]
  beta_col <- paste0("beta_", g_id)
  se_col   <- paste0("se_", g_id)
  if (!(beta_col %in% names(merged) && se_col %in% names(merged))) next
  b_g <- merged[[beta_col]]; s_g <- merged[[se_col]]
  ok <- !is.na(b_g) & !is.na(s_g) & s_g > 0
  if (sum(ok) < 20) next
  b_g[!ok] <- 0; s_g[!ok] <- median(s_g[ok])

  cat(sprintf("\n  --- Gene %s (%s) ---\n", g_name, g_id))
  ds_e <- list(beta = b_g, varbeta = s_g^2, N = N_GTEX_BRAIN, sdY = 1, type = "quant",
                LD = ld_sub, snp = merged$rsid, position = merged$pos_e)
  susie_e <- tryCatch(
    coloc::runsusie(ds_e, suffix = g_name, maxit = 500),
    error = function(e) { cat("    runsusie error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(susie_e)) next
  cs_e <- susie_e$sets$cs
  cat(sprintf("    eQTL credible sets: %d\n", if (is.null(cs_e)) 0 else length(cs_e)))

  cs_res <- tryCatch(
    coloc::coloc.susie(susie_amy, susie_e),
    error = function(e) { cat("    coloc.susie error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(cs_res) || is.null(cs_res$summary)) next
  smry <- as.data.table(cs_res$summary)
  smry[, `:=`(gene = g_name, gene_id = g_id)]
  cat("    coloc.susie summary:\n")
  print(smry[, .(gene, idx1, idx2, hit1, hit2,
                  PP.H4 = round(PP.H4.abf, 4),
                  PP.H3 = round(PP.H3.abf, 4))])
  results[[g_name]] <- smry
}

# --- 8. Save -----------------------------------------------------------------
if (length(results) > 0) {
  out <- rbindlist(results, fill = TRUE)
  out_fp <- file.path(dir_processed, "coloc_susie_chr1q322_eur_formal.parquet")
  write_parquet(out, out_fp)
  fwrite(out, sub("\\.parquet$", ".tsv", out_fp), sep = "\t")
  cat(sprintf("\nWrote %s (%d rows)\n", out_fp, nrow(out)))
} else {
  cat("\nNo coloc.susie results to write.\n")
}
cat(sprintf("\n[%s] Done.\n", Sys.time()))
