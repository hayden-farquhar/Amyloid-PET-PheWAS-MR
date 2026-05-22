#!/usr/bin/env Rscript
# =============================================================================
# Script 19b: Cross-check — run coloc.abf with Ali NHW (EUR-only) sumstats
# to see if PP.H4 is preserved or differs from the multi-ethnic Ali run.
# coloc.abf doesn't use LD, so this isolates the LD-related component of the
# formal coloc.susie result.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

proj_root <- Sys.getenv("AMYLOID_PHEWAS_MR_ROOT", ".")
ali_nhw_fp <- file.path(proj_root, "data/interim/coloc_susie/ali_nhw_chr1q322.tsv")
ali_me_fp  <- file.path(proj_root, "data/interim/coloc_susie/ali_chr1q322.tsv")
eqtl_fp    <- file.path(proj_root, "data/interim/coloc_susie/chr1q32.2_brain_cortex_raw.tsv")

CHR <- 1L
WIN_LO <- 207018704L
WIN_HI <- 208018704L

GENES_OF_INTEREST <- c(
  CR1   = "ENSG00000203710",
  CR1L  = "ENSG00000203709",
  CD46  = "ENSG00000117335"
)

N_ALI_NHW <- 11556
N_ALI_ME  <- 13409
N_GTEX_BRAIN_CORTEX <- 205

# --- Read both Ali sumstats and the eQTL data --------------------------------
ali_nhw <- fread(ali_nhw_fp)
ali_me  <- fread(ali_me_fp)
# Coerce ALT_FREQS to numeric — Ali files sometimes have it as character
ali_nhw[, ALT_FREQS := as.numeric(ALT_FREQS)]
ali_me[, ALT_FREQS := as.numeric(ALT_FREQS)]
# If MAF is missing entirely, set to 0.1 (placeholder; coloc.abf is more
# sensitive to relative MAF than absolute value)
ali_nhw[is.na(ALT_FREQS), ALT_FREQS := 0.1]
ali_me[is.na(ALT_FREQS),  ALT_FREQS := 0.1]
eqtl <- fread(eqtl_fp, header = FALSE, col.names = c(
  "variant_id", "r2", "pvalue", "molecular_trait_id", "molecular_trait_object_id",
  "maf", "gene_id", "median_tpm", "beta", "se",
  "an", "ac", "chr_e", "pos_e", "ref_e", "alt_e", "type", "rsid"
))

# --- ALSO load GTEx significant-pairs (matches original Phase 4 Ali run) ----
# This is the dataset that gave PP.H4 = 0.983 in the original analysis.
gtex_signif_fp <- file.path(proj_root,
  "data/reference/gtex_v8/brain/eqtls/Brain_Cortex.v8.EUR.signif_pairs.txt.gz")
gtex_signif <- fread(gtex_signif_fp)
# variant_id format: chr1_207018728_T_C_b38 → parse to chr, pos, ref, alt
gtex_signif[, c("vid_chr","vid_pos","vid_ref","vid_alt","vid_build") :=
              tstrsplit(variant_id, "_", fixed = TRUE)]
gtex_signif[, vid_chr := as.integer(sub("chr","",vid_chr))]
gtex_signif[, vid_pos := as.integer(vid_pos)]
# Strip Ensembl version suffix (.X) from phenotype_id to match the
# molecular_trait_id format in the full-pairs eQTL data
gtex_signif[, gene_stable := sub("\\..*","", phenotype_id)]
# Restrict to chr1q32.2 window
gtex_signif <- gtex_signif[vid_chr == 1 & vid_pos >= WIN_LO & vid_pos <= WIN_HI]
cat(sprintf("GTEx Brain_Cortex significant-pairs in chr1q32.2: %d rows; %d unique genes\n",
            nrow(gtex_signif), uniqueN(gtex_signif$gene_stable)))

# --- coloc.abf wrapper (using GTEx significant-pairs, matching Phase 4) -----
run_coloc_abf <- function(ali, n_ali, label, gene_id) {
  # USE GTEx SIGNIFICANT-PAIRS — same dataset as the original Phase 4 run
  eqtl_g <- gtex_signif[gene_stable == gene_id,
                         .(pos_e = vid_pos, ref_e = vid_ref, alt_e = vid_alt,
                           beta_e = slope, se_e = slope_se, maf_e = maf)]
  if (nrow(eqtl_g) == 0) {
    cat(sprintf("  [%s × %s] no eQTL rows\n", label, gene_id))
    return(NULL)
  }
  # Merge on chr+pos (both GRCh38)
  m <- merge(
    ali[, .(BP, EA = EFFECT_ALLELE, OA = OTHER_ALLELE,
             BETA_a = BETA, SE_a = SE, EAF = ALT_FREQS)],
    eqtl_g,
    by.x = "BP", by.y = "pos_e", allow.cartesian = TRUE
  )
  # Allele orientation (Ali EA vs eQTL alt)
  m[, orient := fcase(
    EA == alt_e & OA == ref_e, "match",
    EA == ref_e & OA == alt_e, "flip",
    default = "drop"
  )]
  m <- m[orient != "drop"]
  m[orient == "flip", BETA_a := -BETA_a]
  m[orient == "flip", EAF := 1 - EAF]
  m <- m[!is.na(BETA_a) & !is.na(SE_a) & SE_a > 0 &
          !is.na(beta_e) & !is.na(se_e) & se_e > 0]
  # Deduplicate by position — Ali sometimes has multi-allelic SNPs that produce
  # multiple rows at the same BP; keep the row with the strongest amyloid signal
  m[, abs_z := abs(BETA_a / SE_a)]
  setorder(m, BP, -abs_z)
  m <- m[!duplicated(m, by = "BP")]
  if (nrow(m) < 10) {
    cat(sprintf("  [%s × %s] too few SNPs after dedup: %d\n", label, gene_id, nrow(m)))
    return(NULL)
  }
  ds_amy  <- list(beta = m$BETA_a, varbeta = m$SE_a^2, type = "quant", N = n_ali,
                  snp = paste0("s_", m$BP), MAF = pmin(m$EAF, 1 - m$EAF))
  ds_eqtl <- list(beta = m$beta_e, varbeta = m$se_e^2, type = "quant",
                  N = N_GTEX_BRAIN_CORTEX, snp = paste0("s_", m$BP),
                  MAF = pmin(m$maf_e, 1 - m$maf_e))
  res <- coloc::coloc.abf(ds_amy, ds_eqtl)
  s <- res$summary
  cat(sprintf("  [%s × %s] n=%d  PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.4f\n",
              label, gene_id, length(ds_amy$beta),
              s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"],
              s["PP.H3.abf"], s["PP.H4.abf"]))
  invisible(s)
}

cat("=== coloc.abf at chr1q32.2 — Ali NHW (EUR, N=11,556) ===\n")
for (g_name in names(GENES_OF_INTEREST)) {
  g_id <- GENES_OF_INTEREST[[g_name]]
  cat(sprintf("\n--- Gene %s (%s) ---\n", g_name, g_id))
  run_coloc_abf(ali_nhw, N_ALI_NHW, paste0("NHW-", g_name), g_id)
}

cat("\n\n=== coloc.abf at chr1q32.2 — Ali multi-ethnic (N=13,409, for comparison) ===\n")
for (g_name in names(GENES_OF_INTEREST)) {
  g_id <- GENES_OF_INTEREST[[g_name]]
  cat(sprintf("\n--- Gene %s (%s) ---\n", g_name, g_id))
  run_coloc_abf(ali_me, N_ALI_ME, paste0("ME-", g_name), g_id)
}
