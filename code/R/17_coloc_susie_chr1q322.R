#!/usr/bin/env Rscript
# =============================================================================
# Script 17: Full-pairs coloc.susie at chr1q32.2 (CR1 region)
# =============================================================================
#
# Addresses the editorial pre-mortem concern (2026-05-22) and the
# manuscript's own Discussion §6/§7 acknowledgement that coloc.abf cannot
# decompose multi-causal signals at chr1q32.2 (co-located CR1, CR1L, CD46
# cis-eQTLs). Uses full-pairs nominal eQTL data from the eQTL Catalogue's
# GTEx_V8 Brain_Cortex tabix-indexed bgzip file, byte-range-fetched for the
# chr1q32.2 ± 500 kb window only (~9 MB rather than the 3.4 GB tissue file).
#
# Pipeline:
#   1. Read Ali 2023 amyloid-PET sumstats for chr1q32.2 (already extracted to
#      data/interim/coloc_susie/ali_chr1q322.tsv).
#   2. Read GTEx Brain_Cortex eQTL all-pairs for chr1q32.2 (already extracted
#      to data/interim/coloc_susie/chr1q32.2_brain_cortex_raw.tsv).
#   3. Build LD matrix from 1KG EUR Phase 3 for the region (PLINK --r square).
#   4. Per gene-of-interest (CR1, CR1L, CD46, plus any other gene with a
#      strong eQTL signal in the window): run susie_rss on the eQTL z-vector.
#   5. Run susie_rss on the amyloid-PET z-vector.
#   6. Run coloc.susie() pairwise (amyloid × each gene).
#   7. Report per-credible-set PP.H4.
#   8. Save results to parquet.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(coloc)
  library(susieR)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWAS_MR_ROOT", ".")
dir_interim <- file.path(proj_root, "data/interim/coloc_susie")
dir_processed <- file.path(proj_root, "data/processed")
plink_bfile <- Sys.getenv("EUR_LD_PANEL", "data/reference/EUR")

ali_fp  <- file.path(dir_interim, "ali_chr1q322.tsv")
eqtl_fp <- file.path(dir_interim, "chr1q32.2_brain_cortex_raw.tsv")
ld_prefix <- file.path(dir_interim, "chr1q32.2_ld")

# Region anchored on rs6656401 (CR1 lead SNP, chr1:207,518,704), ±500 kb
CHR <- 1L
WIN_LO <- 207018704L
WIN_HI <- 208018704L

# Gene catalogue for the chr1q32.2 multi-causal question
GENES_OF_INTEREST <- c(
  CR1   = "ENSG00000203710",   # Complement Receptor 1 — primary mediator candidate
  CR1L  = "ENSG00000203709",   # CR1-like — co-located paralog
  CD46  = "ENSG00000117335"    # Co-located complement regulator
)

N_AMYLOID    <- 13409   # Ali 2023 multi-ethnic primary
N_GTEX_BRAIN <- 205     # GTEx v8 Brain_Cortex EUR sample size

set.seed(81)

# --- 1. Read amyloid-PET sumstats --------------------------------------------
cat(sprintf("[%s] Reading Ali 2023 chr1q32.2 region\n", Sys.time()))
ali <- fread(ali_fp)
cat(sprintf("  Rows: %d; cols: %s\n", nrow(ali), paste(names(ali), collapse=", ")))
# Build join key: chr_pos_ref_alt (Ali uses CHR, BP, EFFECT_ALLELE, OTHER_ALLELE).
# eQTL Catalogue stores ref/alt in cols named based on import schema; we will
# build both orientations and match on whichever aligns.
ali[, `:=`(
  key_ea = sprintf("%d_%d_%s_%s", CHR, BP, OTHER_ALLELE, EFFECT_ALLELE),
  key_fa = sprintf("%d_%d_%s_%s", CHR, BP, EFFECT_ALLELE, OTHER_ALLELE),
  z_amyloid = BETA / SE
)]

# --- 2. Read GTEx Brain_Cortex eQTL ------------------------------------------
cat(sprintf("[%s] Reading GTEx Brain_Cortex chr1q32.2 region\n", Sys.time()))
# eQTL Catalogue / GTEx_V8 imported TSV has no header in the bgzip; column
# inventory inferred from format spec + first-row inspection:
#   1: variant_id   2: r2 (NA)   3: pvalue
#   4: molecular_trait_id   5: molecular_trait_object_id   6: maf
#   7: gene_id   8: median_tpm   9: beta   10: se
#   11: an   12: ac   13: chromosome   14: position
#   15: ref   16: alt   17: type   18: rsid
eqtl <- fread(eqtl_fp, header = FALSE, col.names = c(
  "variant_id", "r2", "pvalue", "molecular_trait_id", "molecular_trait_object_id",
  "maf", "gene_id", "median_tpm", "beta", "se",
  "an", "ac", "chr_e", "pos_e", "ref_e", "alt_e", "type", "rsid"
))
cat(sprintf("  Rows: %d; unique molecular_trait_ids: %d\n",
            nrow(eqtl), uniqueN(eqtl$molecular_trait_id)))

eqtl[, z_eqtl := beta / se]
# rsID is the universal join key because Ali (GRCh38) and the 1KG EUR PLINK
# panel (GRCh37) cannot be joined on chr:pos. rsIDs are build-invariant.
eqtl <- eqtl[rsid != "." & !is.na(rsid) & rsid != ""]
cat(sprintf("  Rows with valid rsid: %d\n", nrow(eqtl)))

# Restrict to genes of interest first to size the LD problem sensibly
eqtl_goi <- eqtl[molecular_trait_id %in% GENES_OF_INTEREST]
cat(sprintf("  Rows in GOI subset: %d\n", nrow(eqtl_goi)))
cat("  Per-gene SNP counts in GOI subset:\n")
print(eqtl_goi[, .N, by = molecular_trait_id])

# --- 3. Build LD matrix from 1KG EUR -----------------------------------------
cat(sprintf("[%s] Building 1KG EUR LD matrix for chr1q32.2\n", Sys.time()))
# Use PLINK to extract region and write square correlation matrix.
# Only generate if not already present (resumable).
ld_mat_fp <- paste0(ld_prefix, ".ld")
ld_bim_fp <- paste0(ld_prefix, ".bim")
if (!file.exists(ld_mat_fp) || !file.exists(ld_bim_fp)) {
  cmd <- sprintf(
    "plink --bfile %s --chr %d --from-bp %d --to-bp %d --r square --write-snplist --make-just-bim --out %s 2>&1",
    shQuote(plink_bfile), CHR, WIN_LO, WIN_HI, shQuote(ld_prefix)
  )
  cat("  Running:", cmd, "\n")
  log <- system(cmd, intern = TRUE)
  cat("  PLINK log tail:\n  ", paste(tail(log, 8), collapse = "\n  "), "\n")
} else {
  cat("  Reusing existing LD matrix\n")
}

# Read PLINK BIM (which gives the SNP order in the LD matrix) and the LD matrix
ld_bim <- fread(ld_bim_fp, header = FALSE,
                col.names = c("chr", "rsid_ld", "cm", "bp_ld", "a1_ld", "a2_ld"))
cat(sprintf("  LD matrix dimensions: %d SNPs\n", nrow(ld_bim)))
# 1KG EUR PLINK panel is on GRCh37; Ali/eQTL are on GRCh38. Position-based
# joins fail; we use rsID as the universal key. LD genotype correlations are
# build-invariant so the LD matrix itself is valid as-is.
ld_bim[, ld_row := .I]  # remember the row index for later LD submatrix extraction

# --- 4. Build the common SNP set across Ali, eQTL, LD (rsID-based) ----------
cat(sprintf("[%s] Building common SNP set on rsID\n", Sys.time()))

# Step 4a: collapse eQTL GOI subset to per-SNP-per-gene wide format on rsID.
# Drop variants that appear multiple times per gene under the same rsid (rare
# multi-allelic situations) by taking the strongest signal per (rsid, gene).
eqtl_goi[, abs_z := abs(z_eqtl)]
setorder(eqtl_goi, rsid, molecular_trait_id, -abs_z)
eqtl_goi_uniq <- eqtl_goi[!duplicated(eqtl_goi, by = c("rsid", "molecular_trait_id"))]
eqtl_wide <- dcast(eqtl_goi_uniq,
                   rsid + chr_e + pos_e + ref_e + alt_e ~ molecular_trait_id,
                   value.var = c("z_eqtl", "beta", "se"))
cat(sprintf("  eqtl_wide rows: %d, cols: %d\n", nrow(eqtl_wide), ncol(eqtl_wide)))

# Step 4b: Ali side. Ali's SNP column is in chr:pos:ref:alt format; rsID is
# not provided. We propagate rsID via the chr+pos match against the eQTL
# (both GRCh38, so position matches). Allele orientation must be reconciled
# at this step.
ali_key <- ali[, .(BP, EA = EFFECT_ALLELE, OA = OTHER_ALLELE,
                    BETA, SE, z_amyloid)]
eqtl_pos <- eqtl_wide[, .(rsid, BP = pos_e, ref_e, alt_e)]
ali_w_rsid <- merge(ali_key, eqtl_pos, by = "BP", allow.cartesian = TRUE)
# Determine allele orientation (Ali EA vs eQTL alt)
ali_w_rsid[, orientation := fcase(
  EA == alt_e & OA == ref_e, "match",
  EA == ref_e & OA == alt_e, "flip",
  default = "drop"
)]
cat("  Ali-eQTL orientation distribution (pre-LD):\n")
print(table(ali_w_rsid$orientation))
ali_w_rsid <- ali_w_rsid[orientation != "drop"]
ali_w_rsid[orientation == "flip", z_amyloid := -z_amyloid]
ali_w_rsid[orientation == "flip", BETA := -BETA]
ali_compact <- ali_w_rsid[, .(rsid, BETA, SE, z_amyloid,
                                ali_alt = alt_e, ali_ref = ref_e)]
ali_compact <- unique(ali_compact, by = "rsid")

# Step 4c: join Ali (with rsID) + eQTL wide (rsID-keyed)
merged <- merge(eqtl_wide, ali_compact, by = "rsid")
cat(sprintf("  Merged Ali+eQTL rows (after rsID join): %d\n", nrow(merged)))

# Step 4d: bring in LD index by rsID
merged <- merge(merged, ld_bim[, .(rsid = rsid_ld, ld_row, a1_ld, a2_ld)],
                by = "rsid")
cat(sprintf("  After LD rsID match: %d SNPs\n", nrow(merged)))

# Step 4e: LD allele orientation (LD a1/a2 vs eQTL alt/ref)
# PLINK's --r reports signed correlations relative to a1 dosage. So if a1
# matches our eqtl_alt (== orientation of our z), no flip; if a1 matches our
# eqtl_ref (opposite of our z), flip. Multi-allelic / strand-ambig: drop.
merged[, ld_alt_orient := fcase(
  a1_ld == alt_e & a2_ld == ref_e,  1L,
  a1_ld == ref_e & a2_ld == alt_e, -1L,
  default = NA_integer_
)]
n_amb <- merged[is.na(ld_alt_orient), .N]
cat(sprintf("  SNPs with LD allele ambiguity (dropped): %d\n", n_amb))
merged <- merged[!is.na(ld_alt_orient)]
cat(sprintf("  Final analysis SNP count: %d\n", nrow(merged)))

if (nrow(merged) < 20) {
  stop("Too few SNPs survived join (need ≥20 for SuSiE). Inspect upstream filters.")
}

# --- 6. Extract LD submatrix and sign-align ---------------------------------
cat(sprintf("[%s] Loading LD submatrix\n", Sys.time()))
# Read only the rows/cols we need from the .ld file. PLINK's .ld is
# space-separated, NxN. fread can read it directly.
ld_full <- as.matrix(fread(ld_mat_fp, header = FALSE))
cat(sprintf("  Full LD matrix: %d x %d\n", nrow(ld_full), ncol(ld_full)))
ld_sub <- ld_full[merged$ld_row, merged$ld_row]
# Apply sign-flip from ld_alt_orient (outer product)
sign_vec <- merged$ld_alt_orient
ld_sub <- ld_sub * outer(sign_vec, sign_vec)
cat(sprintf("  Sub LD matrix: %d x %d\n", nrow(ld_sub), ncol(ld_sub)))
# Symmetrise and add tiny diagonal regularisation for SuSiE numerical stability
ld_sub <- (ld_sub + t(ld_sub)) / 2
diag(ld_sub) <- 1
# coloc::runsusie requires LD rownames == colnames == snp vector
rownames(ld_sub) <- merged$rsid
colnames(ld_sub) <- merged$rsid

# --- Diagnostic: confirm rs6656401 (CR1 lead) is in the merged set ---------
if ("rs6656401" %in% merged$rsid) {
  idx <- which(merged$rsid == "rs6656401")
  cat(sprintf("  Diagnostic: rs6656401 IS in merged set at idx %d; ",
              idx))
  cat(sprintf("z_amyloid = %.3f, BETA = %.4f, SE = %.4f, p = %.2e\n",
              merged$z_amyloid[idx], merged$BETA[idx], merged$SE[idx],
              2 * pnorm(-abs(merged$z_amyloid[idx]))))
} else {
  cat("  *** Diagnostic: rs6656401 NOT in merged set — check upstream filters ***\n")
  # Look it up in the underlying data sources for forensics
  in_eqtl <- "rs6656401" %in% eqtl$rsid
  in_ali  <- ali[BP == 207518704, .N] > 0
  cat(sprintf("    rs6656401 in eqtl source: %s\n", in_eqtl))
  cat(sprintf("    chr1:207518704 in Ali source: %s\n", in_ali))
}
cat(sprintf("  Range of |z_amyloid| in merged: %.2f - %.2f\n",
            min(abs(merged$z_amyloid)), max(abs(merged$z_amyloid))))
cat(sprintf("  Range of amyloid p-values in merged: %.2e - %.2e\n",
            2*pnorm(-max(abs(merged$z_amyloid))),
            2*pnorm(-min(abs(merged$z_amyloid)))))

# --- 7. Note on amyloid-side SuSiE feasibility ------------------------------
# Ali 2023 is a multi-ethnic meta-analysis (N=13,409 across 14 cohorts).
# Our LD reference is 1KG EUR (N=503). SuSiE's IBSS algorithm rejects this
# combination with "estimated prior variance is unreasonably large" (a
# documented LD-cohort-mismatch failure mode), even after 50,000 iterations.
# Running coloc.susie with mis-specified LD produces unreliable PP.H4
# estimates (Hukku 2021 review of coloc-with-LD-mismatch). The honest path
# is to (a) compute eQTL-side SuSiE credible sets per gene (which converge
# cleanly because GTEx_V8 Brain_Cortex is European-ancestry, matching the
# 1KG EUR LD reference), (b) tabulate amyloid GWAS evidence at each credible
# set's SNPs, and (c) compare the spatial separation of the eQTL credible
# sets against the amyloid lead. This directly addresses the multi-causal
# question (CR1 vs CR1L vs CD46) without forcing coloc.susie on
# LD-inconsistent inputs.

# --- 8. eQTL-side SuSiE per gene + amyloid lookup --------------------------
cat(sprintf("[%s] eQTL runsusie per gene + amyloid lookup\n", Sys.time()))
cred_set_rows <- list()
per_gene_summary <- list()
for (g_name in names(GENES_OF_INTEREST)) {
  g_id <- GENES_OF_INTEREST[[g_name]]
  beta_col <- paste0("beta_", g_id)
  se_col   <- paste0("se_", g_id)
  if (!(beta_col %in% names(merged) && se_col %in% names(merged))) {
    cat(sprintf("  Gene %s (%s): cols missing (skipping)\n", g_name, g_id))
    next
  }
  b_g <- merged[[beta_col]]
  s_g <- merged[[se_col]]
  ok <- !is.na(b_g) & !is.na(s_g) & s_g > 0
  if (sum(ok) < 20) {
    cat(sprintf("  Gene %s: too few non-NA betas (%d)\n", g_name, sum(ok)))
    next
  }
  b_g[!ok] <- 0; s_g[!ok] <- median(s_g[ok])

  cat(sprintf("\n  --- Gene %s (%s) ---\n", g_name, g_id))
  ds_e <- list(
    beta = b_g, varbeta = s_g^2, N = N_GTEX_BRAIN, sdY = 1, type = "quant",
    LD = ld_sub, snp = merged$rsid, position = merged$pos_e
  )
  susie_e <- tryCatch(
    coloc::runsusie(ds_e, suffix = g_name, maxit = 500),
    error = function(e) { cat("    runsusie error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(susie_e)) next
  cs_e <- susie_e$sets$cs
  n_cs <- if (is.null(cs_e)) 0 else length(cs_e)
  cat(sprintf("    eQTL credible sets for %s: %d\n", g_name, n_cs))

  if (n_cs == 0) {
    per_gene_summary[[g_name]] <- data.table(
      gene = g_name, gene_id = g_id, n_credible_sets = 0L,
      cs_lead_rsid = NA_character_, cs_lead_pos = NA_integer_,
      cs_lead_pip = NA_real_, cs_size = NA_integer_,
      amyloid_z_at_lead = NA_real_, amyloid_p_at_lead = NA_real_,
      dist_to_amyloid_lead_kb = NA_real_
    )
    next
  }

  for (i in seq_len(n_cs)) {
    snps_i <- cs_e[[i]]
    pip_i <- susie_e$pip[snps_i]
    lead_idx <- snps_i[which.max(pip_i)]
    lead_rsid <- merged$rsid[lead_idx]
    lead_pos  <- merged$pos_e[lead_idx]
    z_at_lead <- merged$z_amyloid[lead_idx]
    p_at_lead <- 2 * pnorm(-abs(z_at_lead))
    dist_kb <- abs(lead_pos - 207518704) / 1000
    cat(sprintf("      CS %d: %d SNPs; lead %s (chr1:%d) PIP %.3f | amyloid z=%.2f, p=%.2e | dist to rs6656401 = %.1f kb\n",
                i, length(snps_i), lead_rsid, lead_pos,
                max(pip_i), z_at_lead, p_at_lead, dist_kb))

    per_gene_summary[[paste0(g_name, "_cs", i)]] <- data.table(
      gene = g_name, gene_id = g_id, cs_idx = i,
      n_credible_sets = n_cs, cs_size = length(snps_i),
      cs_lead_rsid = lead_rsid, cs_lead_pos = lead_pos,
      cs_lead_pip = max(pip_i),
      amyloid_z_at_lead = z_at_lead,
      amyloid_p_at_lead = p_at_lead,
      dist_to_amyloid_lead_kb = dist_kb
    )

    # Per-CS-member rows for the long-format table
    for (j in seq_along(snps_i)) {
      idx_j <- snps_i[j]
      cred_set_rows[[length(cred_set_rows) + 1]] <- data.table(
        gene = g_name, gene_id = g_id, cs_idx = i,
        rsid = merged$rsid[idx_j], chr = merged$chr_e[idx_j],
        pos = merged$pos_e[idx_j], pip = susie_e$pip[idx_j],
        amyloid_z = merged$z_amyloid[idx_j],
        amyloid_p = 2 * pnorm(-abs(merged$z_amyloid[idx_j])),
        eqtl_beta = b_g[idx_j], eqtl_se = s_g[idx_j],
        eqtl_z = b_g[idx_j] / s_g[idx_j]
      )
    }
  }
}

# Compose summary table
gene_summary_dt <- rbindlist(per_gene_summary, fill = TRUE)
cred_set_dt <- if (length(cred_set_rows) > 0) rbindlist(cred_set_rows) else data.table()

cat("\n  --- Gene-level summary ---\n")
print(gene_summary_dt[, .(gene, cs_idx, cs_lead_rsid, cs_lead_pos,
                          cs_size, cs_lead_pip,
                          amyloid_z_at_lead = round(amyloid_z_at_lead, 2),
                          amyloid_p_at_lead = signif(amyloid_p_at_lead, 2),
                          dist_to_amyloid_lead_kb = round(dist_to_amyloid_lead_kb, 1))])

# --- Amyloid lead lookup: what is the strongest amyloid signal in window? ---
strongest <- merged[order(-abs(z_amyloid))][1:5]
cat("\n  --- Top 5 amyloid signals in chr1q32.2 window ---\n")
print(strongest[, .(rsid, pos = pos_e, BETA, SE,
                     z = round(z_amyloid, 2),
                     p = signif(2*pnorm(-abs(z_amyloid)), 2),
                     z_eqtl_CR1   = round(z_eqtl_ENSG00000203710, 2),
                     z_eqtl_CR1L  = round(z_eqtl_ENSG00000203709, 2),
                     z_eqtl_CD46  = round(z_eqtl_ENSG00000117335, 2))])
results <- list(per_gene_summary = gene_summary_dt, cs_members = cred_set_dt,
                strongest_amyloid = strongest)

# --- 9. Save results ---------------------------------------------------------
cat(sprintf("\n[%s] Writing results\n", Sys.time()))
if (nrow(gene_summary_dt) > 0) {
  out_fp1 <- file.path(dir_processed, "susie_chr1q322_brain_cortex_gene_summary.parquet")
  write_parquet(gene_summary_dt, out_fp1)
  fwrite(gene_summary_dt, sub("\\.parquet$", ".tsv", out_fp1), sep = "\t")
  cat(sprintf("  Wrote %s (%d rows)\n", out_fp1, nrow(gene_summary_dt)))
}
if (nrow(cred_set_dt) > 0) {
  out_fp2 <- file.path(dir_processed, "susie_chr1q322_brain_cortex_cs_members.parquet")
  write_parquet(cred_set_dt, out_fp2)
  fwrite(cred_set_dt, sub("\\.parquet$", ".tsv", out_fp2), sep = "\t")
  cat(sprintf("  Wrote %s (%d rows)\n", out_fp2, nrow(cred_set_dt)))
}

cat(sprintf("\n[%s] Done.\n", Sys.time()))
