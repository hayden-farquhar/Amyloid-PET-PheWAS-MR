# =============================================================================
# Script 25: Reverse-direction MR — AD liability -> amyloid-PET burden
# =============================================================================
#
# Reframe add-on D. Single-SNP Steiger is a weak reverse-causation guard (M2).
# Here we run a genome-wide reverse MR: instrument AD liability with independent
# non-APOE genome-wide-significant loci (Bellenguez 2022) and estimate their
# effect on amyloid-PET burden (Ali 2023). A null AD->amyloid effect supports
# amyloid sitting UPSTREAM of (not downstream of) AD liability; a positive
# effect would indicate bidirectionality.
#
# This is a sensitivity, not a pillar: it asks about whole-AD-liability reverse
# causation, which is related to but not identical with the CR1-specific
# mediation-vs-direct-effect question (that one locus cannot resolve).
#
# Exposure: Bellenguez 2022 AD GWAS, genome-wide-significant SNPs, APOE-excluded,
#           LD-clumped to independence (PLINK, 1KG EUR panel).
# Outcome:  Ali 2023 multi-ethnic amyloid-PET (PRS-CS sumstats; SE from P).
#
# Output: data/processed/reverse_mr_ad_to_amyloid.tsv
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })

proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
procdir <- file.path(proj, "data", "processed")
ref <- file.path(proj, "..", "60 Statin MR Pathogen", "data", "reference", "EUR")
ali_fp <- file.path(proj, "data", "interim", "ali2023_prscs_sumstats.tsv")
tmp <- file.path(proj, "data", "interim", "reverse_mr"); dir.create(tmp, showWarnings = FALSE)

# ---- exposure: Bellenguez GWS, APOE-excluded -------------------------------
bel <- fread(file.path(proj, "data", "reference", "bellenguez_2022",
                       "bellenguez_gws.tsv"))
setnames(bel, c("variant_id","p_value","chromosome","base_pair_location",
                "effect_allele","other_allele","effect_allele_frequency",
                "beta","standard_error"),
         c("SNP","P","chr","pos","ea","oa","eaf","beta","se"), skip_absent = TRUE)
# exclude wide APOE block (chr19:44.4-46.5 Mb GRCh38)
bel <- bel[!(chr == 19 & pos >= 44400000 & pos <= 46500000)]
bel <- bel[grepl("^rs", SNP) & se > 0 & is.finite(beta)]
fwrite(bel[, .(SNP, P)], file.path(tmp, "ad_assoc.tsv"), sep = "\t")

# ---- clump to independent loci ---------------------------------------------
system2("plink", c("--bfile", shQuote(ref), "--clump",
        shQuote(file.path(tmp, "ad_assoc.tsv")),
        "--clump-p1", "5e-8", "--clump-p2", "5e-8", "--clump-r2", "0.001",
        "--clump-kb", "10000", "--clump-snp-field", "SNP", "--clump-field", "P",
        "--out", shQuote(file.path(tmp, "ad_clumped"))),
        stdout = file.path(tmp, "clump.log"), stderr = file.path(tmp, "clump.log"))
cl <- fread(file.path(tmp, "ad_clumped.clumped"))
lead <- cl$SNP
cat(sprintf("AD independent non-APOE instruments after clumping: %d\n", length(lead)))

# ---- outcome: amyloid-PET (Ali) lookup for lead SNPs ------------------------
writeLines(lead, file.path(tmp, "lead_rsids.txt"))
system2("grep", c("-wF", "-f", shQuote(file.path(tmp, "lead_rsids.txt")),
        shQuote(ali_fp)), stdout = file.path(tmp, "ali_sub.tsv"))
ali <- fread(file.path(tmp, "ali_sub.tsv"),
             header = FALSE, col.names = c("SNP","A1","A2","BETA","P"))
ali <- ali[SNP %in% lead]
# SE from two-sided P; guard P==0
ali[, P := pmax(P, 1e-300)]
ali[, SE := abs(BETA) / qnorm(P / 2, lower.tail = FALSE)]
ali <- ali[is.finite(SE) & SE > 0]
cat(sprintf("AD instruments found in Ali amyloid sumstats: %d\n", nrow(ali)))

# ---- harmonise (AD exposure, amyloid outcome) -------------------------------
exp_dat <- format_data(
  as.data.frame(bel[SNP %in% ali$SNP]), type = "exposure",
  snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "ea",
  other_allele_col = "oa", eaf_col = "eaf", pval_col = "P")
out_dat <- format_data(
  as.data.frame(ali), type = "outcome", snp_col = "SNP", beta_col = "BETA",
  se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2",
  pval_col = "P")
h <- harmonise_data(exp_dat, out_dat, action = 2)
h <- h[h$mr_keep, ]
cat(sprintf("Harmonised instruments retained: %d\n", nrow(h)))

# ---- MR: AD -> amyloid ------------------------------------------------------
mr_res <- mr(h, method_list = c("mr_ivw", "mr_weighted_median",
                                "mr_egger_regression"))
egg <- mr_pleiotropy_test(h)
out <- as.data.table(mr_res)[, .(method, nsnp, b, se, pval)]
out[, `:=`(or = exp(b), ci_l = b - 1.96 * se, ci_u = b + 1.96 * se)]
fwrite(out, file.path(procdir, "reverse_mr_ad_to_amyloid.tsv"), sep = "\t")

cat("\n--- Reverse MR: AD liability -> amyloid-PET (Ali) ---\n")
print(out)
cat(sprintf("\nMR-Egger intercept: %.4f (SE %.4f, p = %.3f)\n",
            egg$egger_intercept, egg$se, egg$pval))
cat(sprintf("\nIVW: beta = %.4f (95%% CI %.4f, %.4f), p = %.3g  [%d AD instruments]\n",
            out[method == "Inverse variance weighted", b],
            out[method == "Inverse variance weighted", ci_l],
            out[method == "Inverse variance weighted", ci_u],
            out[method == "Inverse variance weighted", pval],
            out[method == "Inverse variance weighted", nsnp]))
