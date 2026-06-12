# =============================================================================
# Script 22: coloc.abf robustness at chr1q32.2 / CR1
# =============================================================================
#
# Reframe add-on C. Tests how fragile the headline PP.H4 = 0.983 (amyloid x CR1
# brain-cortex eQTL) is to two analyst choices the reviewer flagged (M4):
#   (i)  the coloc.abf prior p12 (default 1e-5), swept 1e-6 .. 1e-4;
#   (ii) the input truncation — GTEx ships SIGNIFICANT-pairs only; Phase 4 used
#        that truncated input. Here we run coloc.abf on FULL-pairs eQTL (the
#        unbiased input) and on a significant-pairs emulation, and compare.
#
# coloc.abf does not use LD; the contrast isolates the per-SNP Bayes-factor
# evidence and the prior. Both datasets are GRCh38, joined on chr:pos.
#
# Input:  data/interim/coloc_susie/ali_chr1q322.tsv          (amyloid region)
#         data/interim/coloc_susie/chr1q32.2_brain_cortex_raw.tsv (full-pairs)
# Output: data/processed/coloc_robustness_cr1.tsv
#         results/figures/fig_power_floor/coloc_robustness_cr1.{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(ggplot2)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
idir <- file.path(proj, "data", "interim", "coloc_susie")
figdir <- file.path(proj, "results", "figures", "fig_power_floor")
procdir <- file.path(proj, "data", "processed")

CR1 <- "ENSG00000203710"
N_AMY <- 13409; N_EQTL <- 205

# ---- amyloid region ---------------------------------------------------------
# Ali 2023 censors some allele frequencies as ">0.05" (common variants);
# coloc.abf uses sdY=1 for the amyloid dataset so MAF is not needed here.
amy <- fread(file.path(idir, "ali_chr1q322.tsv"))
amy <- amy[SE > 0 & is.finite(BETA)]
amy[, key := paste(BP, pmin(EFFECT_ALLELE, OTHER_ALLELE),
                   pmax(EFFECT_ALLELE, OTHER_ALLELE), sep = "_")]
amy <- amy[!duplicated(key)]

# ---- CR1 eQTL full-pairs ----------------------------------------------------
eq <- fread(file.path(idir, "chr1q32.2_brain_cortex_raw.tsv"), header = FALSE,
            col.names = c("variant_id","r2","pvalue","mtid","mtoid","maf",
                          "gene_id","median_tpm","beta","se","an","ac",
                          "chr_e","pos_e","ref_e","alt_e","type","rsid"))
eq <- eq[gene_id == CR1 & se > 0 & is.finite(beta) & maf > 0.01]
eq[, key := paste(pos_e, pmin(ref_e, alt_e), pmax(ref_e, alt_e), sep = "_")]
eq <- eq[!duplicated(key)]

common <- intersect(amy$key, eq$key)
cat(sprintf("Amyloid SNPs: %d  |  CR1 eQTL SNPs: %d  |  common: %d\n",
            nrow(amy), nrow(eq), length(common)))
amy2 <- amy[key %in% common][order(key)]
eq2  <- eq[key %in% common][order(key)]
stopifnot(all(amy2$key == eq2$key))

D_amy <- list(beta = amy2$BETA, varbeta = amy2$SE^2, snp = amy2$key,
              type = "quant", N = N_AMY, sdY = 1)
mk_eqtl <- function(d) list(beta = d$beta, varbeta = d$se^2, snp = d$key,
                            type = "quant", N = N_EQTL, MAF = d$maf)

run_pph4 <- function(eqtl_dt, p12) {
  suppressMessages(
    res <- coloc.abf(D_amy, mk_eqtl(eqtl_dt[key %in% D_amy$snp][order(key)]),
                     p1 = 1e-4, p2 = 1e-4, p12 = p12))
  as.numeric(res$summary["PP.H4.abf"])
}

# ---- (i) prior p12 sweep on full-pairs --------------------------------------
p12_grid <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
full_sweep <- data.table(input = "full_pairs", p12 = p12_grid,
  pph4 = sapply(p12_grid, function(p) run_pph4(eq2, p)))

# ---- (ii) significant-pairs emulation (truncate eQTL to strong signal) ------
# GTEx significant-pairs keeps gene-level FDR-significant SNP-gene pairs only;
# emulate by retaining eQTL SNPs at p < 1e-4 (strong-signal truncation).
eq_sig <- eq2[pvalue < 1e-4]
cat(sprintf("Significant-pairs emulation retains %d of %d common CR1 eQTL SNPs\n",
            nrow(eq_sig), nrow(eq2)))
sig_sweep <- data.table(input = "sig_pairs_emul", p12 = p12_grid,
  pph4 = sapply(p12_grid, function(p) run_pph4(eq_sig, p)))

out <- rbind(full_sweep, sig_sweep)
out[, phase4_reported := 0.983]
fwrite(out, file.path(procdir, "coloc_robustness_cr1.tsv"), sep = "\t")
cat("\n--- coloc.abf PP.H4 robustness (amyloid x CR1) ---\n")
print(out)

# ---- figure -----------------------------------------------------------------
ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999")
p <- ggplot(out, aes(p12, pph4, colour = input, group = input)) +
  geom_hline(yintercept = 0.983, linetype = "dashed",
             colour = ok$grey, linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  annotate("text", x = 1e-5, y = 0.985, label = "Phase-4 reported 0.983",
           size = 2.2, colour = ok$grey, vjust = -0.4) +
  scale_x_log10() + scale_colour_manual(
    values = c(full_pairs = ok$blue, sig_pairs_emul = ok$vermil),
    labels = c(full_pairs = "full-pairs eQTL",
               sig_pairs_emul = "significant-pairs (emulated)")) +
  ylim(0, 1) +
  labs(x = "coloc.abf prior p12 (default 1e-5)", y = "PP.H4 (amyloid x CR1)",
       title = "coloc.abf PP.H4 sensitivity at chr1q32.2 / CR1",
       caption = "Sensitivity to prior and to GTEx significant-pairs truncation") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = c(0.72, 0.22), legend.title = element_blank(),
        plot.caption = element_text(size = 6, colour = ok$grey))
ggsave(file.path(figdir, "coloc_robustness_cr1.pdf"), p, width = 5.2, height = 3.4)
ggsave(file.path(figdir, "coloc_robustness_cr1.png"), p, width = 5.2, height = 3.4,
       dpi = 300)
cat("\nWrote coloc_robustness_cr1.{tsv,pdf,png}\n")
