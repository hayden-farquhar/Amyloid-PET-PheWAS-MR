# =============================================================================
# Script 23: Negative-control panel + coloc.abf evidence-symmetry sweep at CR1
# =============================================================================
#
# Reframe add-on B (novelty lift). The pre-registered single negative control
# (forearm fracture) returned the same PP.H4 = 0.983 at chr1q32.2 as the AD
# signal. This script generalises that one observation into an empirical
# characterisation: a PRE-SPECIFIED panel of biologically-implausible
# phenotypes (fractures, hernias, refraction/eye, dental, varicose veins,
# nasal septum, ingrowing nail) is run through coloc.abf against the CR1
# brain-cortex eQTL, and PP.H4 is plotted against each phenotype's marginal
# association strength at the locus.
#
# Claim tested: at this complement-dense, multi-eQTL, high-LD locus coloc.abf
# PP.H4 tracks the PRESENCE of any locus association in LD with the eQTL, not
# trait-specific shared causation. If implausible phenotypes with a locus
# signal return high PP.H4 identical to AD, PP.H4 is not trait-specific.
#
# Pulls FinnGen R12 regional sumstats (chr1:207.0-208.0 Mb) per phenotype via
# remote tabix, then coloc.abf vs CR1 full-pairs eQTL (default prior).
#
# Output: data/processed/coloc_symmetry_panel.tsv
#         results/figures/fig_power_floor/coloc_symmetry_panel.{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(coloc); library(ggplot2)
})

proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
idir <- file.path(proj, "data", "interim", "coloc_susie")
cache <- file.path(proj, "data", "interim", "coloc_symmetry")
figdir <- file.path(proj, "results", "figures", "fig_power_floor")
procdir <- file.path(proj, "data", "processed")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

# region: rs6656401 chr1:207,518,704 +/- 500 kb (GRCh38)
REGION <- "1:207018704-208018704"
CR1 <- "ENSG00000203710"; N_EQTL <- 205

# ---- panel definition -------------------------------------------------------
neg <- c("ST19_FRACT_FOREA","ST19_FRACT_LOWER_LEG_INCLU_ANKLE",
         "ST19_FRACT_SHOUL_UPPER_ARM","ST19_FRACT_WRIST_HAND_LEVEL",
         "ST19_FRACT_RIBS_STERNUM_THORACIC_SPINE","ST19_FRACT_FEMUR",
         "ST19_FRACT_FOOT_ANKLE","ST19_FRACT_SKULL_FACIAL_BONES",
         "ST19_FRACT_LUMBAR_SPINE_PELVIS",
         "K11_HERNIA","K11_HERING","ABDOM_HERNIA","K11_UMBHER","K11_VENTHER",
         "K11_DIAHER","H7_CATARACTSENILE","H7_MYOPIA","H7_ASTIGMATISM",
         "H7_REFRAACCOMMODIS","I9_VARICVE","J10_SEPTDEV","L12_NAIL_INGROW",
         "M13_RECUDISLOCATIO","K11_ANOMALI_DENTAL_ARCH_RELATIONS1")
pos <- c("F5_DEMENTIA_INCLAVO","F5_DEMENTIA","G6_AD_WIDE","G6_ALZHEIMER",
         "F5_ALZHDEMENT","KRA_PSY_DEMENTIA_EXMORE")
panel <- data.table(phecode = c(neg, pos),
                    class = c(rep("negative_control", length(neg)),
                              rep("AD_positive_control", length(pos))))

meta <- fread(file.path(proj, "osf", "phecode_filter_prespec.csv"))
setnames(meta, c("phenocode","phenotype","category","chapter_prefix",
                 "num_cases","num_controls","url"))
panel <- merge(panel, meta[, .(phecode = phenocode, phenotype, num_cases,
                               num_controls, url)], by = "phecode", all.x = TRUE)

# ---- CR1 eQTL full-pairs (reused) -------------------------------------------
eq <- fread(file.path(idir, "chr1q32.2_brain_cortex_raw.tsv"), header = FALSE,
            col.names = c("variant_id","r2","pvalue","mtid","mtoid","maf",
                          "gene_id","median_tpm","beta","se","an","ac",
                          "chr_e","pos_e","ref_e","alt_e","type","rsid"))
eq <- eq[gene_id == CR1 & se > 0 & is.finite(beta) & maf > 0.01]
eq[, key := paste(pos_e, pmin(ref_e, alt_e), pmax(ref_e, alt_e), sep = "_")]
eq <- eq[!duplicated(key)]
D_eqtl_full <- function(keys) {
  d <- eq[key %in% keys][order(key)]
  list(beta = d$beta, varbeta = d$se^2, snp = d$key, type = "quant",
       N = N_EQTL, MAF = d$maf)
}

# ---- tabix pulls (cache .tbi in a dedicated dir to avoid root spam) ----------
oldwd <- getwd(); setwd(cache)
on.exit({ unlink(file.path(cache, "*.tbi")); setwd(oldwd) }, add = TRUE)
pull <- function(url, phecode) {
  f <- file.path(cache, paste0(phecode, ".tsv"))
  if (file.exists(f) && file.info(f)$size > 0) return(f)
  ok <- tryCatch(system2("tabix", c(shQuote(url), REGION), stdout = f),
                 error = function(e) 1L)
  f
}

# ---- per-phenotype coloc.abf vs CR1 -----------------------------------------
res <- vector("list", nrow(panel))
for (i in seq_len(nrow(panel))) {
  ph <- panel[i]
  f <- pull(ph$url, ph$phecode)
  if (!file.exists(f) || file.info(f)$size == 0) { next }
  d <- fread(f, header = FALSE, col.names = c("chrom","pos","ref","alt","rsids",
             "genes","pval","mlogp","beta","sebeta","af_alt","af_cc","af_ctrl"))
  d <- d[sebeta > 0 & is.finite(beta) & af_alt > 0.01 & af_alt < 0.99]
  d[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
  d <- d[!duplicated(key)]
  common <- intersect(d$key, eq$key)
  if (length(common) < 50) next
  d2 <- d[key %in% common][order(key)]
  N <- ph$num_cases + ph$num_controls; s <- ph$num_cases / N
  D_pt <- list(beta = d2$beta, varbeta = d2$sebeta^2, snp = d2$key,
               type = "cc", N = N, s = s, MAF = pmin(d2$af_alt, 1 - d2$af_alt))
  cc <- suppressMessages(suppressWarnings(
    coloc.abf(D_pt, D_eqtl_full(common))))
  res[[i]] <- data.table(
    phecode = ph$phecode, class = ph$class, phenotype = ph$phenotype,
    num_cases = ph$num_cases, n_common = length(common),
    min_p_locus = min(d2$pval, na.rm = TRUE),
    pp_h4 = as.numeric(cc$summary["PP.H4.abf"]),
    pp_h3 = as.numeric(cc$summary["PP.H3.abf"]),
    pp_h0_h1_h2 = sum(as.numeric(cc$summary[c("PP.H0.abf","PP.H1.abf",
                                              "PP.H2.abf")])))
}
out <- rbindlist(res)
setorder(out, -pp_h4)
fwrite(out, file.path(procdir, "coloc_symmetry_panel.tsv"), sep = "\t")

cat("\n--- coloc.abf symmetry sweep at CR1 (n =", nrow(out), "phenotypes) ---\n")
print(out[, .(phecode, class, min_p_locus = signif(min_p_locus, 2),
              pp_h4 = round(pp_h4, 3))])
sig <- out[min_p_locus < 1e-4]
cat(sprintf("\nPhenotypes with a locus signal (min-p < 1e-4): %d\n", nrow(sig)))
cat(sprintf("  of these, PP.H4 > 0.8: %d (%.0f%%); PP.H4 > 0.5: %d\n",
            sum(sig$pp_h4 > 0.8), 100*mean(sig$pp_h4 > 0.8),
            sum(sig$pp_h4 > 0.5)))
negsig <- sig[class == "negative_control"]
cat(sprintf("  negative-control phenotypes with locus signal: %d; median PP.H4 = %.3f\n",
            nrow(negsig), median(negsig$pp_h4)))

# ---- figure: PP.H4 vs locus signal strength ---------------------------------
ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999")
out[, neglogp := -log10(pmax(min_p_locus, 1e-300))]
p <- ggplot(out, aes(neglogp, pp_h4, colour = class)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", colour = ok$grey,
             linewidth = 0.3) +
  geom_vline(xintercept = 4, linetype = "dotted", colour = ok$grey,
             linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_colour_manual(values = c(AD_positive_control = ok$blue,
                                 negative_control = ok$vermil),
                      labels = c("AD positive control",
                                 "negative control (implausible)")) +
  labs(x = expression(Phenotype~locus~signal~~-log[10](min~p)),
       y = "PP.H4 (phenotype x CR1 eQTL)",
       title = "coloc.abf evidence symmetry at chr1q32.2 / CR1",
       caption = "Each point a phenotype; PP.H4 tracks locus signal, not biological plausibility") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = c(0.72, 0.18), legend.title = element_blank(),
        plot.caption = element_text(size = 6, colour = ok$grey))
ggsave(file.path(figdir, "coloc_symmetry_panel.pdf"), p, width = 5.4, height = 3.6)
ggsave(file.path(figdir, "coloc_symmetry_panel.png"), p, width = 5.4, height = 3.6,
       dpi = 300)
cat("\nWrote coloc_symmetry_panel.{tsv,pdf,png}\n")
