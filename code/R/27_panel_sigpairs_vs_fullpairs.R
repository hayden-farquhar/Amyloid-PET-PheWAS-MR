# =============================================================================
# Script 27: significant-pairs vs full-pairs coloc on the PRE-SPECIFIED panel
# =============================================================================
# Reviewer-round-2 issue 3(a): the artefact (sig-pairs -> PP.H4 ~0.98) was first
# shown on the 16 screening hits, i.e. phenotypes selected for strong locus
# signal. Here we run the sig-pairs-vs-full-pairs contrast on the PRE-SPECIFIED
# 24-phenotype negative-control panel + 6 AD positive controls (R/23), an
# unselected set, and restrict the informative comparison to phenotypes with a
# genuine marginal signal at the CR1 locus (a phenotype with no locus signal
# cannot colocalise under any method, so its low PP.H4 is trivial).
#
# Sig-pairs eQTL = GTEx Brain_Cortex significant-pairs for CR1 (40 SNPs).
# Full-pairs eQTL = already in coloc_symmetry_panel.tsv (Brain_Cortex full-pairs).
#
# Output: data/processed/panel_sigpairs_vs_fullpairs.tsv
#         results/figures/fig_power_floor/panel_sigpairs_vs_fullpairs.{pdf,png}
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(coloc); library(ggplot2) })
proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
fg_cache <- file.path(proj, "data", "interim", "coloc_symmetry")
procdir <- file.path(proj, "data", "processed")
figdir <- file.path(proj, "results", "figures", "fig_power_floor")

# ---- CR1 significant-pairs eQTL (Brain_Cortex) ------------------------------
sp <- fread(cmd = paste0("gunzip -c ", shQuote(file.path(proj, "data", "reference",
  "gtex_v8", "brain", "eqtls", "Brain_Cortex.v8.EUR.signif_pairs.txt.gz")),
  " | awk -F'\\t' 'NR==1 || $1 ~ /ENSG00000203710/'"))
# variant_id chr1_pos_ref_alt_b38
sp[, c("c","pos","ref","alt") := tstrsplit(variant_id, "_", keep = 1:4)]
sp[, pos := as.integer(pos)]
sp[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
sp <- sp[!duplicated(key)]
cat(sprintf("CR1 significant-pairs SNPs: %d\n", nrow(sp)))
D_eqtl_sig <- function(keys) {
  d <- sp[key %in% keys][order(key)]
  list(beta = d$slope, varbeta = d$slope_se^2, snp = d$key, type = "quant",
       N = 205, MAF = d$maf)
}

full <- fread(file.path(procdir, "coloc_symmetry_panel.tsv"))
meta <- fread(file.path(proj, "osf", "phecode_filter_prespec.csv"))
setnames(meta, c("phenocode","phenotype","category","chapter_prefix",
                 "num_cases","num_controls","url"))
fg_cols <- c("chrom","pos","ref","alt","rsids","genes","pval","mlogp",
             "beta","sebeta","af_alt","af_cc","af_ctrl")

res <- list()
for (ph in full$phecode) {
  f <- file.path(fg_cache, sprintf("%s.tsv", ph))
  if (!file.exists(f)) next
  d <- fread(f, header = FALSE, col.names = fg_cols)
  d <- d[sebeta > 0 & is.finite(beta) & af_alt > 0.01 & af_alt < 0.99]
  d[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
  d <- d[!duplicated(key)]
  common <- intersect(d$key, sp$key)
  if (length(common) < 10) next
  d2 <- d[key %in% common][order(key)]
  mr <- meta[phenocode == ph]; N <- mr$num_cases + mr$num_controls
  D_pt <- list(beta = d2$beta, varbeta = d2$sebeta^2, snp = d2$key, type = "cc",
               N = N, s = mr$num_cases / N, MAF = pmin(d2$af_alt, 1 - d2$af_alt))
  cc <- tryCatch(suppressMessages(suppressWarnings(coloc.abf(D_pt, D_eqtl_sig(common)))),
                 error = function(e) NULL)
  if (is.null(cc)) next
  res[[ph]] <- data.table(phecode = ph, n_sig = length(common),
                          pph4_sigpairs = as.numeric(cc$summary["PP.H4.abf"]))
}
sig <- rbindlist(res)
out <- merge(full[, .(phecode, class, phenotype, min_p_locus,
                      pph4_fullpairs = pp_h4)], sig, by = "phecode")
out[, has_locus_signal := min_p_locus < 1e-4]
setorder(out, -pph4_sigpairs)
fwrite(out, file.path(procdir, "panel_sigpairs_vs_fullpairs.tsv"), sep = "\t")

cat("\n--- PRE-SPECIFIED panel: sig-pairs vs full-pairs (phenotypes WITH a genuine locus signal) ---\n")
print(out[has_locus_signal == TRUE,
          .(phecode, class, min_p = signif(min_p_locus, 2),
            sigpairs = round(pph4_sigpairs, 3), fullpairs = round(pph4_fullpairs, 3))])
sigset <- out[has_locus_signal == TRUE & class == "negative_control"]
cat(sprintf("\nNegative controls WITH locus signal (n=%d): median sig-pairs PP.H4 = %.2f -> median full-pairs PP.H4 = %.2f\n",
            nrow(sigset), median(sigset$pph4_sigpairs), median(sigset$pph4_fullpairs)))
cat(sprintf("  inflated by sig-pairs (PP.H4>0.5): %d/%d; corrected by full-pairs (PP.H4<0.5): %d/%d\n",
            sum(sigset$pph4_sigpairs > 0.5), nrow(sigset),
            sum(sigset$pph4_fullpairs < 0.5), nrow(sigset)))

# ---- figure -----------------------------------------------------------------
ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999")
m <- melt(out[has_locus_signal == TRUE], id.vars = c("phecode","class"),
          measure.vars = c("pph4_sigpairs","pph4_fullpairs"),
          variable.name = "input", value.name = "pph4")
m[, input := factor(input, levels = c("pph4_sigpairs","pph4_fullpairs"),
                    labels = c("significant-pairs","full-pairs"))]
p <- ggplot(m, aes(input, pph4, group = phecode, colour = class)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = ok$grey, linewidth = 0.3) +
  geom_line(linewidth = 0.5, alpha = 0.8) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(AD_positive_control = ok$blue,
                                 negative_control = ok$vermil),
                      labels = c("AD positive control", "negative control")) +
  ylim(0, 1) +
  labs(x = NULL, y = "PP.H4 (phenotype x CR1)",
       title = "Pre-specified panel with a genuine CR1-locus signal",
       caption = "Unselected set; significant-pairs inflates negative controls, full-pairs corrects") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.title = element_blank(), plot.caption = element_text(size = 6, colour = ok$grey))
ggsave(file.path(figdir, "panel_sigpairs_vs_fullpairs.pdf"), p, width = 5.0, height = 3.8)
ggsave(file.path(figdir, "panel_sigpairs_vs_fullpairs.png"), p, width = 5.0, height = 3.8, dpi = 300)
cat("\nWrote panel_sigpairs_vs_fullpairs.{tsv,pdf,png}\n")
