# =============================================================================
# Script 30: Combined colocalisation-correction figure for EJHG (4-figure cap)
# =============================================================================
# Merges the two coloc-correction panels into a single 2-panel figure so the
# EJHG package fits the 4-figure limit:
#   A: best PP.H4 per coloc-tested screening hit, sig-pairs vs full-pairs (R/26)
#   B: pre-specified negative-control panel (locus-signal subset), sig vs full (R/27)
# Output: results/figures/fig_ejhg/fig3_coloc_correction_combined.{pdf,png}
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
procdir <- file.path(proj, "data", "processed")
figdir <- file.path(proj, "results", "figures", "fig_ejhg"); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999")
th <- theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.title = element_blank(), legend.text = element_text(size = 6),
        plot.caption = element_text(size = 5.5, colour = ok$grey))

# Panel A: 16 screening hits sig vs full
d <- fread(file.path(procdir, "coloc_sigpairs_vs_fullpairs.tsv"))
ad <- "DEMENTIA|ALZ|^AD_|NEURODEG"
d[, cls := fifelse(grepl(ad, phecode) & !grepl("EO_EXMORE|VASCDEM", phecode),
                   "late-onset AD-spectrum (CR1)", "negative / non-amyloid")]
mA <- melt(d, id.vars = c("phecode","cls"), measure.vars = c("pph4_sigpairs","pph4_fullpairs"),
           variable.name = "input", value.name = "pph4")
mA[, input := factor(input, levels = c("pph4_sigpairs","pph4_fullpairs"),
                     labels = c("sig-pairs\n(+max)", "full-pairs"))]
pA <- ggplot(mA, aes(input, pph4, group = phecode, colour = cls)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = ok$grey, linewidth = 0.3) +
  geom_line(linewidth = 0.5, alpha = 0.8) + geom_point(size = 1.5) +
  scale_colour_manual(values = c("late-onset AD-spectrum (CR1)" = ok$blue, "negative / non-amyloid" = ok$vermil)) +
  ylim(0, 1) + labs(x = NULL, y = "best PP.H4", title = "A  Coloc-tested screening hits") + th

# Panel B: pre-specified panel (locus-signal subset), Brain_Cortex sig vs full
p <- fread(file.path(procdir, "panel_sigpairs_vs_fullpairs.tsv"))
p <- p[has_locus_signal == TRUE]
mB <- melt(p, id.vars = c("phecode","class"), measure.vars = c("pph4_sigpairs","pph4_fullpairs"),
           variable.name = "input", value.name = "pph4")
mB[, input := factor(input, levels = c("pph4_sigpairs","pph4_fullpairs"),
                     labels = c("sig-pairs", "full-pairs"))]
pB <- ggplot(mB, aes(input, pph4, group = phecode, colour = class)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = ok$grey, linewidth = 0.3) +
  geom_line(linewidth = 0.5, alpha = 0.8) + geom_point(size = 1.5) +
  scale_colour_manual(values = c(AD_positive_control = ok$blue, negative_control = ok$vermil),
                      labels = c("AD positive control", "negative control")) +
  ylim(0, 1) + labs(x = NULL, y = "PP.H4 (Brain_Cortex)",
                    title = "B  Pre-specified panel (locus-signal subset)") + th

combo <- pA + pB + plot_layout(widths = c(1.15, 1))
ggsave(file.path(figdir, "fig3_coloc_correction_combined.pdf"), combo, width = 7.0, height = 3.6)
ggsave(file.path(figdir, "fig3_coloc_correction_combined.png"), combo, width = 7.0, height = 3.6, dpi = 300)
cat("Wrote fig3_coloc_correction_combined.{pdf,png}\n")
