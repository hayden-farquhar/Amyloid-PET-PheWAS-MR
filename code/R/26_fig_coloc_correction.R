# =============================================================================
# Script 26: Figure — significant-pairs vs full-pairs coloc correction
# =============================================================================
# Slope/dumbbell figure: best PP.H4 per coloc-tested phecode under the Phase-4
# significant-pairs input (all ~0.98-1.0, indiscriminate) vs the corrected
# full-pairs input (discriminates AD-spectrum CR1 coloc from false positives).
# Source: data/processed/coloc_sigpairs_vs_fullpairs.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
figdir <- file.path(proj, "results", "figures", "fig_power_floor")

d <- fread(file.path(proj, "data", "processed", "coloc_sigpairs_vs_fullpairs.tsv"))
ad <- "DEMENTIA|ALZ|^AD_|NEURODEG"
d[, cls := fifelse(grepl(ad, phecode) & !grepl("EO_EXMORE|VASCDEM", phecode),
                   "late-onset AD-spectrum (CR1 coloc)",
                   "negative / non-amyloid (false positive)")]
m <- melt(d, id.vars = c("phecode","cls"),
          measure.vars = c("pph4_sigpairs","pph4_fullpairs"),
          variable.name = "input", value.name = "pph4")
m[, input := factor(input, levels = c("pph4_sigpairs","pph4_fullpairs"),
                    labels = c("significant-pairs\n(Phase 4)", "full-pairs\n(corrected)"))]
m[, phecode := factor(phecode, levels = d[order(pph4_fullpairs)]$phecode)]

ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999")
p <- ggplot(m, aes(input, pph4, group = phecode, colour = cls)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = ok$grey, linewidth = 0.3) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("late-onset AD-spectrum (CR1 coloc)" = ok$blue,
                                 "negative / non-amyloid (false positive)" = ok$vermil)) +
  ylim(0, 1) +
  labs(x = NULL, y = "best PP.H4 across brain tissues",
       title = "Full-pairs eQTL restores colocalisation specificity at CR1 / APOE",
       caption = "Each line a coloc-tested screening hit; significant-pairs + across-tissue max inflates all to ~0.98") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 6),
        plot.caption = element_text(size = 6, colour = ok$grey))
ggsave(file.path(figdir, "coloc_sigpairs_vs_fullpairs.pdf"), p, width = 5.2, height = 4.0)
ggsave(file.path(figdir, "coloc_sigpairs_vs_fullpairs.png"), p, width = 5.2, height = 4.0, dpi = 300)
cat("Wrote coloc_sigpairs_vs_fullpairs.{pdf,png}\n")
cat(sprintf("Late-onset AD-spectrum median full-pairs PP.H4: %.3f\n",
            median(d[grepl(ad, phecode) & !grepl("EO_EXMORE|VASCDEM", phecode)]$pph4_fullpairs)))
cat(sprintf("False-positive median full-pairs PP.H4: %.3f\n",
            median(d[!grepl(ad, phecode) | grepl("EO_EXMORE|VASCDEM", phecode)]$pph4_fullpairs)))
