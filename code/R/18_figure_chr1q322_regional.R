#!/usr/bin/env Rscript
# =============================================================================
# Script 18: Figure 4 — chr1q32.2 regional plot showing amyloid GWAS signal
#            alongside per-gene SuSiE credible sets
# =============================================================================
#
# Two-panel plot:
#   Top:    amyloid-PET GWAS Manhattan track (-log10 p) for chr1q32.2 ± 500 kb
#           with the lead SNP (rs6656401) annotated
#   Bottom: per-gene SuSiE credible-set member SNPs (CR1, CD46) coloured by gene,
#           sized by PIP, with credible-set leads (rs679515, rs2724384) annotated
#           CR1L is noted as having no credible set
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

proj_root <- "/Users/haydenfarquhar/Documents/Research Projects/81 Amyloid PET PheWAS MR"
ali_fp <- file.path(proj_root, "data/interim/coloc_susie/ali_chr1q322.tsv")
cs_fp  <- file.path(proj_root, "data/processed/susie_chr1q322_brain_cortex_cs_members.parquet")
out_dir <- file.path(proj_root, "results/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Anchor coordinates
AMYLOID_LEAD <- 207518704  # rs6656401 GRCh38
CR1_CS_LEAD  <- 207577223  # rs679515
CD46_CS_LEAD <- 207756858  # rs2724384

# --- Load data ---------------------------------------------------------------
ali <- fread(ali_fp)
ali[, neglog10p := -log10(P)]
ali[, kb_from_lead := (BP - AMYLOID_LEAD) / 1000]

cs <- as.data.table(arrow::read_parquet(cs_fp))
cs[, neglog10p_amyloid := -log10(amyloid_p)]
cs[, kb_from_lead := (pos - AMYLOID_LEAD) / 1000]

# Filter to plot window: ± 350 kb around the amyloid lead (focused view)
WIN_KB <- 350
ali_plot <- ali[BP >= AMYLOID_LEAD - WIN_KB*1000 & BP <= AMYLOID_LEAD + WIN_KB*1000]
cs_plot  <- cs[pos >= AMYLOID_LEAD - WIN_KB*1000 & pos <= AMYLOID_LEAD + WIN_KB*1000]

# --- Annotation points -------------------------------------------------------
annot <- rbind(
  data.table(rsid = "rs6656401", pos = AMYLOID_LEAD, kb_from_lead = 0,
             label = "rs6656401\n(amyloid GWAS lead)",
             colour = "black"),
  data.table(rsid = "rs679515", pos = CR1_CS_LEAD, kb_from_lead = (CR1_CS_LEAD - AMYLOID_LEAD)/1000,
             label = "rs679515\n(CR1 CS lead)",
             colour = "#d62728"),
  data.table(rsid = "rs2724384", pos = CD46_CS_LEAD, kb_from_lead = (CD46_CS_LEAD - AMYLOID_LEAD)/1000,
             label = "rs2724384\n(CD46 CS lead)",
             colour = "#1f77b4")
)
annot[, neglog10p_at_pos := ali_plot$neglog10p[match(pos, ali_plot$BP)]]
# For SNPs not in Ali data exactly, find nearest
for (i in seq_len(nrow(annot))) {
  if (is.na(annot$neglog10p_at_pos[i])) {
    nearest_idx <- which.min(abs(ali_plot$BP - annot$pos[i]))
    annot$neglog10p_at_pos[i] <- ali_plot$neglog10p[nearest_idx]
  }
}

# --- Top panel: amyloid GWAS Manhattan ---------------------------------------
p_top <- ggplot(ali_plot, aes(x = kb_from_lead, y = neglog10p)) +
  geom_point(alpha = 0.4, size = 1.0, colour = "grey50") +
  # Highlight amyloid signal envelope (top 5 SNPs)
  geom_point(data = ali_plot[neglog10p >= -log10(1e-9)],
             colour = "darkorange", size = 1.6, alpha = 0.8) +
  # Vertical lines at the three leads
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.5) +
  geom_vline(xintercept = (CR1_CS_LEAD - AMYLOID_LEAD)/1000,
             linetype = "dashed", colour = "#d62728", alpha = 0.5) +
  geom_vline(xintercept = (CD46_CS_LEAD - AMYLOID_LEAD)/1000,
             linetype = "dashed", colour = "#1f77b4", alpha = 0.5) +
  # Annotations
  geom_label_repel(
    data = annot, aes(x = kb_from_lead, y = neglog10p_at_pos, label = label),
    colour = annot$colour, fontface = "bold", size = 3.2,
    min.segment.length = 0, segment.size = 0.4,
    box.padding = 0.6, point.padding = 0.5,
    nudge_y = c(2, -1.5, 2), seed = 81
  ) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dotted", colour = "red", alpha = 0.6) +
  annotate("text", x = -WIN_KB*0.95, y = -log10(5e-8) + 0.3,
           label = "GWS p = 5e-8", colour = "red", size = 2.8, hjust = 0) +
  scale_x_continuous(limits = c(-WIN_KB, WIN_KB),
                     breaks = seq(-WIN_KB, WIN_KB, 100)) +
  labs(
    x = NULL,
    y = expression(-log[10](italic(p))~"amyloid-PET GWAS"),
    title = "chr1q32.2 region: amyloid-PET GWAS signal (Ali 2023 multi-ethnic)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# --- Bottom panel: SuSiE credible sets ---------------------------------------
# Color and shape mapping per gene
cs_plot[, gene_pretty := factor(gene, levels = c("CR1", "CR1L", "CD46"))]

p_bot <- ggplot(cs_plot, aes(x = kb_from_lead, y = gene_pretty,
                              size = pip, colour = gene_pretty)) +
  geom_point(alpha = 0.85) +
  # Mark the credible-set lead SNPs with a black ring
  geom_point(data = cs_plot[(gene == "CR1" & rsid == "rs679515") |
                              (gene == "CD46" & rsid == "rs2724384")],
             shape = 21, colour = "black", fill = NA, stroke = 1.2,
             size = 5.5) +
  # Verticals matching top panel
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.5) +
  geom_vline(xintercept = (CR1_CS_LEAD - AMYLOID_LEAD)/1000,
             linetype = "dashed", colour = "#d62728", alpha = 0.5) +
  geom_vline(xintercept = (CD46_CS_LEAD - AMYLOID_LEAD)/1000,
             linetype = "dashed", colour = "#1f77b4", alpha = 0.5) +
  # Add a "no credible set" annotation for CR1L
  annotate("text", x = 0, y = 2, label = "(no credible set in Brain_Cortex)",
           colour = "grey40", size = 3.0, fontface = "italic") +
  scale_x_continuous(limits = c(-WIN_KB, WIN_KB),
                     breaks = seq(-WIN_KB, WIN_KB, 100)) +
  scale_colour_manual(values = c(CR1 = "#d62728", CR1L = "grey60", CD46 = "#1f77b4"),
                      guide = "none") +
  scale_size_continuous(range = c(2, 6), name = "PIP",
                        breaks = c(0.05, 0.10, 0.15, 0.20)) +
  labs(
    x = "Distance from amyloid GWAS lead (rs6656401, chr1:207,518,704) — kb",
    y = "SuSiE credible-set\nmember per gene",
    title = "Per-gene SuSiE credible-set members (GTEx_V8 Brain_Cortex)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 11))

# --- Combine panels ----------------------------------------------------------
fig4 <- p_top / p_bot + plot_layout(heights = c(2, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = NULL,  # Caption goes under the figure in the manuscript; do not bake a "Figure N" title into the PNG (avoids submission-vs-figure label drift)
    caption = sprintf("Ali 2023 amyloid-PET multi-ethnic GWAS in the chr1q32.2 region (top panel; n=%d SNPs in ±%d kb window).\nPer-gene SuSiE credible-set member SNPs from GTEx_V8 Brain_Cortex (bottom panel; CR1 in red, CD46 in blue; CR1L has no credible set in this tissue).\nPoint size in panel B is proportional to SuSiE posterior inclusion probability (PIP). Black-ringed points are the credible-set leads (highest PIP per credible set).\nDashed vertical lines mark the amyloid GWAS lead rs6656401 (black), CR1 credible-set lead rs679515 (red), and CD46 credible-set lead rs2724384 (blue).\nThe CR1 credible set co-localises spatially with the amyloid signal envelope; the CD46 credible set sits 238 kb distant from the amyloid lead with no amyloid signal.",
                       nrow(ali_plot), WIN_KB),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.caption = element_text(size = 8, hjust = 0,
                                              margin = margin(t = 10)))
  )

# --- Save --------------------------------------------------------------------
png_fp <- file.path(out_dir, "figure_4_chr1q322_regional.png")
pdf_fp <- file.path(out_dir, "figure_4_chr1q322_regional.pdf")
ggsave(png_fp, fig4, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(pdf_fp, fig4, width = 10, height = 8, bg = "white")

cat(sprintf("Wrote %s\n", png_fp))
cat(sprintf("Wrote %s\n", pdf_fp))
cat(sprintf("Panel A: %d amyloid GWAS SNPs in ±%d kb window\n", nrow(ali_plot), WIN_KB))
cat(sprintf("Panel B: %d credible-set member SNPs across %d genes\n",
            nrow(cs_plot), length(unique(cs_plot$gene))))
