# =============================================================================
# Script 16: Manuscript figures
# =============================================================================
#
# Generates the three core manuscript figures from the atlas + replication
# data already on disk:
#
#   Figure 2 — CR1 locus zoom: amyloid-PET -log10(p) and CR1 Brain_Cortex
#              eQTL -log10(p) on the same x-axis (chr1:207 Mb ± 500 kb).
#              Visualises the PP.H4 = 0.98 colocalisation.
#
#   Figure 3 — Cross-ancestry concordance: Ali EUR β vs Ali ASN β scatter
#              with the 16 unified IVs, coloured by panel membership.
#              Visualises the within-study cross-ancestry replication.
#
#   Figure 4 — Forest plot of top robust hits: top 15 phecodes × IVW β with
#              95% CI bars, panel as colour, robust+coloc-strengthened flagged.
#
# Outputs PNG (300 DPI for inspection) + PDF (vector for publication) to
# results/figures/.
#
# Run as:
#   Rscript R/16_manuscript_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(R.utils)
})

proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_figs <- file.path(proj_root, "results", "figures")
dir.create(dir_figs, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Figure 2 — CR1 locus zoom
# =============================================================================

cat(sprintf("[%s] Figure 2 — CR1 locus zoom\n", Sys.time()))

# CR1 lead SNP: rs6656401 at chr1:207,518,704 (GRCh38)
LEAD_CHR <- 1L
LEAD_BP  <- 207518704L
WINDOW_KB <- 500L

# 1. Read Ali 2023 amyloid sumstats within window
ali_fp <- file.path(proj_root, "data/raw/ali_2023",
  "AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt")
cmd <- sprintf("awk -v c=%d -v l=%d -v h=%d 'NR==1 || ($1==c && $2>=l && $2<=h)' %s",
                LEAD_CHR, LEAD_BP - WINDOW_KB*1000L, LEAD_BP + WINDOW_KB*1000L,
                shQuote(ali_fp))
ali <- fread(cmd = cmd)
cat(sprintf("  Ali SNPs in window: %d\n", nrow(ali)))

# 2. Read GTEx Brain_Cortex CR1 (ENSG00000203710) eQTLs within window
gtex_fp <- file.path(proj_root, "data/reference/gtex_v8/brain/eqtls",
                      "Brain_Cortex.v8.EUR.signif_pairs.txt.gz")
gtex <- fread(gtex_fp)
# Filter to CR1 (gene-level; ENSG ID may have version suffix)
cr1_gtex <- gtex[grepl("^ENSG00000203710", phenotype_id)]
# Parse variant_id chr_pos_ref_alt_b38
cr1_gtex[, c("var_chr","var_pos","var_ref","var_alt","build") := tstrsplit(variant_id, "_")]
cr1_gtex[, var_pos := as.integer(var_pos)]
cr1_gtex <- cr1_gtex[var_chr == sprintf("chr%d", LEAD_CHR) &
                     var_pos >= LEAD_BP - WINDOW_KB*1000L &
                     var_pos <= LEAD_BP + WINDOW_KB*1000L]
cat(sprintf("  CR1 brain-cortex eQTL SNPs in window: %d\n", nrow(cr1_gtex)))

# Build long-format data for the dual plot
df_amyloid <- ali[, .(track = "Amyloid-PET (Ali 2023, N=13,409)",
                      pos = as.integer(BP),
                      neglogp = -log10(as.numeric(P)))]
df_eqtl    <- cr1_gtex[, .(track = "CR1 brain-cortex eQTL (GTEx v8, N=205)",
                           pos = var_pos,
                           neglogp = -log10(as.numeric(pval_nominal)))]
df <- rbindlist(list(df_amyloid, df_eqtl))
df <- df[!is.na(neglogp) & is.finite(neglogp)]
# Cap extreme p-values for visualisation
df[, neglogp := pmin(neglogp, 50)]

fig2 <- ggplot(df, aes(pos / 1e6, neglogp)) +
  geom_point(aes(colour = track), alpha = 0.6, size = 1) +
  geom_vline(xintercept = LEAD_BP / 1e6, linetype = "dashed", colour = "grey50") +
  annotate("text", x = LEAD_BP / 1e6 + 0.01, y = max(df$neglogp) * 0.95,
            label = "rs6656401\n(CR1 lead)", hjust = 0, size = 3, colour = "grey30") +
  facet_wrap(~ track, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#3366CC", "#CC3333"), guide = "none") +
  labs(
    title = sprintf("CR1 locus zoom (chr%d:%.1f–%.1f Mb)",
                     LEAD_CHR, (LEAD_BP - WINDOW_KB*1000L)/1e6, (LEAD_BP + WINDOW_KB*1000L)/1e6),
    subtitle = expression("Colocalisation PP.H"[4]*" = 0.98 between amyloid-PET signal and CR1 brain-cortex eQTL"),
    x = sprintf("Position (Mb, GRCh38, chr%d)", LEAD_CHR),
    y = expression(-log[10](p))
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(dir_figs, "fig2_CR1_locuszoom.png"), fig2, width = 8, height = 7, dpi = 300)
ggsave(file.path(dir_figs, "fig2_CR1_locuszoom.pdf"), fig2, width = 8, height = 7)
cat(sprintf("  Wrote fig2_CR1_locuszoom.{png,pdf}\n"))

# =============================================================================
# Figure 3 — Cross-ancestry concordance (Ali EUR vs Ali ASN)
# =============================================================================

cat(sprintf("\n[%s] Figure 3 — Cross-ancestry concordance\n", Sys.time()))

rep_iv <- as.data.table(read_parquet(
  file.path(proj_root, "data/processed/replication_ali_asn_per_iv.parquet")))
rep_iv <- rep_iv[!is.na(beta_eur) & !is.na(beta_asn)]
cat(sprintf("  IVs in concordance plot: %d\n", nrow(rep_iv)))

# Annotate panel membership (highest-rank panel for label colour)
rep_iv[, panel_label := fcase(
  grepl("tier1_apoe_wide_excl", panels, fixed = TRUE), "T1 wide-excl (Pivot A primary)",
  grepl("tier1_apoe_narrow_excl", panels, fixed = TRUE), "T1 narrow-excl",
  grepl("tier1_apoe_inclusive", panels, fixed = TRUE),  "T1 APOE-inclusive",
  grepl("tier2_apoe_wide_excl", panels, fixed = TRUE),  "T2 wide-excl",
  grepl("tier2_apoe_narrow_excl", panels, fixed = TRUE),"T2 narrow-excl",
  default = "other"
)]
rep_iv[, label := fifelse(h4_replicates | rsid == "rs6656401", rsid, NA_character_)]

fig3 <- ggplot(rep_iv, aes(beta_eur, beta_asn)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(colour = panel_label, size = -log10(pmax(replication_pval, 1e-10))),
              alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 25, box.padding = 0.5) +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](p[ASN]))) +
  scale_colour_brewer(palette = "Set2", name = "Highest-priority panel") +
  labs(
    title = "Cross-ancestry instrument concordance",
    subtitle = "Ali 2023 EUR effect vs Ali 2023 East Asian effect (within-study, harmonised on EUR effect-allele); dashed diagonal = perfect cross-ancestry agreement",
    x = expression(beta[EUR]*" — Ali 2023 multi-ethnic (N=11,556 NHW)"),
    y = expression(beta[EAS]*" — Ali 2023 East Asian arm (N=1,494)")
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"))
ggsave(file.path(dir_figs, "fig3_cross_ancestry_concordance.png"), fig3,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(dir_figs, "fig3_cross_ancestry_concordance.pdf"), fig3,
       width = 10, height = 8)
cat(sprintf("  Wrote fig3_cross_ancestry_concordance.{png,pdf}\n"))

# =============================================================================
# Figure 4 — Forest plot of top robust hits
# =============================================================================

cat(sprintf("\n[%s] Figure 4 — Forest plot of top robust hits\n", Sys.time()))

atlas <- as.data.table(read_parquet(
  file.path(proj_root, "data/processed/amyloid_pet_causal_atlas.parquet")))

# Filter to robust hits; deduplicate phecodes by keeping the strongest panel
robust <- atlas[robust_hit == TRUE & !is.na(pval_ivw)]
robust <- robust[order(pval_ivw)]
# Pick top 15 unique phecodes (keep all panels for each)
top_phecodes <- unique(robust$phecode)[1:min(15, uniqueN(robust$phecode))]
top <- robust[phecode %in% top_phecodes]
top[, panel_label := fcase(
  panel == "tier1_apoe_inclusive",   "T1 APOE-inclusive",
  panel == "tier1_apoe_narrow_excl", "T1 narrow-excl",
  panel == "tier1_apoe_wide_excl",   "T1 wide-excl (Pivot A primary)",
  panel == "tier2_apoe_narrow_excl", "T2 narrow-excl",
  panel == "tier2_apoe_wide_excl",   "T2 wide-excl"
)]
top[, phecode_label := paste0(substr(phenotype, 1, 35),
                               ifelse(nchar(phenotype) > 35, "...", ""),
                               " (", phecode, ")")]
# Order phecodes by min-p across panels (most significant on top)
phecode_order <- top[, .(min_p = min(pval_ivw, na.rm = TRUE)),
                     by = phecode_label][order(-min_p), phecode_label]
top[, phecode_label := factor(phecode_label, levels = phecode_order)]

top[, ci_lo := b_ivw - 1.96 * se_ivw]
top[, ci_hi := b_ivw + 1.96 * se_ivw]

fig4 <- ggplot(top, aes(b_ivw, phecode_label, colour = panel_label)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                  position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = coloc_strengthened), size = 2.5,
              position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                      name = "Coloc-strengthened\n(PP.H4 > 0.5)") +
  scale_colour_brewer(palette = "Set1", name = "Instrument panel") +
  labs(
    title = "Top robust hits across instrument panels",
    subtitle = "β_IVW with 95% CI per (phecode × panel). Robust hits pass: IVW p<1e-4, weighted-median direction match, clean Egger intercept (p>0.05), clean Cochran Q (p>0.05). Triangles indicate brain-eQTL colocalisation evidence (PP.H4 > 0.5).",
    x = expression(beta[IVW]*" (log-odds per 1 SD genetically-predicted amyloid-PET)"),
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text.y = element_text(size = 8))
ggsave(file.path(dir_figs, "fig4_forest_top_hits.png"), fig4,
       width = 12, height = 8, dpi = 300)
ggsave(file.path(dir_figs, "fig4_forest_top_hits.pdf"), fig4,
       width = 12, height = 8)
cat(sprintf("  Wrote fig4_forest_top_hits.{png,pdf}\n"))

cat(sprintf("\n[%s] All 3 manuscript figures generated successfully.\n", Sys.time()))
