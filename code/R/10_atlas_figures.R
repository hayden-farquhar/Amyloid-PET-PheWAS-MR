# =============================================================================
# Script 10: Amyloid-PET Causal Atlas — figures
# =============================================================================
#
# Reads the Causal Atlas parquet and produces:
#
#   Figure 1A — Volcano plot per panel: effect size (β_IVW) vs -log10(p_IVW),
#               coloured by phecode category. Lines: Bonferroni @ p<3e-5,
#               FDR-q<0.05.
#   Figure 1B — Cross-panel comparison: Tier 1 wide-excluded (Pivot A primary)
#               vs Tier 2 narrow-excluded. Concordant non-APOE hits at the top
#               right corner.
#   Figure 1C — Top-hits table (text/PDF): top 30 rows of the atlas filtered to
#               robust_hit == TRUE, ranked by IVW p-value, with within-panel
#               FDR-q.
#
# Each figure is written as both PNG (300 DPI for inspection) and PDF
# (vector for publication). All saved to results/figures/.
#
# Run as:
#   Rscript R/10_atlas_figures.R [atlas_path]
# Where atlas_path defaults to data/processed/amyloid_pet_causal_atlas.parquet
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_figs  <- file.path(proj_root, "results", "figures")
dir_tabs  <- file.path(proj_root, "results", "tables")
dir.create(dir_figs, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tabs, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
atlas_fp <- if (length(args) >= 1) args[1] else
            file.path(proj_root, "data/processed/amyloid_pet_causal_atlas.parquet")
if (!file.exists(atlas_fp)) {
  stop(sprintf("Atlas parquet not found: %s\nRun R/09_build_causal_atlas.R first.",
               atlas_fp))
}

atlas <- as.data.table(read_parquet(atlas_fp))
cat(sprintf("[%s] Atlas: %d rows × %d cols\n", Sys.time(), nrow(atlas), ncol(atlas)))

# --- Panel labels + colors ---------------------------------------------------

PANEL_LABELS <- c(
  "tier1_apoe_inclusive"    = "Tier 1 APOE-inclusive\n(descriptive only)",
  "tier1_apoe_narrow_excl"  = "Tier 1 APOE-narrow excluded\n(registered primary, n=5)",
  "tier1_apoe_wide_excl"    = "Tier 1 APOE-wide excluded\n(Pivot A primary, n=3)",
  "tier2_apoe_narrow_excl"  = "Tier 2 narrow excluded\n(MR-RAPS primary, n=13)",
  "tier2_apoe_wide_excl"    = "Tier 2 wide excluded\n(MR-RAPS sens, n=11)"
)

# Bonferroni threshold across phecodes within panel (rough)
n_phecodes <- uniqueN(atlas$phecode[!is.na(atlas$phecode)])
bonferroni_p <- 0.05 / n_phecodes

# --- Figure 1A — Volcano per panel -------------------------------------------

mk_volcano <- function(panel_data, panel_name) {
  d <- copy(panel_data[!is.na(pval_ivw) & !is.na(b_ivw)])
  if (nrow(d) == 0) return(NULL)
  d[, neglog10p := -log10(pval_ivw)]
  d[, significant := pval_ivw < bonferroni_p]
  d[, label_for_plot := ifelse(significant, as.character(phenotype), NA_character_)]

  ggplot(d, aes(b_ivw, neglog10p)) +
    geom_point(aes(colour = category, size = significant), alpha = 0.7) +
    geom_hline(yintercept = -log10(bonferroni_p), linetype = "dashed",
               colour = "grey40") +
    geom_text_repel(aes(label = label_for_plot), size = 2.5, max.overlaps = 20,
                    box.padding = 0.4) +
    scale_size_manual(values = c(`TRUE` = 2.5, `FALSE` = 1.2), guide = "none") +
    labs(
      title = paste("Panel:", PANEL_LABELS[panel_name] %||% panel_name),
      x = expression("Causal effect estimate "*beta[IVW]),
      y = expression(-log[10](p[IVW])),
      colour = "ICD-10 chapter"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          plot.title = element_text(size = 9, face = "bold"))
}

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

panels <- names(PANEL_LABELS)
volcanoes <- list()
for (pn in panels) {
  p <- mk_volcano(atlas[panel == pn], pn)
  if (!is.null(p)) volcanoes[[pn]] <- p
}

if (length(volcanoes) > 0) {
  combined <- wrap_plots(volcanoes, ncol = 2, guides = "collect") +
    plot_annotation(
      title = "Amyloid-PET Causal Atlas — locus-level MR (IVW) per instrument panel",
      subtitle = sprintf(
        "%d phecodes × %d panels; dashed line: Bonferroni p<%.1e",
        n_phecodes, length(panels), bonferroni_p),
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
  ggsave(file.path(dir_figs, "fig1A_volcano_per_panel.png"),
         combined, width = 14, height = 18, dpi = 300)
  ggsave(file.path(dir_figs, "fig1A_volcano_per_panel.pdf"),
         combined, width = 14, height = 18)
  cat(sprintf("[%s] Wrote Figure 1A: volcano per panel (%d panels)\n",
              Sys.time(), length(volcanoes)))
}

# --- Figure 1B — Tier 1 wide-excl vs Tier 2 narrow-excl concordance --------

t1_wide <- atlas[panel == "tier1_apoe_wide_excl",
                 .(phecode, t1_b = b_ivw, t1_pval = pval_ivw,
                   phenotype, category)]
t2_narrow <- atlas[panel == "tier2_apoe_narrow_excl",
                   .(phecode, t2_b = b_ivw, t2_pval = pval_ivw)]
concord <- merge(t1_wide, t2_narrow, by = "phecode")
concord <- concord[!is.na(t1_b) & !is.na(t2_b)]

if (nrow(concord) > 0) {
  concord[, concordant_strong := !is.na(t1_pval) & t1_pval < 0.01 &
                                   !is.na(t2_pval) & t2_pval < 0.01 &
                                   sign(t1_b) == sign(t2_b)]
  fig1B <- ggplot(concord, aes(t1_b, t2_b)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    geom_point(aes(colour = category, size = -log10(pmin(t1_pval, t2_pval))),
               alpha = 0.7) +
    geom_text_repel(
      data = concord[concordant_strong == TRUE],
      aes(label = phenotype), size = 2.5, max.overlaps = 25,
      box.padding = 0.4
    ) +
    scale_size_continuous(range = c(0.5, 5), name = expression(-log[10](min(p)))) +
    labs(
      title = "Cross-panel concordance: Tier 1 wide-excluded vs Tier 2 narrow-excluded",
      subtitle = "Non-APOE hits at top-right (positive in both) and bottom-left (negative in both) are robust; off-diagonal points are panel-specific.",
      x = expression(beta[IVW]*" — Tier 1 wide-excluded (3 IVs; Pivot A primary)"),
      y = expression(beta[IVW]*" — Tier 2 narrow-excluded (13 IVs; MR-RAPS primary)"),
      colour = "ICD-10 chapter"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  ggsave(file.path(dir_figs, "fig1B_concordance_t1wide_t2narrow.png"),
         fig1B, width = 11, height = 9, dpi = 300)
  ggsave(file.path(dir_figs, "fig1B_concordance_t1wide_t2narrow.pdf"),
         fig1B, width = 11, height = 9)
  cat(sprintf("[%s] Wrote Figure 1B: T1-wide vs T2-narrow concordance (%d phecodes)\n",
              Sys.time(), nrow(concord)))
}

# --- Top hits table ----------------------------------------------------------

top_n <- 30
top_robust <- atlas[robust_hit == TRUE][order(pval_ivw)][1:min(top_n, .N)]
top_overall <- atlas[!is.na(pval_ivw)][order(pval_ivw)][1:min(top_n, .N)]

# Format for human reading
fmt_top <- function(d, label) {
  d[, .(phecode, panel, phenotype, category, n_iv,
        beta_IVW = round(b_ivw, 4),
        SE_IVW   = round(se_ivw, 4),
        p_IVW    = signif(pval_ivw, 3),
        FDR_q    = signif(ivw_fdr_q_within_panel, 3),
        Q_p      = signif(q_pval, 3),
        Egger_int_p = signif(egger_intercept_pval, 3),
        coloc_PP_H4 = if ("best_pp_h4" %in% names(d)) round(best_pp_h4, 3) else NA,
        robust_hit, coloc_strengthened
  )][1:min(top_n, .N)]
}

fwrite(fmt_top(top_robust, "robust"),
       file.path(dir_tabs, "top_robust_hits.tsv"), sep = "\t")
fwrite(fmt_top(top_overall, "overall"),
       file.path(dir_tabs, "top_overall_hits.tsv"), sep = "\t")
cat(sprintf("[%s] Wrote top-hits tables to results/tables/\n", Sys.time()))

cat("\nTop 20 robust hits (most likely manuscript Table 2):\n")
print(fmt_top(top_robust, "robust")[1:min(20, .N)])
