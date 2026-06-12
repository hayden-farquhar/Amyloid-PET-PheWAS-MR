# =============================================================================
# Script 21: Quantifying the scan — power floor + empirical-null calibration
# =============================================================================
#
# Reframe add-on A. Converts the "phenome-wide null is uninformative" limitation
# into two quantified results for the conservative wide-APOE-excluded panel
# (n_IV = 3; the Pivot A per-locus primary):
#
#   (i)  Per-phecode minimum detectable odds ratio at 80% power
#        (Burgess 2014, Int J Epidemiol — binary-outcome MR power). Puts a
#        hard, citable floor under the "no non-AD signal" claim.
#   (ii) Empirical-null calibration of the 1,574 wide-panel IVW p-values
#        (genomic-inflation lambda + QQ), with and without the CR1/AD-spectrum
#        hits, testing whether the scan is well-calibrated or inflated.
#
# Power model (Burgess 2014): for a binary outcome the expected chi-square_1 of
# the IVW estimate is  ncp = N * R2 * beta1^2 * K * (1-K), where N = total
# sample, K = case fraction, R2 = exposure variance explained by the panel,
# beta1 = log-OR per SD of exposure. Minimum detectable beta1 at power (1-b),
# two-sided alpha:  beta1_min = (z_{1-a/2} + z_{1-b}) / sqrt(N*R2*K*(1-K)).
#
# Input:  data/processed/amyloid_pet_causal_atlas.parquet
#         data/interim/instruments/tier1_apoe_wide_excl_index_snps.tsv
# Output: data/processed/power_floor_wide_panel.{parquet,tsv}
#         results/tables/calibration_summary.tsv
#         results/figures/fig_power_floor/{power_floor,qq_calibration}.{pdf,png}
# =============================================================================

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(ggplot2); library(patchwork)
})

proj   <- normalizePath(file.path(dirname(sub("--file=", "",
            grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
procdir <- file.path(proj, "data", "processed")
figdir  <- file.path(proj, "results", "figures", "fig_power_floor")
tabdir  <- file.path(proj, "results", "tables")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)

ok <- list(blue = "#0072B2", vermil = "#D55E00", grey = "#999999",
           green = "#009E73", orange = "#E69F00")
base_theme <- theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e6e6", linewidth = 0.3),
        plot.caption = element_text(size = 6, colour = ok$grey),
        legend.title = element_blank(), legend.key.size = unit(0.4, "cm"))

# ---- instrument R2 for the wide-APOE-excluded panel (3 IVs) ------------------
ivs <- fread(file.path(proj, "data", "interim", "instruments",
                       "tier1_apoe_wide_excl_index_snps.tsv"))
# per-SNP R2 = 2 p (1-p) beta^2  (exposure standardised to SD units; the IVW
# beta on FinnGen is log-OR per 1 SD genetically-predicted amyloid-PET)
ivs[, r2_snp := 2 * eaf * (1 - eaf) * beta^2]
R2 <- sum(ivs$r2_snp)
cat(sprintf("Wide-APOE-excluded panel R2 = %.5f (%.3f%%) across %d IVs\n",
            R2, 100 * R2, nrow(ivs)))

# ---- per-phecode power floor ------------------------------------------------
atlas <- as.data.table(read_parquet(file.path(procdir,
                                              "amyloid_pet_causal_atlas.parquet")))
w <- atlas[panel == "tier1_apoe_wide_excl"]
w[, N := num_cases + num_controls]
w[, K := num_cases / N]

mdor <- function(N, K, R2, power = 0.80, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2) + qnorm(power)
  b1 <- z / sqrt(N * R2 * K * (1 - K))
  exp(b1)
}
power_at <- function(OR, N, K, R2, alpha = 0.05) {
  ncp <- N * R2 * (log(OR))^2 * K * (1 - K)
  pchisq(qchisq(1 - alpha, 1), df = 1, ncp = ncp, lower.tail = FALSE)
}

alpha_bonf <- 0.05 / nrow(w)            # within-panel Bonferroni
w[, `:=`(
  mdor80_nominal = mdor(N, K, R2, 0.80, 0.05),
  mdor80_bonf    = mdor(N, K, R2, 0.80, alpha_bonf),
  power_or1p5_nominal = power_at(1.5, N, K, R2, 0.05),
  power_or2_nominal   = power_at(2.0, N, K, R2, 0.05),
  power_or2_bonf      = power_at(2.0, N, K, R2, alpha_bonf)
)]

out <- w[, .(phecode, phenotype, category, num_cases, num_controls, N, K,
             b_ivw, se_ivw, pval_ivw, robust_hit,
             mdor80_nominal, mdor80_bonf,
             power_or1p5_nominal, power_or2_nominal, power_or2_bonf)]
setorder(out, mdor80_bonf)
write_parquet(out, file.path(procdir, "power_floor_wide_panel.parquet"))
fwrite(out, file.path(procdir, "power_floor_wide_panel.tsv"), sep = "\t")

# ---- headline power statements ----------------------------------------------
n_tot <- nrow(out)
frac <- function(x) sprintf("%d/%d (%.1f%%)", sum(x), n_tot, 100 * mean(x))
cat("\n--- POWER FLOOR (wide-APOE-excluded panel, R2 =", round(R2, 4), ") ---\n")
cat("Median min detectable OR @80% power, nominal a=0.05 :",
    round(median(out$mdor80_nominal), 2), "\n")
cat("Median min detectable OR @80% power, Bonferroni     :",
    round(median(out$mdor80_bonf), 2), "\n")
cat("Phecodes powered (>=80%) to detect OR>=1.5 @nominal  :",
    frac(out$power_or1p5_nominal >= 0.80), "\n")
cat("Phecodes powered (>=80%) to detect OR>=2.0 @nominal  :",
    frac(out$power_or2_nominal   >= 0.80), "\n")
cat("Phecodes powered (>=80%) to detect OR>=2.0 @Bonf     :",
    frac(out$power_or2_bonf      >= 0.80), "\n")

# ---- empirical-null calibration ---------------------------------------------
# genomic inflation lambda from the 1,574 IVW p-values
pv <- w[!is.na(pval_ivw) & pval_ivw > 0, pval_ivw]
chisq_obs <- qchisq(1 - pv, df = 1)
lambda_all <- median(chisq_obs) / qchisq(0.5, 1)

# AD-spectrum / CR1-driven phecodes to exclude for the "non-AD bulk" lambda
ad_pat <- "DEMENT|ALZ|AD_|NEURODEG|F5_VASCDEM"
w[, is_ad := grepl(ad_pat, phecode, ignore.case = TRUE)]
pv_noad <- w[is_ad == FALSE & !is.na(pval_ivw) & pval_ivw > 0, pval_ivw]
lambda_noad <- median(qchisq(1 - pv_noad, 1)) / qchisq(0.5, 1)

cal <- data.table(
  metric = c("n_phecodes_tested", "lambda_GC_all",
             "n_nonAD", "lambda_GC_nonAD", "panel_R2"),
  value  = c(length(pv), round(lambda_all, 4),
             length(pv_noad), round(lambda_noad, 4), round(R2, 5)))
fwrite(cal, file.path(tabdir, "calibration_summary.tsv"), sep = "\t")
cat("\n--- CALIBRATION ---\n")
cat("lambda_GC (all 1,574)      :", round(lambda_all, 3), "\n")
cat("lambda_GC (non-AD bulk)    :", round(lambda_noad, 3),
    "  [n =", length(pv_noad), "]\n")

# ---- figures ----------------------------------------------------------------
# (1) power-floor: min detectable OR vs case count, AD hits highlighted
pf <- out[is.finite(mdor80_bonf)]
p_floor <- ggplot(pf, aes(num_cases, mdor80_bonf)) +
  geom_point(colour = ok$grey, alpha = 0.35, size = 0.7) +
  geom_point(data = pf[robust_hit == TRUE],
             colour = ok$vermil, size = 1.4) +
  geom_hline(yintercept = 2, linetype = "dashed",
             colour = ok$blue, linewidth = 0.4) +
  scale_x_log10(labels = scales::comma) + scale_y_log10() +
  labs(x = "FinnGen R12 cases (log scale)",
       y = "Min detectable OR @80% power\n(within-panel Bonferroni)",
       title = "A  Detection floor of the 3-IV wide-APOE-excluded scan",
       caption = "Dashed line OR = 2; red = robust hits (CR1 locus)") +
  base_theme

# (2) QQ plot of observed vs expected -log10(p)
qq <- data.table(p = sort(pv))
qq[, `:=`(obs = -log10(p), exp = -log10(ppoints(length(p))))]
qq_noad <- data.table(p = sort(pv_noad))
qq_noad[, `:=`(obs = -log10(p), exp = -log10(ppoints(length(p))))]
p_qq <- ggplot(qq, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0, colour = ok$grey,
              linetype = "dashed", linewidth = 0.4) +
  geom_point(colour = ok$blue, size = 0.7, alpha = 0.5) +
  geom_point(data = qq_noad, colour = ok$green, size = 0.6, alpha = 0.5) +
  labs(x = expression(Expected~-log[10](p)),
       y = expression(Observed~-log[10](p)),
       title = "B  Empirical-null calibration (wide panel)",
       caption = sprintf("blue = all (lambda=%.2f); green = non-AD (lambda=%.2f)",
                         lambda_all, lambda_noad)) +
  base_theme

combo <- p_floor + p_qq + plot_layout(widths = c(1, 1))
ggsave(file.path(figdir, "power_floor_calibration.pdf"), combo,
       width = 9.5, height = 3.6)
ggsave(file.path(figdir, "power_floor_calibration.png"), combo,
       width = 9.5, height = 3.6, dpi = 300)

cat("\nWrote:\n  ", file.path(procdir, "power_floor_wide_panel.parquet"),
    "\n  ", file.path(tabdir, "calibration_summary.tsv"),
    "\n  ", file.path(figdir, "power_floor_calibration.pdf"), "\n")
