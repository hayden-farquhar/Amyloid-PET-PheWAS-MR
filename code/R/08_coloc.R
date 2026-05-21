# =============================================================================
# Script 08: Phase 4 — Colocalisation (coloc.abf) at lead loci
# =============================================================================
#
# For each Phase 2 phecode hit (default: panel-IVW p < 1e-4 in any Tier 2 or
# Tier 1 wide-excl panel), identify the lead instrument SNP, define a ±500 kb
# window, and test whether the amyloid-PET signal colocalises with eQTL
# signal at any gene in that window. Tissues queried:
#
#   - GTEx v8 brain (13 tissues; EUR-only signif_pairs)
#   - eQTLGen Phase 1 whole blood (significant cis-eQTLs FDR<0.05)
#
# Per pre-reg §17 colocalisation framework. Primary method: coloc.abf
# (Wakefield 2009 / Giambartolomei 2014). Multi-causal SuSiE extension (§17.7)
# is deferred to a second script once we have hits worth deep-investigating.
#
# Input:
#   - data/processed/locus_mr_results.parquet (or checkpoint while Phase 2 running)
#   - data/raw/ali_2023/AmyloidPET_GWAS_multiEthnic_metaAnalysis_*.txt
#   - data/reference/gtex_v8/brain/eqtls/*.signif_pairs.txt.gz
#   - data/reference/eqtlgen/cis-eQTLsFDR0.05.txt.gz
#   - data/interim/instruments/unified_iv_list.tsv
#
# Output:
#   - data/processed/coloc_results.parquet
#
# DATA CAVEAT: GTEx ships signif_pairs only (FDR<0.05 per tissue). Full-pairs
# data are 50-200 GB per tissue and not downloaded. This means coloc.abf
# priors will be approximate — the manuscript Methods explicitly reports this
# limitation and treats high-PP.H4 hits with significant-pairs data as a
# *screen*, not a definitive coloc claim. coloc.susie with full-pairs data
# would be the rigorous follow-up for confirmed hits.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(R.utils)
  library(coloc)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)
dir_interim   <- file.path(proj_root, "data", "interim")
dir_processed <- file.path(proj_root, "data", "processed")
dir_ref       <- file.path(proj_root, "data", "reference")

ALI_SUMSTATS <- file.path(proj_root, "data/raw/ali_2023",
  "AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt")
GTEX_DIR <- file.path(dir_ref, "gtex_v8/brain/eqtls")
EQTLGEN_FP <- file.path(dir_ref, "eqtlgen/cis-eQTLsFDR0.05.txt.gz")

SEED <- 81
set.seed(SEED)

# --- Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
p_threshold <- if (length(args) >= 1) as.numeric(args[1]) else 1e-4
window_kb   <- if (length(args) >= 2) as.integer(args[2]) else 500L
N_amyloid   <- 13409  # Ali 2023 multi-ethnic
N_gtex_brain <- list(  # approximate per-tissue sample sizes (GTEx v8 EUR)
  Brain_Cortex = 205, Brain_Hippocampus = 165,
  Brain_Anterior_cingulate_cortex_BA24 = 147, Brain_Amygdala = 129,
  Brain_Hypothalamus = 170, Brain_Substantia_nigra = 114,
  Brain_Caudate_basal_ganglia = 194, Brain_Putamen_basal_ganglia = 170,
  Brain_Cerebellum = 209, Brain_Cerebellar_Hemisphere = 175,
  Brain_Frontal_Cortex_BA9 = 175,
  Brain_Nucleus_accumbens_basal_ganglia = 202,
  Brain_Spinal_cord_cervical_c.1 = 126
)
N_eqtlgen <- 31684

# --- Identify lead loci from Phase 2 -----------------------------------------

# Prefer full results; fall back to checkpoint
phase2_fp <- file.path(dir_processed, "locus_mr_results.parquet")
phase2_ck <- file.path(dir_processed, "locus_mr_checkpoint.parquet")
ck_size <- if (file.exists(phase2_ck)) file.size(phase2_ck) else 0
full_size <- if (file.exists(phase2_fp)) file.size(phase2_fp) else 0
src <- if (full_size > ck_size) phase2_fp else phase2_ck
cat(sprintf("[%s] Reading Phase 2 from %s\n", Sys.time(), basename(src)))
phase2 <- as.data.table(read_parquet(src))

# Identify candidate phecodes with IVW p < threshold in panels of interest
target_panels <- c("tier1_apoe_wide_excl",
                    "tier2_apoe_narrow_excl", "tier2_apoe_wide_excl")
candidates <- phase2[!is.na(pval) & pval < p_threshold &
                       panel %in% target_panels & method == "ivw",
                     .(phecode, panel, n_iv, b, se, pval, phenotype, category)]
cat(sprintf("[%s] Candidate phecodes (IVW p < %.0e in %d panels): %d unique\n",
            Sys.time(), p_threshold, length(target_panels),
            uniqueN(candidates$phecode)))

if (nrow(candidates) == 0) {
  cat("No candidate hits found at this threshold; exiting.\n")
  cat("(Re-run with --p_threshold relaxed if Phase 2 finished and screen is empty.)\n")
  quit(status = 0)
}

# Pick lead IVs by intersecting top candidates with the unified IV list
iv_unified <- fread(file.path(dir_interim, "instruments/unified_iv_list.tsv"))

# Define lead loci: take the strongest panel-method estimate per phecode,
# then map to the IV(s) in that panel
lead_loci <- list()
for (pid in unique(candidates$phecode)) {
  top_row <- candidates[phecode == pid][order(pval)][1]
  panel <- top_row$panel
  ivs_for_panel <- iv_unified[grepl(panel, panels, fixed = TRUE)]
  # Use the strongest IV in the panel (highest F-stat) as the locus anchor
  if (nrow(ivs_for_panel) == 0) next
  anchor <- ivs_for_panel[order(-f_stat)][1]
  lead_loci[[length(lead_loci) + 1]] <- data.table(
    phecode = pid, panel = panel,
    phecode_pval = top_row$pval,
    lead_rsid = anchor$rsid, lead_chr = as.integer(anchor$chr),
    lead_bp_hg38 = as.integer(anchor$hg38_pos),
    win_lo = as.integer(anchor$hg38_pos) - window_kb * 1000L,
    win_hi = as.integer(anchor$hg38_pos) + window_kb * 1000L
  )
}
lead_loci <- rbindlist(lead_loci)
cat(sprintf("[%s] Lead loci to coloc: %d\n", Sys.time(), nrow(lead_loci)))
print(lead_loci[, .(phecode, panel, lead_rsid, lead_chr, lead_bp_hg38)])

# --- Read Ali 2023 sumstats lazily per locus ---------------------------------

read_ali_locus <- function(chr, lo, hi) {
  # Ali 2023 columns: CHR, BP, EFFECT_ALLELE, OTHER_ALLELE, BETA, SE, P, ...
  # Stream the file via awk to filter to the locus window
  cmd <- sprintf(
    "awk -v c=%d -v l=%d -v h=%d 'NR==1 || ($1==c && $2>=l && $2<=h)' %s",
    chr, lo, hi, shQuote(ALI_SUMSTATS)
  )
  fread(cmd = cmd, na.strings = c("NA", ""))
}

# --- Read GTEx tissue eQTL within locus --------------------------------------

read_gtex_locus <- function(tissue, chr, lo, hi) {
  fp <- file.path(GTEX_DIR, sprintf("Brain_%s.v8.EUR.signif_pairs.txt.gz", tissue))
  if (!file.exists(fp)) return(NULL)
  # variant_id format: chr1_14677_G_A_b38
  d <- fread(fp)
  d[, c("var_chr", "var_pos", "var_ref", "var_alt", "build") :=
       tstrsplit(variant_id, "_", fixed = TRUE)]
  d[, var_chr := as.integer(sub("chr", "", var_chr))]
  d[, var_pos := as.integer(var_pos)]
  d[var_chr == chr & var_pos >= lo & var_pos <= hi]
}

# --- Read eQTLGen within locus -----------------------------------------------

# Build hg37-coords from rsIDs if needed; eQTLGen uses GRCh37 positions
# while Ali 2023 uses GRCh38. Match by rsID when possible.
read_eqtlgen_locus <- function(rsid_list) {
  if (length(rsid_list) == 0) return(NULL)
  # Use awk to filter the .gz to matching rsIDs (very fast)
  rsid_grep <- paste(rsid_list, collapse = "|")
  cmd <- sprintf("gunzip -c %s | awk -F'\\t' 'NR==1 || $2 ~ /%s/'",
                 shQuote(EQTLGEN_FP), rsid_grep)
  d <- tryCatch(fread(cmd = cmd, na.strings = c("NA", "")),
                error = function(e) NULL)
  d
}

# --- Coloc.abf wrapper -------------------------------------------------------
#
# For each (locus × tissue × gene), build two coloc input lists:
#   ds1 (exposure: amyloid-PET): SNP-level beta, varbeta, MAF, type='quant', N
#   ds2 (outcome: eQTL): same
# Then coloc::coloc.abf returns PP.H0/H1/H2/H3/H4.
#   PP.H4 = probability of shared causal variant (the coloc signal of interest).
coloc_one <- function(ali_d, eqtl_d, gene_id, tissue, n_eqtl) {
  # Restrict to SNPs in both datasets — match by chr+pos (Ali is GRCh38)
  if (is.null(ali_d) || nrow(ali_d) == 0 ||
      is.null(eqtl_d) || nrow(eqtl_d) == 0) return(NULL)
  # eqtl_d from GTEx has var_chr, var_pos, slope, slope_se, maf
  m <- merge(ali_d[, .(BP, BETA, SE, EAF = as.numeric(ALT_FREQS),
                       A1 = EFFECT_ALLELE, A2 = OTHER_ALLELE)],
             eqtl_d[, .(var_pos, slope, slope_se, maf, var_ref, var_alt)],
             by.x = "BP", by.y = "var_pos")
  m <- m[!is.na(BETA) & !is.na(SE) & !is.na(EAF) & SE > 0]
  if (nrow(m) < 10) return(NULL)
  ds1 <- list(beta = m$BETA, varbeta = m$SE^2, MAF = pmin(m$EAF, 1 - m$EAF),
              type = "quant", N = N_amyloid, snp = paste0("snp_", m$BP))
  ds2 <- list(beta = m$slope, varbeta = m$slope_se^2, MAF = m$maf,
              type = "quant", N = n_eqtl, snp = paste0("snp_", m$BP))
  res <- tryCatch(coloc::coloc.abf(ds1, ds2),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  data.table(
    gene = gene_id, tissue = tissue, n_shared_snps = nrow(m),
    PP.H0 = res$summary["PP.H0.abf"],
    PP.H1 = res$summary["PP.H1.abf"],
    PP.H2 = res$summary["PP.H2.abf"],
    PP.H3 = res$summary["PP.H3.abf"],
    PP.H4 = res$summary["PP.H4.abf"]
  )
}

# --- Driver per lead locus ---------------------------------------------------

cat(sprintf("[%s] Starting coloc sweep over %d lead loci × ~13 brain tissues + blood\n",
            Sys.time(), nrow(lead_loci)))

all_results <- list()
brain_tissues <- names(N_gtex_brain)

for (i in seq_len(nrow(lead_loci))) {
  L <- lead_loci[i]
  cat(sprintf("\n[%s] Locus %d/%d: %s (%s) — chr%d:%d±%dkb\n",
              Sys.time(), i, nrow(lead_loci), L$phecode, L$lead_rsid,
              L$lead_chr, L$lead_bp_hg38, window_kb))

  ali_d <- tryCatch(read_ali_locus(L$lead_chr, L$win_lo, L$win_hi),
                    error = function(e) NULL)
  if (is.null(ali_d) || nrow(ali_d) < 10) {
    cat("  Skip: too few Ali SNPs in window\n"); next
  }
  cat(sprintf("  Ali SNPs in window: %d\n", nrow(ali_d)))

  # GTEx brain tissues
  for (tissue in brain_tissues) {
    tissue_short <- sub("^Brain_", "", tissue)
    gtex_d <- tryCatch(read_gtex_locus(tissue_short, L$lead_chr, L$win_lo, L$win_hi),
                       error = function(e) NULL)
    if (is.null(gtex_d) || nrow(gtex_d) == 0) next
    gene_list <- unique(gtex_d$phenotype_id)
    for (g in gene_list) {
      gtex_gene <- gtex_d[phenotype_id == g]
      cr <- coloc_one(ali_d, gtex_gene, g, tissue, N_gtex_brain[[tissue]])
      if (!is.null(cr)) {
        cr[, `:=`(phecode = L$phecode, lead_rsid = L$lead_rsid,
                  lead_chr = L$lead_chr, lead_bp = L$lead_bp_hg38)]
        all_results[[length(all_results) + 1]] <- cr
      }
    }
  }

  # eQTLGen whole blood — match by rsID (build-agnostic). For now skip,
  # since the Ali-eQTLGen position match requires lift. Place a TODO.
  # TODO: implement eQTLGen colocalisation via rsID matching once we
  # confirm Ali rsIDs are available in this locus window.
}

if (length(all_results) == 0) {
  cat("\nNo coloc results computed (no overlapping SNPs between Ali and eQTL data).\n")
  quit(status = 0)
}

results_dt <- rbindlist(all_results, fill = TRUE)
out_fp <- file.path(dir_processed,
                    sprintf("coloc_results_%s.parquet",
                             format(Sys.time(), "%Y%m%d-%H%M%S")))
write_parquet(results_dt, out_fp)
write_parquet(results_dt, file.path(dir_processed, "coloc_results.parquet"))
cat(sprintf("\n[%s] Wrote %s (%d rows)\n", Sys.time(), out_fp, nrow(results_dt)))

# Summary: top PP.H4 hits
cat("\nTop colocalisations (PP.H4 > 0.5):\n")
hits <- results_dt[PP.H4 > 0.5][order(-PP.H4)]
if (nrow(hits) > 0) {
  print(hits[1:min(20, .N),
             .(phecode, lead_rsid, gene, tissue, n_shared_snps,
               PP.H4 = round(PP.H4, 3), PP.H3 = round(PP.H3, 3))])
} else {
  cat("  (no colocalisations PP.H4 > 0.5 at this threshold)\n")
}
