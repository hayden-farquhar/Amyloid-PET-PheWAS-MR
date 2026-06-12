# =============================================================================
# Script 24: Full-pairs coloc regeneration (corrects the Phase-4 significant-
#            pairs artefact across the whole coloc layer)
# =============================================================================
#
# Phase 4 ran coloc.abf against GTEx SIGNIFICANT-pairs eQTL data (only FDR-
# significant SNP-gene pairs in the window). At tight-LD loci this truncation
# strips the variants that distinguish nearby signals and inflates PP.H4 toward
# false-positive colocalisation (demonstrated for the forearm-fracture negative
# control: 37 sig-pairs SNPs -> PP.H4 0.983; 2,465 full-pairs SNPs -> 0.107).
#
# This script re-runs every (phecode x tissue x gene) combination recorded in
# coloc_results.parquet using FULL-pairs nominal eQTL data (byte-range tabix
# from the EBI eQTL Catalogue GTEx_V8 mirror) and full FinnGen R12 regional
# sumstats, then rebuilds the atlas coloc columns.
#
# Scope: 16 phecodes x 2 loci (CR1 chr1:207,518,704; APOE chr19:44,888,997)
#        x 12 brain tissues x 18 genes = 527 coloc runs.
#
# coloc.abf is direction- and LD-agnostic (per-SNP lABF from beta/varbeta under
# the single-causal-variant model), so SNPs are matched on position+alleles
# only; no allele-sign alignment is required.
#
# Output: data/processed/coloc_results_fullpairs.parquet
#         data/processed/amyloid_pet_causal_atlas.parquet (coloc cols updated;
#           prior version preserved with timestamp suffix)
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(coloc); library(arrow) })

proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
procdir <- file.path(proj, "data", "processed")
cache   <- file.path(proj, "data", "interim", "coloc_fullpairs")
fg_cache <- file.path(proj, "data", "interim", "coloc_symmetry")  # reuse FinnGen pulls
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
dir.create(fg_cache, showWarnings = FALSE, recursive = TRUE)

EQTL_BASE <- "http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge"
WIN <- 500000L
loci <- data.table(
  locus = c("CR1", "APOE"), chr = c(1L, 19L),
  bp = c(207518704L, 44888997L))
loci[, `:=`(lo = bp - WIN, hi = bp + WIN)]
loci[, region := sprintf("%d:%d-%d", chr, lo, hi)]

prev <- as.data.table(read_parquet(file.path(procdir, "coloc_results.parquet")))
prev[, gene_base := sub("\\..*$", "", gene)]
combos <- unique(prev[, .(phecode, tissue, gene, gene_base, lead_chr, lead_bp)])
phecodes <- unique(prev$phecode)
tissues  <- unique(prev$tissue)
cat(sprintf("Regenerating %d coloc combos: %d phecodes x %d tissues x %d genes\n",
            nrow(combos), length(phecodes), length(tissues),
            uniqueN(prev$gene)))

meta <- fread(file.path(proj, "osf", "phecode_filter_prespec.csv"))
setnames(meta, c("phenocode","phenotype","category","chapter_prefix",
                 "num_cases","num_controls","url"))

# ---- fetch helpers (cache .tbi in dedicated cwd to avoid root spam) ----------
oldwd <- getwd(); setwd(cache)
on.exit({ unlink(file.path(cache, "*.tbi")); setwd(oldwd) }, add = TRUE)

fetch_eqtl <- function(tissue) {
  f <- file.path(cache, sprintf("eqtl_%s.tsv", tissue))
  if (file.exists(f) && file.info(f)$size > 1000) return(f)
  url <- sprintf("%s/%s.tsv.gz", EQTL_BASE, tissue)
  # pull both locus windows in one tabix call
  args <- c(shQuote(url), loci$region)
  system2("tabix", args, stdout = f)
  f
}
fetch_finngen <- function(ph) {
  f <- file.path(fg_cache, sprintf("%s.tsv", ph))
  if (file.exists(f) && file.info(f)$size > 100) return(f)
  url <- meta[phenocode == ph, url]
  lc <- combos[combos$phecode == ph, ]$lead_chr[1]
  reg <- loci[chr == lc, region]
  system2("tabix", c(shQuote(url), reg), stdout = f)
  f
}

eqtl_cols <- c("variant_id","r2","pvalue","mtid","mtoid","maf","gene_id",
               "median_tpm","beta","se","an","ac","chr_e","pos_e","ref_e",
               "alt_e","type","rsid")
fg_cols <- c("chrom","pos","ref","alt","rsids","genes","pval","mlogp",
             "beta","sebeta","af_alt","af_cc","af_ctrl")

cat("Fetching full-pairs eQTL per tissue (both loci windows)...\n")
eqtl_files <- setNames(sapply(tissues, fetch_eqtl), tissues)
cat("Fetching FinnGen regional sumstats per phecode...\n")
fg_files <- setNames(sapply(phecodes, fetch_finngen), phecodes)

# ---- load + index eQTL by (tissue) ------------------------------------------
load_eqtl <- function(f) {
  d <- fread(f, header = FALSE, col.names = eqtl_cols)
  d <- d[se > 0 & is.finite(beta) & maf > 0.01]
  d[, key := paste(pos_e, pmin(ref_e, alt_e), pmax(ref_e, alt_e), sep = "_")]
  d
}
eqtl_dt <- lapply(eqtl_files, load_eqtl)

load_fg <- function(f) {
  d <- fread(f, header = FALSE, col.names = fg_cols)
  d <- d[sebeta > 0 & is.finite(beta) & af_alt > 0.01 & af_alt < 0.99]
  d[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
  d[!duplicated(key)]
}
fg_dt <- lapply(fg_files, load_fg)

# ---- coloc per (phecode, tissue, gene) --------------------------------------
res <- vector("list", nrow(combos))
for (i in seq_len(nrow(combos))) {
  cb <- combos[i]
  ph <- fg_dt[[cb$phecode]]
  eqt <- eqtl_dt[[cb$tissue]]
  if (is.null(ph) || is.null(eqt)) next
  g <- eqt[gene_id == cb$gene_base]
  g <- g[!duplicated(key)]
  common <- intersect(ph$key, g$key)
  if (length(common) < 20) { next }
  p2 <- ph[key %in% common][order(key)]
  g2 <- g[key %in% common][order(key)]
  mrow <- meta[phenocode == cb$phecode]
  N <- mrow$num_cases + mrow$num_controls; s <- mrow$num_cases / N
  N_eq <- max(round(median(g2$an, na.rm = TRUE) / 2), 50)
  D1 <- list(beta = p2$beta, varbeta = p2$sebeta^2, snp = p2$key,
             type = "cc", N = N, s = s, MAF = pmin(p2$af_alt, 1 - p2$af_alt))
  D2 <- list(beta = g2$beta, varbeta = g2$se^2, snp = g2$key, type = "quant",
             N = N_eq, MAF = g2$maf)
  cc <- tryCatch(suppressMessages(suppressWarnings(coloc.abf(D1, D2))),
                 error = function(e) NULL)
  if (is.null(cc)) next
  res[[i]] <- data.table(
    gene = cb$gene, tissue = cb$tissue, n_shared_snps = length(common),
    PP.H0 = cc$summary["PP.H0.abf"], PP.H1 = cc$summary["PP.H1.abf"],
    PP.H2 = cc$summary["PP.H2.abf"], PP.H3 = cc$summary["PP.H3.abf"],
    PP.H4 = cc$summary["PP.H4.abf"], phecode = cb$phecode,
    lead_chr = cb$lead_chr, lead_bp = cb$lead_bp)
}
fullres <- rbindlist(res)
setnames(fullres, c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"),
         c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4"))  # keep names
write_parquet(fullres, file.path(procdir, "coloc_results_fullpairs.parquet"))
fwrite(fullres, file.path(procdir, "coloc_results_fullpairs.tsv"), sep = "\t")

# ---- compare best PP.H4 per phecode: sig-pairs vs full-pairs -----------------
best_old <- prev[, .(pph4_sigpairs = max(PP.H4)), by = phecode]
best_new <- fullres[, .(pph4_fullpairs = max(PP.H4),
                        top_gene = gene[which.max(PP.H4)],
                        top_tissue = tissue[which.max(PP.H4)],
                        n_gt0p5 = sum(PP.H4 > 0.5),
                        n_gt0p8 = sum(PP.H4 > 0.8)), by = phecode]
cmp <- merge(best_old, best_new, by = "phecode", all = TRUE)
setorder(cmp, -pph4_sigpairs)
cat("\n--- best PP.H4 per phecode: significant-pairs (Phase 4) vs full-pairs ---\n")
print(cmp[, .(phecode, pph4_sigpairs = round(pph4_sigpairs, 3),
              pph4_fullpairs = round(pph4_fullpairs, 3),
              top_gene = sub("\\..*$","",top_gene), top_tissue)])
fwrite(cmp, file.path(procdir, "coloc_sigpairs_vs_fullpairs.tsv"), sep = "\t")

# ---- update atlas coloc columns ---------------------------------------------
ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
atlas <- as.data.table(read_parquet(file.path(procdir,
                                              "amyloid_pet_causal_atlas.parquet")))
write_parquet(atlas, file.path(procdir,
              sprintf("amyloid_pet_causal_atlas_PRE_FULLPAIRS_%s.parquet", ts)))

bn <- best_new[, .(phecode, fp_best = pph4_fullpairs, fp_gene = top_gene,
                   fp_tissue = top_tissue, fp_n05 = n_gt0p5, fp_n08 = n_gt0p8)]
atlas <- merge(atlas, bn, by = "phecode", all.x = TRUE)
# only phecodes that were coloc-tested get updated; others keep NA->0
atlas[!is.na(fp_best), `:=`(
  best_pp_h4 = fp_best, top_coloc_gene = fp_gene, top_coloc_tissue = fp_tissue,
  n_coloc_pph4_gt0p5 = fp_n05, n_coloc_pph4_gt0p8 = fp_n08)]
atlas[, coloc_strengthened := (robust_hit == TRUE) & !is.na(best_pp_h4) &
        best_pp_h4 > 0.5]
atlas[, c("fp_best","fp_gene","fp_tissue","fp_n05","fp_n08") := NULL]

write_parquet(atlas, file.path(procdir, "amyloid_pet_causal_atlas.parquet"))
fwrite(atlas, file.path(proj, "results", "amyloid_pet_causal_atlas.tsv"), sep = "\t")

n_cs_after <- atlas[coloc_strengthened == TRUE, .N]
cat(sprintf("\nColoc-strengthened hits: was 20 (sig-pairs) -> %d (full-pairs)\n",
            n_cs_after))
cat("Atlas updated; prior version preserved with suffix _PRE_FULLPAIRS_", ts, "\n")
