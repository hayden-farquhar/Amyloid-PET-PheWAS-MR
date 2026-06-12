# =============================================================================
# Script 29: Brain_Cortex-primary colocalisation (resolves the across-tissue-max
#            anti-conservativeness; reviewer-round-3 issue 1)
# =============================================================================
# The coloc_strengthened flag and Table 2 used the MAXIMUM PP.H4 across 13 brain
# tissues, which is anti-conservative (a max over correlated, truncation-noisy
# estimates catches outliers). Brain_Cortex is the pre-specified, biologically-
# relevant tissue. This script reports Brain_Cortex as primary for all 13
# CR1-locus coloc-tested phecodes, with one INTERNALLY CONSISTENT pipeline for
# both significant-pairs and full-pairs (type='cc' phenotype, type='quant' eQTL),
# resolving the prior Phase-4-vs-R/27 sig-pairs discrepancy. Across-tissue max is
# retained as a secondary column.
#
# Output: data/processed/coloc_brain_cortex_primary.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(coloc); library(arrow) })
proj <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
if (length(proj) == 0 || is.na(proj)) proj <- getwd()
idir <- file.path(proj, "data", "interim", "coloc_susie")
fg <- file.path(proj, "data", "interim", "coloc_symmetry")
procdir <- file.path(proj, "data", "processed")
CR1 <- "ENSG00000203710"; N_EQTL <- 205

phe <- c("F5_DEMENTIA_INCLAVO","F5_DEMENTIA","G6_AD_WIDE","G6_ALZHEIMER",
         "KRA_PSY_DEMENTIA_EXMORE","F5_ALZHDEMENT","AD_LO_EXMORE",
         "G6_NEURODEGENERATIVE","NEURODEGOTH","AD_EO_EXMORE","F5_VASCDEM",
         "ST19_FRACT_FOREA","Q17_LVOTO_BROAD")
meta <- fread(file.path(proj, "osf", "phecode_filter_prespec.csv"))
setnames(meta, c("phenocode","phenotype","category","chapter_prefix",
                 "num_cases","num_controls","url"))
fg_cols <- c("chrom","pos","ref","alt","rsids","genes","pval","mlogp",
             "beta","sebeta","af_alt","af_cc","af_ctrl")

# ---- full-pairs CR1 eQTL (Brain_Cortex) -------------------------------------
ef <- fread(file.path(idir, "chr1q32.2_brain_cortex_raw.tsv"), header = FALSE,
            col.names = c("variant_id","r2","pvalue","mtid","mtoid","maf",
                          "gene_id","median_tpm","beta","se","an","ac",
                          "chr_e","pos_e","ref_e","alt_e","type","rsid"))
ef <- ef[gene_id == CR1 & se > 0 & is.finite(beta) & maf > 0.01]
ef[, key := paste(pos_e, pmin(ref_e, alt_e), pmax(ref_e, alt_e), sep = "_")]
ef <- ef[!duplicated(key)]

# ---- significant-pairs CR1 eQTL (Brain_Cortex; 40 SNPs) ---------------------
sp <- fread(cmd = paste0("gunzip -c ", shQuote(file.path(proj, "data", "reference",
  "gtex_v8", "brain", "eqtls", "Brain_Cortex.v8.EUR.signif_pairs.txt.gz")),
  " | awk -F'\\t' 'NR==1 || $1 ~ /ENSG00000203710/'"))
sp[, c("c","pos","ref","alt") := tstrsplit(variant_id, "_", keep = 1:4)]
sp[, pos := as.integer(pos)]
sp[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
sp <- sp[!duplicated(key)]

eqtl_quant <- function(d, beta, se, maf, keys) {
  x <- d[key %in% keys][order(key)]
  list(beta = x[[beta]], varbeta = x[[se]]^2, snp = x$key, type = "quant",
       N = N_EQTL, MAF = x[[maf]])
}
run <- function(d2, edt, beta, se, maf) {
  common <- intersect(d2$key, edt$key)
  if (length(common) < 10) return(NA_real_)
  p2 <- d2[key %in% common][order(key)]
  mr <- meta[phenocode == d2$ph[1]]; N <- mr$num_cases + mr$num_controls
  D1 <- list(beta = p2$beta, varbeta = p2$sebeta^2, snp = p2$key, type = "cc",
             N = N, s = mr$num_cases / N, MAF = pmin(p2$af_alt, 1 - p2$af_alt))
  cc <- tryCatch(suppressMessages(suppressWarnings(
    coloc.abf(D1, eqtl_quant(edt, beta, se, maf, common)))), error = function(e) NULL)
  if (is.null(cc)) NA_real_ else as.numeric(cc$summary["PP.H4.abf"])
}

# across-tissue max (full-pairs) from the regenerated layer
fullres <- as.data.table(read_parquet(file.path(procdir, "coloc_results_fullpairs.parquet")))
fmax <- fullres[grepl(CR1, gene), .(maxTiss_full = max(PP.H4)), by = phecode]

out <- rbindlist(lapply(phe, function(p) {
  f <- file.path(fg, paste0(p, ".tsv"))
  d <- fread(f, header = FALSE, col.names = fg_cols)
  d <- d[sebeta > 0 & is.finite(beta) & af_alt > 0.01 & af_alt < 0.99]
  d[, key := paste(pos, pmin(ref, alt), pmax(ref, alt), sep = "_")]
  d <- d[!duplicated(key)]; d[, ph := p]
  data.table(phecode = p,
             BC_sigpairs = run(d, sp, "slope", "slope_se", "maf"),
             BC_fullpairs = run(d, ef, "beta", "se", "maf"))
}))
out <- merge(out, fmax, by = "phecode", all.x = TRUE)
out[, cls := fifelse(BC_fullpairs > 0.5, "CR1-colocalising (late-onset AD-spectrum)",
                     "non-colocalising")]
setorder(out, -BC_fullpairs)
fwrite(out, file.path(procdir, "coloc_brain_cortex_primary.tsv"), sep = "\t")
cat("--- Brain_Cortex-primary coloc (consistent pipeline) ---\n")
print(out[, .(phecode, BC_sigpairs = round(BC_sigpairs,3),
              BC_fullpairs = round(BC_fullpairs,3),
              maxTiss_full = round(maxTiss_full,3))])
cat(sprintf("\nColoc-strengthened by Brain_Cortex full-pairs PP.H4 > 0.5: %d phecodes\n",
            sum(out$BC_fullpairs > 0.5, na.rm = TRUE)))
cat(sprintf("AD-spectrum (BC>0.5): full-pairs %.3f-%.3f; under sig-pairs %.3f-%.3f\n",
            min(out[BC_fullpairs>0.5]$BC_fullpairs), max(out[BC_fullpairs>0.5]$BC_fullpairs),
            min(out[BC_fullpairs>0.5]$BC_sigpairs), max(out[BC_fullpairs>0.5]$BC_sigpairs)))
cat(sprintf("Forearm fracture: BC sig-pairs %.3f -> BC full-pairs %.3f (across-tissue max %.3f)\n",
            out[phecode=='ST19_FRACT_FOREA']$BC_sigpairs,
            out[phecode=='ST19_FRACT_FOREA']$BC_fullpairs,
            out[phecode=='ST19_FRACT_FOREA']$maxTiss_full))
