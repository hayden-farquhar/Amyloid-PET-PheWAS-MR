#!/usr/bin/env Rscript
# Build a non-APOE wide-PGS IV list: top N PRS-CS SNPs outside the wide-APOE
# LD block (chr19:44.4-46.5 Mb). Output is sorted by chr, pos for tabix.

suppressPackageStartupMessages({
  library(data.table); library(arrow)
})

setwd(Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here()))

args <- commandArgs(trailingOnly = TRUE)
N_TOP <- if (length(args) >= 1) as.integer(args[1]) else 50L

files <- list.files("data/interim/prscs",
                    pattern = "_chr[0-9]+\\.txt$", full.names = TRUE)
prscs <- rbindlist(lapply(files, function(fp) {
  fread(fp, col.names = c("chr","rsid","bp","A1","A2","beta_pgs"),
        colClasses = c("integer","character","integer","character","character","numeric"))
}))
cat(sprintf("Total PRS-CS SNPs: %d\n", nrow(prscs)))

# Exclude wide-APOE LD block: chr19:44,400,000 - 46,500,000 (GRCh37 coords used by PRS-CS)
# Note: PRS-CS uses GRCh37 bim positions (from 1000G EUR.bim symlinked from P60).
# The Tier 1 wide-APOE GRCh38 exclusion is chr19:44,400,000-46,500,000.
# At chr19:44-46 Mb, GRCh37 and GRCh38 are within ~500 kb of each other; the
# wide boundary spans 2.1 Mb so the build mismatch is negligible at this scope.
prscs_noapoe <- prscs[!(chr == 19 & bp >= 44400000 & bp <= 46500000)]
cat(sprintf("After wide-APOE exclusion: %d (-%d)\n",
            nrow(prscs_noapoe), nrow(prscs) - nrow(prscs_noapoe)))

top_n <- prscs_noapoe[order(-abs(beta_pgs))][1:N_TOP]
cat(sprintf("Top %d non-APOE SNPs:\n", N_TOP))
cat(sprintf("  Weight range: %.2e to %.2e (median %.2e)\n",
            min(abs(top_n$beta_pgs)), max(abs(top_n$beta_pgs)),
            median(abs(top_n$beta_pgs))))
cat(sprintf("  Chromosomes covered: %d (%s)\n",
            uniqueN(top_n$chr), paste(sort(unique(top_n$chr)), collapse=",")))

# Sort by chr, pos for tabix
top_n <- top_n[order(chr, bp)]

# Reorder columns to match sweep_finngen_r12.sh expectation:
# rsid first (col 1), then chr (col 2), then bp/hg38_pos (col 3)
top_n <- top_n[, .(rsid, chr, hg38_pos = bp, A1, A2, beta_pgs)]

out_fp <- "data/interim/instruments/pgs_wide_iv_list_noAPOE.tsv"
fwrite(top_n, out_fp, sep = "\t")
cat(sprintf("Wrote: %s (%d SNPs)\n", out_fp, nrow(top_n)))
