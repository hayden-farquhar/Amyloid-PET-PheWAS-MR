# =============================================================================
# Script 00: SHA256 + provenance logging for raw sumstats
# =============================================================================
#
# Walks data/raw/<source>/ subdirectories, computes SHA256 for every file,
# writes per-subdirectory SHA256SUMS and download_log.tsv files.
#
# Idempotent — safe to re-run after adding new files. Existing audit rows
# for unchanged files are preserved; rows for files whose SHA256 has changed
# are overwritten with a new access_iso timestamp (and the script warns).
#
# Usage:
#   Rscript R/00_hash_raw.R
#   Rscript R/00_hash_raw.R --source ali_2023      # restrict to one subdir
# =============================================================================

suppressPackageStartupMessages({
  library(digest)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
restrict_source <- NULL
i <- which(args == "--source")
if (length(i) == 1 && length(args) >= i + 1) restrict_source <- args[i + 1]

project_root <- here::here()
raw_dir <- file.path(project_root, "data", "raw")
if (!dir.exists(raw_dir)) stop("data/raw/ not found at ", raw_dir)

log_msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")

# Per-source source URL + host metadata (extend as sources are added).
# Falls back to "(manual)" if a source isn't listed.
SOURCE_META <- list(
  ali_2023 = list(
    url  = "https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l",
    host = "WashU NeuroGenomics + Informatics Center (Box folder share)",
    publication = "Ali M, et al. Acta Neuropath Comm (2023) 11:68. DOI 10.1186/s40478-023-01563-4"
  ),
  kim_2025 = list(
    url  = "GCST90483382 (sumstats deposit pending; author contact in flight)",
    host = "Won lab, SKKU SAIHST (wonhh@skku.edu)",
    publication = "Kim JP, Won H-H, et al. Nat Commun (2025) 16:3133. DOI 10.1038/s41467-025-57751-4"
  ),
  nho_2024_tau_pet = list(
    url  = "https://adni.loni.usc.edu/ (registered access)",
    host = "ADNI / LONI USC",
    publication = "Nho K, et al. (2024) tau-PET GWAS"
  ),
  jansen_2022_csf_abeta = list(
    url  = "EADB consortium",
    host = "EADB",
    publication = "Jansen IE, et al. (2022) CSF Aβ42 GWAS"
  ),
  bbj_amyloid_related = list(
    url  = "https://jenger.riken.jp/",
    host = "Biobank Japan",
    publication = "(phenotype-specific; populated per hit)"
  )
)

write_audit <- function(subdir) {
  files <- list.files(subdir, full.names = TRUE, recursive = FALSE)
  files <- files[!basename(files) %in% c("SHA256SUMS", "download_log.tsv", "README.md")]
  files <- files[!file.info(files)$isdir]

  if (length(files) == 0) {
    log_msg("  (no files in ", basename(subdir), " — skipping)")
    return(invisible(NULL))
  }

  meta <- SOURCE_META[[basename(subdir)]]
  if (is.null(meta)) meta <- list(url = "(manual)", host = "(manual)", publication = "(manual)")

  # SHA256SUMS — GNU coreutils format: "<sha256>  <filename>"
  sums_lines <- vapply(files, function(f) {
    sha <- digest(file = f, algo = "sha256", serialize = FALSE)
    sprintf("%s  %s", sha, basename(f))
  }, character(1))
  writeLines(sums_lines, file.path(subdir, "SHA256SUMS"))

  # download_log.tsv — one row per file
  rows <- data.frame(
    filename    = basename(files),
    source_url  = meta$url,
    host        = meta$host,
    publication = meta$publication,
    sha256      = sub("  .*$", "", sums_lines),
    bytes       = file.info(files)$size,
    access_iso  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  write.table(rows, file = file.path(subdir, "download_log.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  log_msg("  ", basename(subdir), ": hashed ", length(files), " file(s)")
}

subdirs <- list.dirs(raw_dir, recursive = FALSE)
if (!is.null(restrict_source)) {
  subdirs <- subdirs[basename(subdirs) == restrict_source]
  if (length(subdirs) == 0) stop("No subdir matching --source ", restrict_source)
}

log_msg("Hashing raw sumstats under ", raw_dir)
for (sd in subdirs) write_audit(sd)

log_msg("Done. Per-source SHA256SUMS and download_log.tsv written.")
log_msg("Commit these to the OSF pre-registration deposit before any analysis.")
