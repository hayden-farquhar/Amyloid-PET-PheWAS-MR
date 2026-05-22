# R/20_adni_cascade_pgs.R
#
# ADNI tau-PET cascade: tests whether the pre-registered amyloid PGS
# (PRS-CS-shrunk Ali 2023 multi-ethnic weights, Pivot A primary
# inferential layer) predicts AD-region tau-PET burden in ADNI,
# providing biological-validation evidence for the manuscript's
# chr1q32.2/CR1 mediator finding.
#
# Pre-registered primary endpoint: META_TEMPORAL_SUVR
# (Berkeley FTP pipeline, 6mm; first scan per subject)
#
# Pre-registered model:
#   META_TEMPORAL_SUVR ~ amyloid_PGS + age_at_scan + sex
#                      + APOE_e4_carrier + PC1 + PC2 + PC3 + PC4
#
# Sensitivity / secondary:
#   - Tau-PVC (META_TEMPORAL_SUVR from UCBERKELEY_TAUPVC_6MM)
#   - Cross-modality: CENTILOIDS ~ amyloid_PGS (replicates Ali in ADNI)
#   - APOE-stratified (carrier vs non-carrier)
#   - Longitudinal LMM (all scans, random subject intercept)
#   - EUR-ancestry restricted (matches PRS-CS LD reference)
#
# Data sources (all DUA-restricted — never copy out of data/raw/adni/):
#   - data/raw/adni/WGS_Omni25_BIN_wo_ConsentsIssues.{bed,bim,fam}
#   - data/raw/adni/adnimerge2/ADNIMERGE2/data/*.rda
#
# PRS-CS weights (already-computed):
#   - data/interim/prscs/amyloid_pgs_pst_eff_a1_b0.5_phi1e-02_chr*.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(lme4)
  library(lmerTest)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PROJ      <- Sys.getenv("AMYLOID_PHEWAS_MR_ROOT", ".")
ADNI_GENO <- file.path(PROJ, "data/raw/adni/WGS_Omni25_BIN_wo_ConsentsIssues")
ADNI_RDA  <- file.path(PROJ, "data/raw/adni/adnimerge2/ADNIMERGE2/data")
PRSCS_DIR <- file.path(PROJ, "data/interim/prscs")
OUT_DIR   <- file.path(PROJ, "data/processed/adni_cascade")
TMP_DIR   <- file.path(PROJ, "data/interim/adni_cascade")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1 — Build PLINK score file from PRS-CS posterior
# (one file, all chromosomes, three columns: SNP A1 BETA)
# -----------------------------------------------------------------------------
score_file <- file.path(TMP_DIR, "amyloid_prscs_score.tsv")
if (!file.exists(score_file)) {
  message("Step 1: combining PRS-CS posterior weights into a PLINK score file...")
  prscs_files <- list.files(PRSCS_DIR,
                            pattern = "amyloid_pgs_pst_eff_.*\\.txt$",
                            full.names = TRUE)
  prscs <- bind_rows(lapply(prscs_files, function(f) {
    read_tsv(f, col_names = c("CHR","SNP","BP","A1","A2","BETA"),
             col_types = "iciccd", progress = FALSE)
  }))
  write_tsv(prscs |> select(SNP, A1, BETA), score_file, col_names = FALSE)
  message("  wrote ", nrow(prscs), " variants to ", basename(score_file))
} else {
  message("Step 1 skipped: ", basename(score_file), " exists.")
}

# -----------------------------------------------------------------------------
# Step 2 — PLINK --score on the ADNI Omni2.5M
# -----------------------------------------------------------------------------
score_out <- file.path(TMP_DIR, "adni_amyloid_pgs")
if (!file.exists(paste0(score_out, ".profile"))) {
  message("Step 2: running PLINK --score on ADNI Omni2.5M...")
  cmd <- sprintf(
    "plink --bfile %s --score %s sum --out %s",
    shQuote(ADNI_GENO), shQuote(score_file), shQuote(score_out)
  )
  message("  ", cmd)
  system(cmd, intern = FALSE)
} else {
  message("Step 2 skipped: ", basename(score_out), ".profile exists.")
}
pgs <- read.table(paste0(score_out, ".profile"), header = TRUE) |>
  as_tibble() |>
  transmute(PTID = IID, n_scored = CNT, n_alleles = CNT2, pgs_sum = SCORESUM)
message("  PGS computed for ", nrow(pgs), " subjects; median CNT = ",
        median(pgs$n_scored))

# -----------------------------------------------------------------------------
# Step 3 — Ancestry PCs via PLINK --pca (top 10)
# -----------------------------------------------------------------------------
pca_out <- file.path(TMP_DIR, "adni_pca")
if (!file.exists(paste0(pca_out, ".eigenvec"))) {
  message("Step 3: PCA on ADNI Omni2.5M (LD-pruned, MAF>0.05, HWE>1e-6)...")
  prune_out <- file.path(TMP_DIR, "adni_prune")
  if (!file.exists(paste0(prune_out, ".prune.in"))) {
    system(sprintf(
      "plink --bfile %s --maf 0.05 --hwe 1e-6 --indep-pairwise 200 50 0.2 --out %s",
      shQuote(ADNI_GENO), shQuote(prune_out)
    ))
  }
  system(sprintf(
    "plink --bfile %s --extract %s.prune.in --pca 10 header --out %s",
    shQuote(ADNI_GENO), shQuote(prune_out), shQuote(pca_out)
  ))
}
pcs <- read.table(paste0(pca_out, ".eigenvec"), header = TRUE) |>
  as_tibble() |>
  select(PTID = IID, starts_with("PC"))
message("  PCs loaded for ", nrow(pcs), " subjects")

# -----------------------------------------------------------------------------
# Step 4 — Load ADNIMERGE2 phenotype + covariate tables
# -----------------------------------------------------------------------------
message("Step 4: loading ADNIMERGE2 tables...")
load_rda <- function(name) {
  env <- new.env()
  load(file.path(ADNI_RDA, paste0(name, ".rda")), envir = env)
  get(name, envir = env)
}
tau    <- load_rda("UCBERKELEY_TAU_6MM")     |> as_tibble()
taupvc <- load_rda("UCBERKELEY_TAUPVC_6MM")  |> as_tibble()
amy    <- load_rda("UCBERKELEY_AMY_6MM")     |> as_tibble()
apoe   <- load_rda("APOERES")                |> as_tibble()
dem    <- load_rda("PTDEMOG")                |> as_tibble()
dx     <- load_rda("DXSUM")                  |> as_tibble()

# APOE-ε4 carrier status (GENOTYPE field is e.g. "3/3", "3/4", "4/4")
apoe_pt <- apoe |>
  filter(!is.na(GENOTYPE)) |>
  group_by(PTID) |> slice(1) |> ungroup() |>
  mutate(APOE_e4_n      = stringr::str_count(GENOTYPE, "4"),
         APOE_e4_carrier = as.integer(APOE_e4_n >= 1)) |>
  select(PTID, APOE_GENOTYPE = GENOTYPE, APOE_e4_n, APOE_e4_carrier)

# Demographics — one row per PTID. PTDOB is character "MM/YYYY";
# PTDOBYY is the year as numeric with identical coverage. Use PTDOBYY.
dem_pt <- dem |>
  group_by(PTID) |> slice(1) |> ungroup() |>
  mutate(birth_year = as.integer(PTDOBYY)) |>
  select(PTID, sex = PTGENDER, birth_year,
         PTRACCAT, PTETHCAT, PTEDUCAT)

# Diagnosis at first available visit (uses scan-date alignment downstream
# rather than baseline-only)
dx_first <- dx |>
  filter(!is.na(DIAGNOSIS)) |>
  arrange(PTID, EXAMDATE) |>
  group_by(PTID) |> slice(1) |> ungroup() |>
  select(PTID, DIAGNOSIS_first = DIAGNOSIS)

# -----------------------------------------------------------------------------
# Step 5 — Build the analysis frame
# -----------------------------------------------------------------------------
build_pheno <- function(suvr_table, tracer_filter, outcome_col,
                        label = "tau_meta_temp") {
  out <- suvr_table |>
    filter(TRACER %in% tracer_filter, !is.na(.data[[outcome_col]])) |>
    mutate(scan_date = as.Date(SCANDATE)) |>
    select(PTID, scan_date, outcome = all_of(outcome_col), TRACER) |>
    inner_join(pgs,      by = "PTID") |>
    inner_join(pcs,      by = "PTID") |>
    inner_join(apoe_pt,  by = "PTID") |>
    inner_join(dem_pt,   by = "PTID") |>
    left_join (dx_first, by = "PTID") |>
    mutate(age_at_scan = as.numeric(format(scan_date, "%Y")) - birth_year,
           sex_m       = as.integer(sex == "Male"),
           pgs_z       = scale(pgs_sum)[, 1])
  attr(out, "outcome_label") <- label
  out
}

tau_long    <- build_pheno(tau,    "FTP", "META_TEMPORAL_SUVR", "tau_meta_temp")
taupvc_long <- build_pheno(taupvc, "FTP", "META_TEMPORAL_SUVR", "tau_meta_temp_pvc")
amy_long    <- build_pheno(amy, c("FBB","FBP","NAV"), "CENTILOIDS",      "amy_centiloid")

first_per_subject <- function(d) d |> arrange(PTID, scan_date) |>
  group_by(PTID) |> slice(1) |> ungroup()

tau_xs    <- first_per_subject(tau_long)
taupvc_xs <- first_per_subject(taupvc_long)
amy_xs    <- first_per_subject(amy_long)

message("\n=== Analysis frame sizes (subjects with all covariates) ===")
message("  tau-PET FTP cross-sectional:    n = ", nrow(tau_xs))
message("  tau-PET FTP longitudinal scans: n = ", nrow(tau_long))
message("  tau-PVC cross-sectional:        n = ", nrow(taupvc_xs))
message("  amyloid-PET cross-sectional:    n = ", nrow(amy_xs))
message("  tau ∩ amyloid (within-subject): n = ",
        length(intersect(tau_xs$PTID, amy_xs$PTID)))

# -----------------------------------------------------------------------------
# Step 6 — Pre-registered primary cascade test
# -----------------------------------------------------------------------------
message("\n=== Step 6: pre-registered primary cascade ===")
formula_primary <- outcome ~ pgs_z + age_at_scan + sex_m +
                             APOE_e4_carrier + PC1 + PC2 + PC3 + PC4

fit_primary <- lm(formula_primary, data = tau_xs)
print(summary(fit_primary)$coefficients)
ci_primary <- confint(fit_primary, "pgs_z")
message(sprintf(
  "  PGS effect on META_TEMPORAL_SUVR (cross-sectional, n=%d): β = %.4f (95%% CI %.4f, %.4f), p = %.3g",
  nrow(tau_xs),
  coef(fit_primary)["pgs_z"], ci_primary[1], ci_primary[2],
  summary(fit_primary)$coefficients["pgs_z", "Pr(>|t|)"]
))

# -----------------------------------------------------------------------------
# Step 7 — Sensitivity analyses (results saved to one tidy table)
# -----------------------------------------------------------------------------
message("\n=== Step 7: sensitivity analyses ===")

run_lm <- function(d, label, formula = formula_primary) {
  fit <- lm(formula, data = d)
  s   <- summary(fit)$coefficients["pgs_z", ]
  ci  <- confint(fit, "pgs_z")
  tibble(
    analysis    = label,
    n           = nrow(d),
    beta_pgs    = s["Estimate"],
    se_pgs      = s["Std. Error"],
    ci_lo       = ci[1],
    ci_hi       = ci[2],
    p_pgs       = s["Pr(>|t|)"],
    r2          = summary(fit)$r.squared,
    adj_r2      = summary(fit)$adj.r.squared
  )
}

run_lmm <- function(d, label) {
  fit <- lmer(outcome ~ pgs_z + age_at_scan + sex_m + APOE_e4_carrier +
              PC1 + PC2 + PC3 + PC4 + (1 | PTID), data = d)
  s   <- summary(fit)$coefficients["pgs_z", ]
  ci  <- confint(fit, "pgs_z", method = "Wald")
  tibble(
    analysis    = label,
    n           = nrow(d),
    n_subjects  = n_distinct(d$PTID),
    beta_pgs    = s["Estimate"],
    se_pgs      = s["Std. Error"],
    ci_lo       = ci[1],
    ci_hi       = ci[2],
    p_pgs       = s["Pr(>|t|)"]
  )
}

results <- bind_rows(
  run_lm(tau_xs,                                 "primary_tau_meta_temp_xs"),
  run_lm(taupvc_xs,                              "sens_tau_meta_temp_pvc_xs"),
  run_lm(amy_xs,                                 "sens_amy_centiloid_xs"),
  run_lm(tau_xs |> filter(APOE_e4_carrier == 1), "sens_tau_apoe_e4_carrier"),
  run_lm(tau_xs |> filter(APOE_e4_carrier == 0), "sens_tau_apoe_e4_noncarrier"),
  run_lm(tau_xs |> filter(PTRACCAT == "White" & PTETHCAT == "Not Hispanic or Latino"),
         "sens_tau_eur_only_xs"),
)

# Longitudinal LMM rows (different schema — keep separate then row-bind)
lmm_rows <- bind_rows(
  run_lmm(tau_long, "sens_tau_long_lmm"),
  run_lmm(amy_long, "sens_amy_long_lmm")
)

print(results)
cat("\n--- Longitudinal LMM ---\n")
print(lmm_rows)

# -----------------------------------------------------------------------------
# Step 8 — Persist
# -----------------------------------------------------------------------------
write_tsv(results,  file.path(OUT_DIR, "adni_cascade_lm_results.tsv"))
write_tsv(lmm_rows, file.path(OUT_DIR, "adni_cascade_lmm_results.tsv"))
saveRDS(fit_primary, file.path(OUT_DIR, "fit_primary_lm.rds"))
message("\nWrote results to ", OUT_DIR)
message("  - adni_cascade_lm_results.tsv")
message("  - adni_cascade_lmm_results.tsv")
message("  - fit_primary_lm.rds")
