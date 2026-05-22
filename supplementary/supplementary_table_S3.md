# Supplementary Table S3: ADNI cascade arm (Amendment 5)

## ADNI cascade arm (Amendment 5): cohort details + per-analysis full results

**Source:** Post-hoc biological-validation arm under OSF Amendment 5 (filed 2026-05-22T10:29:53 UTC at OSF parent project [osf.io/7u8a9](https://osf.io/7u8a9/); SHA256 57f00825…22692f77); model specification frozen on the public reproducibility repository before any cascade regression results were viewed. Data source: ADNI Omni2.5M genotypes and the ADNIMERGE2 v0.1.1 R package (UCBERKELEY_TAU_6MM, UCBERKELEY_TAUPVC_6MM, UCBERKELEY_AMY_6MM, APOERES, PTDEMOG, DXSUM tables). PRS-CS posterior weights (Ali 2023 multi-ethnic; φ = 10⁻²; seed = 81; 1,119,894 SNPs total; 486,274 intersect Omni2.5M — 43% coverage, above the asymptotic-coverage threshold for PGS performance). PLINK 1.9 --score sum for PGS; PLINK 1.9 --pca 10 on LD-pruned (200 / 50 / 0.2; MAF > 0.05; HWE > 10⁻⁶) variants for ancestry control.

**Cohort assembly:**

| Cohort | Subjects | Scans |
|---|---:|---:|
| ADNI Omni2.5M genotyped (total) | 812 | — |
| Berkeley FTP tau-PET subjects (any genotype status) | 1,204 | 2,171 |
| Berkeley amyloid-PET subjects (any tracer; any geno status) | 2,177 | 4,582 |
| Operational: FTP tau-PET ∩ Omni2.5M | **218** | 465 |
| Operational: amyloid-PET ∩ Omni2.5M | 719 | 1,971 |
| Both modalities ∩ Omni2.5M (within-subject path) | 217 | — |
| Longitudinal scans per subject in FTP cohort (range) | 1–7 (median 2) | — |

**Pre-registered primary cascade model (frozen on the public reproducibility repository before execution):**

`META_TEMPORAL_SUVR ~ amyloid_PGS_z + age_at_scan + sex + APOE_e4_carrier + PC1 + PC2 + PC3 + PC4`

estimated by ordinary least squares with the PGS expressed in standard-deviation units, on the cross-sectional first-scan-per-subject dataset.

**Primary + six pre-specified sensitivity results (full coefficients):**

| Analysis | n (subj / scans) | β_PGS | SE_PGS | 95% CI | p | R² | adj R² |
|---|---:|---:|---:|---|---|---:|---:|
| **Primary: tau META-TEMP SUVR (xs)** | 218 / 218 | **0.0810** | 0.0190 | 0.0436, 0.1183 | **2.98×10⁻⁵** | 0.200 | 0.169 |
| Sensitivity 1: tau META-TEMP SUVR (PVC) | 204 / 204 | 0.1379 | 0.0380 | 0.0632, 0.2127 | 3.56×10⁻⁴ | 0.179 | 0.145 |
| Sensitivity 2: amyloid centiloid (cross-modality) | 719 / 719 | 33.27 | 1.08 | 31.15, 35.39 | 1.28×10⁻¹³³ | 0.663 | 0.659 |
| Sensitivity 3: APOE-ε4 carriers | 72 / 72 | 0.1474 | 0.0448 | 0.0580, 0.2368 | 1.62×10⁻³ | 0.218 | 0.133 |
| Sensitivity 4: APOE-ε4 non-carriers | 146 / 146 | 0.0435 | 0.0151 | 0.0137, 0.0733 | 4.52×10⁻³ | 0.105 | 0.060 |
| Sensitivity 5: EUR-only (White, Non-Hisp) | 195 / 195 | 0.0964 | 0.0219 | 0.0531, 0.1397 | 1.85×10⁻⁵ | 0.208 | 0.174 |
| Longitudinal LMM: tau (random subj intercept) | 218 / 465 | 0.0900 | 0.0191 | 0.0526, 0.1274 | 4.40×10⁻⁶ | — | — |
| Longitudinal LMM: amyloid centiloid | 719 / 1,971 | 33.83 | 1.16 | 31.55, 36.10 | 1.49×10⁻¹²³ | — | — |
| **PGS × APOE-ε4 interaction term** | 218 / 218 | **0.1170** | 0.0377 | 0.0426, 0.1913 | **2.19×10⁻³** | 0.236 | 0.203 |

Primary-model concordant covariate effects: APOE-ε4 carrier β = 0.130 (p = 2.3×10⁻³); age at scan β = 0.0012 (p = 0.65); sex (male) β = −0.076 (p = 0.05). Ancestry PCs were not individually significant. Longitudinal LMM uses `lme4` + `lmerTest`; Satterthwaite df.

**Interaction model spec:** `META_TEMPORAL_SUVR ~ pgs_z * APOE_e4_carrier + age_at_scan + sex_m + PC1 + PC2 + PC3 + PC4` (n=218; OLS). The interaction row above is the cross-product coefficient; both main effects are also in the model. In this same model the pgs_z main effect is β = 0.0328 (SE 0.0242, p = 0.18) — this is the non-carrier-stratum slope (since APOE_e4_carrier is binary 0/1); the carrier-stratum slope is 0.0328 + 0.1170 = 0.1498, consistent with the carrier-only Sensitivity 3 row (β = 0.147). The significant interaction (p = 2.19×10⁻³) tells us the polygenic-amyloid effect on tau-PET is approximately 4.6-fold larger in APOE-ε4 carriers than in non-carriers; the non-carrier effect is independently significant in the stratified Sensitivity 4 model (β = 0.044, p = 4.52×10⁻³). APOE-ε4 amplifies, but does not create, the polygenic-amyloid → tau-PET relationship.

**Pre-specified decision rule (frozen in Amendment 5):** "Supportive of cascade hypothesis if β_PGS > 0 with p < 0.05 in the primary model AND consistent direction (β_PGS > 0) in ≥4 of 6 sensitivity analyses." Observed: primary positive at p = 2.98×10⁻⁵; all 6 sensitivity analyses positive and consistent in direction. The decision rule is satisfied.

**DUA-compliance note:** ADNI individual-level data are not redistributed; the values above are group-level summary statistics only. Researchers can apply for ADNI access at adni.loni.usc.edu and reproduce these results using the cascade-analysis script in the public reproducibility deposit (see main manuscript Code Availability).
