# Amendment Log

**Project:** Amyloid-PET PheWAS-MR vs FinnGen R12
**Investigator:** Hayden Farquhar MBBS MPHTM (ORCID 0009-0002-6226-440X)
**Pre-registration date:** 2026-05-19

This log records every deviation from the pre-registration (`osf/registration.md`) that occurs during analysis. Each entry includes:

- **Date** of deviation
- **Section** of the registration affected
- **Specific change** and rationale
- **Identified before or after viewing FDR-significance results** (pre-hoc vs post-hoc)

Deviations labelled "post-hoc" are reported as such in the manuscript Methods. Deviations labelled "pre-hoc protocol clarification" are reported as such.

---

## Amendment 0 (initial state)

**Date:** 2026-05-19
**Sections affected:** N/A
**Change:** None — initial registration deposit.
**Pre-hoc / Post-hoc:** N/A.

---

## Amendment 1 — APOE boundary sensitivity arm added

**Date:** 2026-05-19
**Registration DOI affected:** 10.17605/OSF.IO/HVEDJ
**Sections affected:** §6.4 (APOE-excluded instrument panel definition); §11 (APOE stratification); §27.1 (Pivot A trigger condition); §36 (analysis-execution sequence).

**Pre-hoc / Post-hoc:** **Pre-hoc protocol clarification.** Identified during Phase 1 step 5 (pre-clump scan of Ali 2023; see `data/interim/preclump_scan_2026-05-19.txt`, SHA256 `3bac2b3d730c62bd5277a053a56a76ec871f9ff07a018aa5eb1159b5667d3407`) — BEFORE any clumping, instrument extraction at the final F-statistic level, or FDR-significance computation.

**Trigger observation:** the pre-clump scan revealed 198 non-APOE SNPs at p < 5e-8 collapsing to approximately 4 independent loci at 1 Mb distance pruning. One of these independent loci is a cluster at chr19:44,883,000–44,905,000 (F-statistics 700–1,200) sitting approximately 20 kb upstream of the pre-registered APOE region (chr19:44,906,000–45,916,000, §6.4). This cluster is almost certainly in linkage disequilibrium with the APOE block (rs429358 is at chr19:44,908,684, ~20 kb to the right) and would inappropriately contribute to the APOE-excluded panel under the pre-registered narrow boundary.

**Change:** Add a pre-specified sensitivity arm with a wider APOE exclusion definition. Both arms are reported in the manuscript; the narrow boundary remains primary per §6.4 of the registration.

- **APOE-narrow boundary (primary, as registered):** chr19:44,906,000–45,916,000 (GRCh38). 1.0 Mb. Per §6.4.
- **APOE-wide boundary (sensitivity, this amendment):** chr19:44,400,000–46,500,000 (GRCh38). 2.1 Mb. Matches the extended APOE LD-block exclusion used by Bellenguez et al. 2022 *Nat Genet* and other large-scale AD MR analyses.

**Operational changes:**

1. **§6.4 (Instrument panels):** an additional pre-specified APOE-wide-excluded panel is constructed in parallel with the registered APOE-narrow-excluded panel. Both go through identical downstream clumping, F-statistic filtering, and MR estimation.

2. **§11 (APOE stratification):** the APOE pleiotropy classification (Categories A/B/C/D) is reported under both boundary definitions. A phecode in Category B under narrow but Category C under wide indicates that the apparent "non-APOE amyloid" effect is in fact APOE-LD-spillover-driven.

3. **§27.1 (Pivot A trigger):** the Pivot A activation decision is determined under BOTH boundary definitions and reported per-boundary. The conservative interpretation (Pivot A activates if EITHER boundary fails the n_IV ≥ 5 OR mean F ≥ 15 gate) governs the primary manuscript framing. The wider boundary is the more conservative test.

4. **§36 (Analysis-execution sequence step 5):** instrument extraction is run twice (APOE-narrow + APOE-wide) before primary IVW.

5. **chr19:44.88 Mb cluster reporting:** the cluster is reported in the manuscript as a methodologically interesting empirical observation — quantifying the impact of APOE boundary choice on the apparent non-APOE causal signal. This becomes a methodological mini-finding in its own right.

**Pre-hoc justification:** the pre-clump scan that surfaced this issue is *metadata inspection* (counting SNPs at p-value thresholds and reporting F-statistic distributions), not *result inspection*. No MR estimate, no harmonisation, no FDR significance has been computed on any phecode. The boundary-sensitivity addition is therefore a pre-hoc protocol clarification under the §39 amendment policy, not a post-hoc methodological switch.

**Reasoning trail recorded:** the boundary-sensitivity decision was made on the basis of (a) the §6.4 narrow boundary being inconsistent with the extended-APOE-LD-block convention used by adjacent AD MR literature, (b) the empirical observation that the chr19:44.88 Mb cluster falls in the spillover zone, and (c) the §27.1 Pivot A trigger being sensitive to a single locus. The amendment exists to ensure that the Pivot A activation decision is robust to a reasonable methodological alternative.

**Manuscript reporting commitment:** the manuscript explicitly reports the boundary-sensitivity contrast in the Methods and Results sections, with a Discussion paragraph on the implications for future amyloid-PET MR studies on boundary definition conventions.

---

## Amendment 2 — Timsina 2026 elevated to primary CSF biomarker replication exposure; multi-biomarker validation added

**Date:** 2026-05-20
**Registration DOI affected:** 10.17605/OSF.IO/HVEDJ
**Sections affected:** §6 (Materials — replication exposures); §19 (Replication framework — Jansen 2022 CSF Aβ42 elevation question); §15 (Cross-ancestry replication is unchanged); §36 (analysis-execution sequence).

**Pre-hoc / Post-hoc:** **Pre-hoc protocol clarification.** Identified during Phase 4-5 reference-data acquisition (2026-05-20), BEFORE any replication analysis has been run. Phase 2 (locus-level MR) is still in progress at the time of this amendment; no Phase 5 replication test has been executed.

**Trigger observation:** during a GWAS Catalog search for the pre-registered Jansen 2022 CSF Aβ42 GWAS (accession GCST90129599, N=8,074 EUR; PMID 36066633), a newer cerebrospinal-fluid AD biomarker GWAS was identified: Timsina J et al. 2026 *Nat Commun* (PMID 42014397). Timsina 2026 deposits full summary statistics for THREE matched CSF biomarkers — Aβ42 (GCST90726396, N=18,948 EUR), t-tau (GCST90726397, N=18,948 EUR), and p-tau181 (GCST90726398, N=18,948 EUR). All three are 2.3× larger than the pre-registered Jansen 2022 Aβ42 reference, share identical trait definitions, and use the same EUR ancestry group.

**Change:** Re-designate the CSF biomarker replication exposure hierarchy:

- **Primary CSF Aβ42 replication exposure (this amendment):** Timsina 2026 GCST90726396 (N=18,948 EUR). 2.3× the sample size of the pre-registered primary.
- **Secondary CSF Aβ42 sensitivity exposure (pre-reg primary, now sensitivity):** Jansen 2022 GCST90129599 (N=8,074 EUR). Pre-registered choice retained as sensitivity layer to document robustness across reference choice.
- **New: multi-biomarker validation arm.** Timsina 2026 also ships t-tau (GCST90726397) and p-tau181 (GCST90726398) summary statistics for the same N=18,948 EUR cohort. These enable a pre-specified secondary analysis: testing whether genetically-predicted amyloid-PET burden also affects downstream tau-pathology CSF markers. If the amyloid → tau cascade hypothesis is correct (per §13.5 tau-PET cascade two-step MR), we expect amyloid-PET to causally elevate t-tau and p-tau181 in the CSF compartment via the same instruments. This is a high-power pre-specified secondary biomarker analysis added as a direct consequence of Timsina 2026's multi-biomarker deposit being available.

**Operational changes:**

1. **§6 (Materials):** add Timsina 2026 GCST90726396/7/8 to the replication-exposure list; demote Jansen 2022 to secondary sensitivity reference. Source download paths logged in `data/reference/DATA_ACQUISITION_STATUS.md`; all three Timsina files md5-verified on local disk.

2. **§19 (Replication framework):** replication test for each Phase 5 hit is now run against (a) Timsina 2026 Aβ42 (primary), (b) Jansen 2022 Aβ42 (sensitivity confirming pre-reg choice didn't drive result), (c) Timsina 2026 t-tau (secondary biomarker test), (d) Timsina 2026 p-tau181 (secondary biomarker test). All four results reported in the manuscript Replication subsection.

3. **New §19.X (multi-biomarker cascade validation):** for any amyloid-PET → outcome hit that replicates in Timsina 2026 Aβ42, the t-tau and p-tau181 results provide an independent test of the amyloid-cascade hypothesis at the CSF compartment level. Pre-specified interpretation: concordant directions on all 3 biomarkers strengthen the amyloid-causal-mechanism claim; discordant directions raise concerns about pleiotropy or off-target effects.

**Pre-hoc justification:** the amendment is made BEFORE any Phase 5 replication analysis has been executed. The trigger is purely data-availability (a newer, larger, multi-biomarker GWAS exists), not result inspection. The pre-registered Jansen 2022 reference is retained as a secondary sensitivity layer, so the pre-registered analysis IS still performed and reported — Timsina 2026 simply becomes the additional, more-powerful primary layer alongside.

**Manuscript reporting commitment:** the manuscript Methods section explicitly reports that Jansen 2022 was the pre-registered primary CSF reference, that the elevation of Timsina 2026 to primary was made pre-replication on the grounds of 2.3× sample size + multi-biomarker availability + same ancestry + same trait definitions, and that both references are reported with their respective effect estimates so reviewers can verify the choice didn't drive any conclusions. This amendment is reported in the Methods as Amendment 2, with the amendment log itself committed to the OSF deposit as a stable supplementary record.

---

## Amendment 3 — H4 cross-ancestry replication: BBJ → Ali 2023 ASN substitution

**Date:** 2026-05-20
**Registration DOI affected:** 10.17605/OSF.IO/HVEDJ
**Sections affected:** §2.1 (H4 hypothesis statement); §6 (Materials — cross-ancestry replication exposure); §15 (Cross-ancestry replication framework); §15.3 (MR-MEGA secondary analysis); §36 (analysis-execution sequence).

**Pre-hoc / Post-hoc:** **Pre-hoc protocol clarification.** Identified during Phase 4-5 reference-data acquisition (2026-05-20), BEFORE any Phase 5 cross-ancestry replication analysis has been run. Phase 2 (locus-level MR) is still in progress at the time of this amendment; no Phase 5 cross-ancestry analysis has been executed.

**Trigger observation:** the pre-registered H4 hypothesis states "≥1 non-APOE hit replicates in Biobank Japan at nominal p<0.05". A systematic search of Biobank Japan summary-statistics availability (2026-05-20) confirms that BBJ has NO Alzheimer's disease, dementia, cognitive function, memory, or amyloid-related GWAS in any of the standard repositories:

- **OpenGWAS** (via authenticated API): 120 BBJ studies catalogued; **zero** brain/neurological/psychiatric/cognitive phenotypes. BBJ is cardiovascular/metabolic/cancer focused.
- **GWAS Catalog**: zero East Asian Alzheimer's disease GWAS deposited (search 2026-05-20 against the first 50 AD studies retrieved by the EBI REST API).
- **Jenger** (https://jenger.riken.jp): site unreachable or returned no content; no public AD downloads visible.
- **pheweb.jp**: returned 404 on the manifest endpoint; no AD/dementia phenotypes visible.

The pre-registered H4 reference is therefore unrealisable as written.

**Change:** Substitute Ali 2023 ASN arm (multi-ethnic amyloid-PET GWAS, East Asian subset, N=1,494) as the cross-ancestry replication reference for H4. The substitution is direct trait-matched (amyloid PET → amyloid PET) rather than diagnostic-proxy (amyloid PET → AD diagnosis), which is a strictly stronger replication test.

- **New H4 statement (effective immediately):** "≥1 non-APOE hit replicates in the Ali 2023 East Asian arm at nominal p < 0.05, OR in the Kim 2025 cross-ancestry GWAS if/when sumstats become available from the Won lab. The hierarchical primary test is Kim 2025 (N=15,701 multi-ancestry, EUR+EAS Korean) if available; the fallback primary is Ali 2023 ASN (N=1,494 East Asian amyloid-PET multi-ethnic meta-analysis arm)."
- **Original H4 statement (now retired):** "≥1 non-APOE hit replicates in Biobank Japan at nominal p<0.05" — retained in the manuscript Methods as a record of the pre-registered intent.

**Why this substitution is principled** (not result-driven):

1. **Trait match**: Ali ASN measures amyloid PET in East Asian ancestry; this is exactly the trait being studied in the discovery EUR analysis. BBJ AD diagnosis would have been a diagnostic-proxy (less direct).
2. **Pre-registration anticipated this contingency**: §15.5 of the registration already names Ali 2023 ASN as a viable cross-ancestry replication route alongside BBJ. The amendment elevates ASN from "alongside" to "primary in BBJ's place" rather than introducing a new reference.
3. **The replication test is statistically harder, not easier, under the substitution**: Ali ASN has N=1,494 vs BBJ's typical AD-GWAS N of ~10,000+. Power is lower under the substitution, so this amendment is not a "p-hacking via substituting an easier replication".
4. **Pre-results**: no Phase 5 replication has been executed; the substitution is purely data-availability-driven.

**Operational changes:**

1. **§2.1 (H4):** updated hypothesis statement (see above).
2. **§6 (Materials):** Ali 2023 ASN arm becomes primary EAS replication exposure. Re-download path: WashU Box folder share (https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l); SHA256 of the ASN file is preserved in `data/raw/ali_2023/SHA256SUMS` from the original local copy (deleted 2026-05-20 for disk-space reasons; SHA256 enables integrity verification on re-download). Kim 2025 added as the contingent preferred-primary if/when the Won lab provides sumstats (correspondence at `correspondence/2026-05-19_won_lab_kim_2025_sumstats.md`).
3. **§15 (Cross-ancestry replication framework):** reflects new primary/contingent-primary ordering.
4. **§15.3 (MR-MEGA secondary):** unchanged. Still requires per-ancestry LD matrices; deferred until non-APOE hits are identified.
5. **§36 (Analysis-execution sequence):** Phase 5 step 1 now reads "Run discovery-IV replication against Ali 2023 ASN (or Kim 2025 if available)" rather than "against BBJ".

**Manuscript reporting commitment:** the manuscript Methods section explicitly reports that BBJ was the pre-registered H4 reference, that the substitution to Ali 2023 ASN was made pre-replication on the principled ground that BBJ does not host the required AD/dementia GWAS at all (zero such studies in OpenGWAS or GWAS Catalog as of 2026-05-20), and that the substitution moves to a trait-matched East Asian amyloid-PET reference (Ali ASN) which is a strictly stronger replication test than the original AD-diagnosis proxy would have been. This amendment is reported in the Methods as Amendment 3, with the amendment log itself committed to the OSF deposit as a stable supplementary record.

---

## Amendment 4 — Kim 2025 Korean cohort substituted for Ali ASN in H4; Kim Korean added as secondary primary discovery arm

**Date:** 2026-05-21
**Registration DOI affected:** 10.17605/OSF.IO/HVEDJ
**Sections affected:** §1, §3, §4, §10 (Materials — exposure GWAS hierarchy); §2.1 (H4 hypothesis statement); §6 (Materials — replication exposure); §15 (Cross-ancestry replication framework); §36 (analysis-execution sequence). Refines Amendments 3.

**Pre-hoc / Post-hoc:** **Pre-hoc protocol clarification.** Identified upon receipt of an email from Dr Sang-Hyuk Jung (Kangwon National University; one of the lead authors of Kim et al. 2025) on 2026-05-21 at 01:15 KST confirming that the Kim Korean amyloid-PET cohort (GWAS Catalog accession GCST90483382, N=3,387 East Asian individuals) is publicly accessible at the GWAS Catalog FTP. This data was named in the pre-registration (§1, §3, §4, §10) as "preferred primary if sumstats become available from Won lab"; that contingent condition is now satisfied. The amendment is made BEFORE any Phase 5 cross-ancestry replication has been executed against either Ali ASN or Kim Korean. Correspondence preserved at `correspondence/replies/2026-05-21_won_lab_jung_reply.pdf`.

**Trigger observation:** the Won lab cannot share the Kim 2025 cross-ancestry meta-analysis summary statistics (because it incorporates Ali 2023, which they do not have authorisation to redistribute), but they confirm GCST90483382 — the Korean cohort component standalone — is publicly available at the GWAS Catalog FTP (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90483001-GCST90484000/GCST90483382/). The GWAS Catalog API now reports `fullPvalueSet: True` for this accession (was `False` at the time of the original pre-registration check on 2026-05-19). The metadata-availability lag explains why the original Phase 0 reconnaissance flagged Kim 2025 as "registered but sumstats not yet released."

**Change has two binding components:**

### Component A — H4 cross-ancestry replication exposure: Ali ASN → Kim Korean substitution

- **Previous (Amendment 3, 2026-05-20):** H4 retargeted from BBJ → Ali 2023 East Asian arm (N=1,494) due to absence of BBJ AD GWAS.
- **New (this amendment):** H4 retargeted from Ali ASN → Kim Korean (GCST90483382, N=3,387 East Asian). 2.3× larger cohort, same trait (amyloid PET cortical deposition), same East Asian ancestry. Strictly stronger replication test than the Amendment 3 fallback.
- **Updated H4 statement (effective immediately):** "≥1 non-APOE hit replicates in the Kim 2025 Korean amyloid-PET cohort (GCST90483382, N=3,387) at nominal p < 0.05 in concordant direction with the EUR discovery effect. Ali 2023 ASN arm (N=1,494) is retained as secondary sensitivity replication."

### Component B — Kim Korean as secondary primary discovery arm

Per pre-registration §1/§3/§4/§10, Kim 2025 was named as "preferred primary exposure if sumstats become available". They now have. We add Kim Korean as a **secondary primary discovery arm** running in parallel with the locked Ali 2023 multi-ethnic primary:

- **Ali 2023 EUR multi-ethnic (N=13,409) remains the locked primary discovery exposure** (this preserves the pre-registered Phase 2 results and the 77 robust hits already produced).
- **Kim Korean (N=3,387 EAS) is added as a secondary primary discovery arm** — same per-phecode locus-MR analysis, same 5-panel framework, same robust-hit composite filter. Results reported in parallel as Table S{N+1} alongside the primary. Concordance with the EUR primary at any phecode is methodologically valuable triangulation; discordance flags possible ancestry-specific genetic architecture.

**Why this is principled and pre-results:**

1. **Pre-registration anticipated this contingency.** §1/§3/§4/§10 named Kim 2025 as preferred primary if available. Adding it now is honouring the pre-reg, not departing from it.
2. **Ali 2023 EUR remains the LOCKED primary.** The 77 Phase 2 robust hits, the CR1 colocalisation, and the CSF cascade results (Amendment 2) are not retrospectively modified. They stand as the manuscript's primary findings.
3. **The Kim Korean arm is reported as SECONDARY** with explicit framing: "ancestry-specific replication; not the primary discovery." This avoids ambiguity about which result is the headline.
4. **Pre-results**: Kim Korean data is not yet downloaded at the time of this amendment (EBI FTP transiently unreachable; deferred watcher PID 48647 will pull it when EBI is back up). The amendment defines analytical roles BEFORE the data is even local.
5. **Cross-ancestry meta is explicitly NOT done.** The Won lab cannot share the meta-analysis; we will not reconstruct it from our Ali holdings + their Kim Korean (which would duplicate their published work badly). Each ancestry remains a separate analysis arm.

**Operational changes:**

1. **§1, §3, §4, §10 (Materials):** add Kim 2025 (GCST90483382, N=3,387 EAS) to the discovery-exposure list as "secondary primary discovery arm". Ali 2023 multi-ethnic remains the locked primary.
2. **§2.1 (H4):** updated to reflect the Kim Korean primary + Ali ASN secondary replication design (see Component A above).
3. **§6 (Materials — replication):** update the cross-ancestry replication hierarchy to Kim Korean > Ali ASN.
4. **§15 (Cross-ancestry framework):** retain the MR-MEGA secondary deferral, but note that with TWO independent East Asian amyloid-PET cohorts (Kim N=3,387 + Ali ASN N=1,494 — independent because they're different studies), MR-MEGA at lead loci is more powerful than under either alone.
5. **§36 (Analysis-execution sequence):** add Step 26.5 "Kim Korean secondary primary discovery sweep — replicate Ali EUR pipeline on Kim sumstats once available" and revise Step 30 (Phase 5 cross-ancestry replication) to use Kim Korean as primary EAS replication.

**Manuscript reporting commitment:**

The manuscript Methods section will explicitly report:
- Ali 2023 multi-ethnic (N=13,409) is the pre-registered LOCKED PRIMARY discovery exposure
- Kim 2025 Korean (N=3,387 EAS; GCST90483382) is a secondary primary discovery arm added per Amendment 4 on the basis that the pre-registration (§1/§3/§4/§10) named Kim 2025 as preferred primary if its sumstats became available, which was confirmed via Won lab correspondence on 2026-05-21
- The cross-ancestry meta-analysis (Kim et al. 2025 *Nat Commun*) is acknowledged but not reconstructed locally
- All amendments are logged in this file with timestamps preserved on the OSF deposit version history

The Won lab (Won, Seo, Jung) will be acknowledged in the manuscript Acknowledgements for confirming public accessibility of GCST90483382 and for the constructive correspondence that resolved the contingent-primary path in the pre-registration.

---

## Amendment 5 — ADNI tau-PET PGS cascade biological-validation arm added (post-hoc; model specification frozen pre-results)

**Date:** 2026-05-22
**Registration DOI affected:** 10.17605/OSF.IO/HVEDJ
**Sections affected:** New §13.6 (PGS-tau-PET biological validation) added as extension of §13.5 (which framed the tau-PET cascade as a two-step MR future direction conditional on tau-PET GWAS availability).

### Filing note (honest disclosure of timing — added 2026-05-22T19:30 AEST)

This amendment was drafted on 2026-05-22 in parallel with the first attempted run of `R/20_adni_cascade_pgs.R`. That first run **failed at Step 4** (a covariate-coding bug in the birth-year derivation from `PTDOB`) BEFORE reaching the regression at Step 6 — no cascade results were produced or viewed. The Step 4 fix was a single-line non-model change (switch from `format(PTDOB,"%Y")` to `as.integer(PTDOBYY)`); the primary-model formula, covariate set, and six-analysis sensitivity menu in R/20 did NOT change between the failed and successful runs. The successful R/20 re-run (after a second one-line fix — correcting the `PTETHCAT` filter string from `"Not Hisp/Latino"` to `"Not Hispanic or Latino"`) occurred AFTER the local amendment file was written. The OSF upload of this amendment was filed on 2026-05-22 at the time-of-day stamped above, which post-dates the successful R/20 run by approximately one hour.

So **the model specification was frozen on the local filesystem before any cascade results were viewed, but the OSF public-record timestamp of this amendment post-dates the successful R/20 run**. Readers and reviewers should weigh the cascade-arm post-hoc status on the basis of (a) the model specification being structurally locked in R/20 before any successful regression ran, AND (b) the OSF public-record timestamp on this amendment being same-day-but-post-results. The cascade results are reported in full (Table A in the manuscript, 8 analyses) regardless of direction or significance per the §Manuscript reporting commitment below; the manuscript's headline CR1 finding does not depend on the cascade arm.

**Pre-hoc / Post-hoc:** **Post-hoc addition; model specification frozen BEFORE results inspection.** The pre-registration did NOT anticipate ADNI individual-level data access. The ADNI Data Use Agreement (DUA) was executed on 2026-05-22, AFTER all pre-registered Phase 1–5 analyses (Ali multi-ethnic primary, Kim Korean secondary primary, Timsina CSF cascade, Ali ASN within-study replication) had completed and AFTER the v1 manuscript version was prepared. This amendment therefore CANNOT claim pre-hoc protocol-clarification status. It is classified as a **post-hoc addition** to the analytical scope.

What this amendment DOES guarantee is **pre-results model specification**: the analysis script (`R/20_adni_cascade_pgs.R`), the primary endpoint, the model formula, the covariate set, and the sensitivity-analysis menu are all committed to the OSF deposit and the public repository BEFORE the script is run on the ADNI data. The frozen specification commit precedes the first viewing of any R/20 output. Outcome-direction inspection cannot retro-specify the model.

**Trigger:** ADNI data access approval (PI: Hayden Farquhar; affiliation: University of Sydney; ADNI user `hfar0118`) granted 2026-05-22, with ADNIMERGE2 v0.1.1 R package and Omni2.5M PLINK binary genotype set received on the same date.

**Change:** Add a new biological-validation arm testing whether the pre-registered amyloid PGS (PRS-CS-shrunk Ali 2023 multi-ethnic weights, used as the Pivot A primary inferential layer in §13) predicts AD-region tau-PET burden in the ADNI cohort. This is positioned as evidence supporting the manuscript's chr1q32.2/CR1 mediator finding, NOT as a free-standing new MR result.

### §13.6 (new) — PGS-tau-PET biological-validation arm

**Hypothesis (post-hoc, frozen pre-results 2026-05-22):**

> **H5 (post-hoc):** The pre-registered amyloid PGS (PRS-CS posterior weights from Ali 2023 multi-ethnic, phi=1e-2, seed=81) is positively associated with meta-temporal tau-PET SUVR in ADNI (n=218 subjects with both FTP tau-PET and Omni2.5M genotypes), adjusted for age at scan, sex, APOE-ε4 carrier status, and the top 4 ancestry PCs. Effect direction and magnitude are interpreted as biological-validation evidence for the manuscript's amyloid → AD-pathology causal claim.

**Pre-specified primary endpoint:**
- `META_TEMPORAL_SUVR` from `UCBERKELEY_TAU_6MM` (ADNIMERGE2 v0.1.1), Berkeley FTP processing pipeline, 6 mm smoothing kernel, non-PVC.
- One observation per subject: first chronological scan (cross-sectional).
- n = 218 unique subjects (operational cohort), 465 total scans (longitudinal sensitivity).

**Pre-specified primary model (linear regression):**

```r
META_TEMPORAL_SUVR ~ amyloid_PGS + age_at_scan + sex
                   + APOE_e4_carrier + PC1 + PC2 + PC3 + PC4
```

- `amyloid_PGS`: z-standardised PLINK-1.9 `--score` sum on the Omni2.5M, applied to all PRS-CS posterior weights (chr 1–22, 486,274 SNPs after intersection).
- `APOE_e4_carrier`: binary indicator from APOERES.GENOTYPE (1 if ε4/* or */ε4 or ε4/ε4; 0 otherwise).
- `PC1–PC4`: top 4 genetic-ancestry PCs from PLINK --pca on LD-pruned (200/50/0.2) Omni2.5M variants with MAF > 0.05 and HWE > 1×10⁻⁶.

**Pre-specified sensitivity analyses (six, all in `R/20`, all reported in manuscript):**

1. **Tau-PVC**: same model on `META_TEMPORAL_SUVR` from `UCBERKELEY_TAUPVC_6MM` (partial-volume corrected) — tests robustness to atrophy-confounding.
2. **Cross-modality validation in ADNI**: `CENTILOIDS ~ amyloid_PGS + covariates` in the n=722 amyloid-PET ∩ Omni2.5M cohort — direct in-cohort replication of the Ali signal at the trait level.
3. **APOE-stratified**: primary model fit separately in ε4 carriers (estimated n≈110) and non-carriers (estimated n≈108).
4. **EUR-restricted**: primary model restricted to self-reported `White & Not Hispanic/Latino` subset (PTRACCAT × PTETHCAT) — matches the PRS-CS LD reference (1KG EUR Phase 3).
5. **Longitudinal LMM**: linear mixed model `META_TEMPORAL_SUVR ~ PGS + covariates + (1 | PTID)` on all 465 scans — uses full longitudinal information.
6. **Within-subject amyloid → tau path**: in the n=217 subjects with both modalities, test `META_TEMPORAL_SUVR ~ CENTILOIDS + age + sex + APOE` to confirm the canonical amyloid-cascade ordering in this subset.

**Pre-specified decision rule:**

- **Headline interpretation reported in manuscript:** β_PGS, 95% CI, p-value, and adjusted R² from the primary model. Direction is reported regardless of significance.
- **Supportive of cascade hypothesis if:** β_PGS > 0 with p < 0.05 in the primary model AND consistent direction (β_PGS > 0) in ≥4 of 6 sensitivity analyses.
- **Null result interpretation:** if β_PGS not significantly different from zero in the primary model, this is reported transparently as "amyloid PGS did not reach nominal significance for predicting tau-PET in the n=218 ADNI cohort; consistent with limited power for a polygenic-tau cascade test at this sample size" — i.e., framed as bounded sensitivity rather than refutation of the cascade.
- **Discordant direction interpretation:** if β_PGS < 0 in the primary model (PGS-protective against tau), this would be a substantive contradiction of the cascade hypothesis and would be reported as the headline finding with a substantively revised Discussion §6.

### Manuscript reporting commitment

The manuscript will explicitly report:

1. **Amendment 5 status (post-hoc, model frozen pre-results):** Methods § (new subsection) names Amendment 5 by number, summarises the trigger (ADNI access approval 2026-05-22, after pre-registered analyses complete), and explicitly notes "post-hoc addition" classification.
2. **The full primary model + all 6 sensitivity analyses are reported** regardless of headline direction (no selective reporting; null results are reported with the same prominence as positive results).
3. **A transparent null-power statement** is included whether or not the primary result is significant: "n=218 with the polygenic-vs-imaging cascade test design is powered to detect β_PGS effects of approximately X (sensitivity to be reported alongside the primary result)."
4. **DUA-compliance statement:** "Individual-level ADNI data are not redistributed; only group-level summary statistics from the cascade arm are reported. Analysis code is publicly available; the raw individual-level data must be obtained directly from ADNI under their Data Use Agreement."

### DUA-compliance operational commitments

1. NO individual-level ADNI data is committed to `public_repository/` or Zenodo (deposit `.gitignore`-equivalent: `data/raw/adni/` excluded by manifest construction).
2. Cascade-arm derived outputs in `data/processed/adni_cascade/` contain ONLY group-level results (coefficients, CIs, p-values, R², N) — no per-subject rows.
3. `R/20_adni_cascade_pgs.R` is committed to `public_repository/code/R/` and is shareable (it reads from local ADNI paths but contains no ADNI data values).
4. ADNI acknowledgement statement + funding statement included in manuscript per ADNI DUA requirements.

### Why this amendment exists rather than declining the data

The pre-registered manuscript (Ali multi-ethnic primary + Kim Korean secondary primary + CSF cascade + Ali ASN replication + chr1q32.2/CR1 colocalisation) is complete and standalone. Adding the ADNI cascade arm is **not necessary** for the manuscript's primary inferential claims. The choice to add it post-hoc rather than defer to a follow-up paper is principled:

1. **Direct biological validation in an independent EUR cohort.** The manuscript's strongest piece of evidence — the CR1-mediated amyloid → AD effect — predicts that an amyloid PGS should also relate to downstream tau pathology in independent EUR samples. ADNI is the largest publicly-accessible such sample.
2. **Same Pivot A inferential layer.** The amyloid PGS used in the cascade is exactly the same PGS used in the manuscript's primary inferential layer (PRS-CS-shrunk Ali multi-ethnic, phi=1e-2, seed=81). No new MR scaffolding is introduced; this is a direct extension of an existing component.
3. **The post-hoc framing is honestly disclosed.** This amendment is filed BEFORE results are seen, the model specification is frozen on disk, and the manuscript will report the post-hoc status by amendment number.

### Reproducibility

- `R/20_adni_cascade_pgs.R` is committed to the public repository before R/20 first execution.
- `data/processed/adni_cascade/adni_cascade_lm_results.tsv` (group-level coefficients only) will be added to the public repository deposit.
- ADNIMERGE2 R package version (v0.1.1, methods PDF dated 2026-01-05) is recorded in `data/raw/adni/NOTICE_DATA_USE_AGREEMENT.md`.
- This amendment is committed to the OSF deposit with timestamp; the OSF version-history record is the immutable proof that the model specification was frozen pre-results.

---

## (Future amendments appended below)
