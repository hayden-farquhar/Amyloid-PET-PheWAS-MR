# OSF Pre-Registration: Amyloid-PET Polygenic Burden as Exposure in a FinnGen R12 Phenome-Wide Mendelian Randomization

**Title:** Phenome-wide Mendelian randomization of genetically-proxied cortical amyloid-β deposition against 1,574 FinnGen R12 disease endpoints, with APOE stratification, ancestry-stratified instrument analyses, sex-stratified sub-analyses, colocalisation at lead loci, reverse MR for AD-related phecodes, mediation analysis through cognitive function, and equivalence testing for safety-relevant outcomes.

**Investigator:** Hayden Farquhar MBBS MPHTM
**ORCID:** 0009-0002-6226-440X
**Affiliation:** Independent researcher; library access via University of Sydney
**Email:** hayden.farquhar@icloud.com

**Date of pre-registration:** 2026-05-19
**Registration template:** OSF Prereg (Open-Ended)
**Status:** FROZEN upon OSF deposition — deviations logged in the amendment log

**Reporting standards:** STROBE-MR (Skrivankova et al. 2021); PRISMA-MR where applicable; equator-network-compliant.

**Document structure.** Sections are organised into eleven thematic parts:

- Part A — Context (§§1–4)
- Part B — Materials (§§5–8)
- Part C — Primary analysis (§§9–10)
- Part D — Decomposition analyses (§§11–18)
- Part E — Quality controls (§§19–20)
- Part F — Safety and clinical translation (§§21–23)
- Part G — Inference framework (§§24–28)
- Part H — Reporting (§§29–30)
- Part I — Limitations (§§31–34)
- Part J — Process and governance (§§35–39)
- Part K — Manifest and sign-off (§§40–41)

---

**Part A — Context**

---

## 1. Study Information

### 1.1 Funding

None. Sole-investigator project.

### 1.2 Conflicts of interest

None declared.

### 1.3 Ethics

HREC-exempt — all data sources are public, aggregate-level summary statistics. No individual-level data are accessed at any stage. See the ethics-determination companion document.

### 1.4 Code licence

MIT. Code deposited at GitHub at manuscript submission.

### 1.5 Data sharing

All processed results parquets, harmonised instrument tables, MR estimates, colocalisation outputs, and analysis code published to OSF and GitHub at submission. Raw GWAS summary statistics remain at their original public sources; SHA256 hashes and access dates of all input files are logged in the data-snapshot companion document.

---

## 2. Research Questions

### 2.1 Primary research question

Using cortical amyloid-PET deposition as the Mendelian-randomization exposure (Ali et al. 2023; *Acta Neuropath Comm*; DOI 10.1186/s40478-023-01563-4; N=13,409 multi-ethnic meta-analysis), which of the 1,574 FinnGen R12 phecode outcomes (case-count ≥ 500; pre-specified at the companion phecode-filter document) show evidence of a causal effect of genetically-predicted amyloid burden, robust to APOE-stratified sensitivity analyses, with concordance across pleiotropy-robust MR methods?

### 2.2 Secondary research questions

1. **APOE pleiotropy decomposition.** After excluding the APOE region (GRCh38 chr19:44,906,000–45,916,000), how many phecodes show evidence of a causal effect of non-APOE-driven amyloid burden? Multivariable MR (MVMR) with APOE ε4 dosage as a second exposure further decomposes APOE-region effects into ε4-mediated versus non-ε4-region pleiotropy.

2. **Cross-ancestry transportability.** Using the Ali 2023 ancestry-stratified sub-cohorts (Non-Hispanic White N=11,556; East Asian N=1,494; African N=359), do MR estimates show ancestry-specific heterogeneity (Cochran's Q across ancestries), and does a formal trans-ancestry meta-regression (MR-MEGA) identify ancestry-axis–correlated effects? Biobank Japan and Pan-UKB AFR phecode sumstats provide outcome-side cross-ancestry replication for FDR-significant hits.

3. **Sex-stratified causal effects.** Using the Ali 2023 sex-stratified sub-cohorts (female N=5,195; male N=4,625), do MR estimates show sex-specific causal effects (Cochran's Q across sexes) consistent with the documented sex differences in AD progression and amyloid-PET-positivity-to-dementia conversion?

4. **Directionality and AD-mediation.** For phecodes within the dementia/AD chapter (ICD-10 F00, F03, G30, G31), does Steiger directionality and formal reverse MR (using Bellenguez 2022 AD GWAS as a reverse exposure) support amyloid → AD rather than reverse causation? For non-AD phecodes with evidence of a causal effect, does multivariable MR with AD (or general cognitive function, Davies et al. 2018) as a mediator attenuate the direct amyloid → outcome effect?

5. **Amyloid-direct vs AD-mediated decomposition.** For every phecode, is the amyloid-PET causal effect distinct from the AD causal effect when both are estimated jointly via multivariable MR? The amyloid-on-outcome coefficient from an MVMR model that includes AD-risk as a second exposure (β_amyloid|AD) is the principled "amyloid-direct" estimator; the univariate difference β_residual = β_amyloid − β_AD is reported as a descriptive sanity check that ignores estimator covariance. This addresses the question "does amyloid burden have effects beyond what is captured by clinical AD diagnosis?"

6. **Safety equivalence.** For phecodes that map to known anti-amyloid therapeutic adverse events (intracerebral haemorrhage, cerebral oedema, anaphylaxis), is the MR effect statistically equivalent to no effect at a pre-specified safety margin (TOST framework)?

7. **APOE ε4 dose-response.** Beyond binary APOE-included vs APOE-excluded panels, does the amyloid → outcome MR effect scale with APOE ε4 dosage (0 / 1 / 2 ε4 alleles)? This is genuinely novel: no prior amyloid-PET MR has stratified by ε4 dosage tier, even though APOE ε4 is a well-established dosage-response modifier of amyloid accumulation kinetics.

8. **Pre-symptomatic vs symptomatic amyloid effects.** Is the amyloid → outcome effect detectable in the pre-symptomatic (cognitively-normal amyloid-positive) state, separately from the clinical-AD state? This directly addresses the central anti-amyloid-therapeutic question: should pre-symptomatic amyloid-positive individuals be treated, or should treatment wait for clinical AD onset?

9. **Tau-PET cascade decomposition.** For amyloid → outcome effects, what fraction is mediated through downstream tau-PET pathology (the canonical AD pathological cascade) versus what fraction is amyloid-direct independent of tau? A pre-registered two-step MR (amyloid → tau-PET → outcome) extends the GWAS-by-subtraction framework to the established pathological cascade.

10. **Non-linear amyloid dose-response.** Is the amyloid → outcome MR effect linear across the amyloid burden range, or does the effect have a threshold or curvilinear shape? This directly addresses the anti-amyloid-therapy linearity-of-effect assumption and informs whether partial amyloid clearance produces partial benefit (linear) versus all-or-nothing benefit (threshold).

---

## 3. Hypotheses

All hypotheses are directional and pre-specified. Significance is assessed at family-wise error rate α = 0.05 with Benjamini-Hochberg FDR control within hypothesis families.

| # | Hypothesis | Test | Decision rule |
|---|---|---|---|
| **H1** | At least 5 FinnGen R12 phecodes show a significant MR causal effect of genetically-proxied amyloid deposition at FDR q < 0.05 when APOE-inclusive instruments are used | IVW random-effects across the 1,574-phecode universe, BH-FDR within phecode chapter | H1 supported if ≥ 5 phecodes pass FDR q < 0.05 AND ≥ 3 of those show concordant direction across IVW + weighted median + weighted mode + MR-PRESSO (outlier-corrected) AND effect-size tier ≥ "weak support" (§26.3) |
| **H2** | At least 2 phecodes show a significant MR effect after APOE-region instruments are excluded, indicating non-APOE-driven amyloid causal effects | IVW random-effects on APOE-excluded instrument panel | H2 supported if ≥ 2 phecodes pass FDR q < 0.05 with effect-size tier ≥ "weak support" AND APOE-excluded instrument mean F-statistic ≥ 10 AND n_IV ≥ 5 (so the test is interpretable rather than uninformative) |
| **H3** | For clinical-AD phecodes (F00 Vascular dementia, F03 Unspecified dementia, G30 Alzheimer's disease, G31 Other degenerative diseases of the nervous system), Steiger directionality and formal reverse MR support amyloid → AD rather than reverse | TwoSampleMR Steiger filter per SNP; reverse MR using Bellenguez 2022 AD GWAS (Nat Genet; DOI 10.1038/s41588-022-01024-z; N=487,511 EUR) as a reverse exposure against amyloid-PET as outcome | H3 supported if Steiger directionality favours amyloid → AD in ≥ 80% of instruments AND reverse-MR IVW effect estimate for AD → amyloid is non-significant or smaller in magnitude than forward MR |
| **H4** | At least 1 non-APOE significant hit replicates in a cross-ancestry outcome panel (Biobank Japan matched phenotypes or Pan-UKB AFR) at nominal p < 0.05 with concordant direction | IVW on outcome-side cross-ancestry sumstats for the same instrument set | H4 supported if ≥ 1 non-APOE hit replicates at nominal p < 0.05 with concordant direction in at least one of {BBJ matched, Pan-UKB AFR matched} |
| **H5** | Ancestry-stratified instrument analyses (Ali 2023 NHW, ASN, AFR arms) demonstrate cross-ancestry transportability of the primary findings. The primary verdict is operationalised as **directional consistency** (sign of the IVW estimate concordant between the EUR-heavy primary and each ancestry arm where n_IV ≥ 2); the secondary verdict is formal trans-ancestry heterogeneity testing where powered | Sign-concordance check per FDR-significant non-APOE hit across NHW / ASN / AFR arms; secondary: MR-MEGA trans-ancestry meta-regression where powered | H5 (primary, directional) supported if direction agrees in ≥ 2 of 3 ancestry arms for ≥ 1 non-APOE hit. H5 (secondary, heterogeneity) supported if MR-MEGA additionally identifies ≥ 1 phecode with ancestry-axis–correlated effect (p_anc < 0.05 after BH-FDR). The primary directional verdict is more achievable than the heterogeneity test given the small ASN (N=1,494) and AFR (N=359) sample sizes |
| **H6** | Sex-stratified instrument analyses (Ali 2023 female vs male arms) show that at least one phecode demonstrates significantly heterogeneous effects across sexes (Cochran's Q across-sex p < 0.05 after FDR), with greater effect magnitude in females consistent with the documented sex differences in amyloid-positivity-to-AD conversion | Sex-stratified IVW; Cochran's Q test for between-sex heterogeneity per phecode | H6 supported if ≥ 1 phecode shows significant between-sex heterogeneity AND direction of female-male difference is consistent with female-amplified effect for AD-spectrum and psychiatric phecodes |

### 3.1 Hypothesis families for multiple-testing control

- **Family A (primary):** H1 (APOE-inclusive sweep across 1,574 phecodes × IVW). BH-FDR within each of 19 ICD-10–derived chapters; familywise α = 0.05.
- **Family B (primary):** H2 (APOE-excluded sweep across 1,574 phecodes × IVW). BH-FDR within chapter; familywise α = 0.05.
- **Family C (secondary):** H3, H4, H5, H6 each tested on their own hit-restricted sub-families. Family-C Holm step-down across the 4 secondary hypotheses at α = 0.05.

### 3.2 Null response plan and pre-specified pivots

If H1 returns < 5 phecodes at FDR q < 0.05, OR if APOE-excluded instruments fail the gate (n_IV < 5 OR mean F < 15), the project pivots to one of three pre-specified secondary deliverables. The pivots are not framed by submission venue; they are framed by what scientific deliverable replaces the primary finding:

1. **Pivot A — Weak-instrument atlas + power-calculation toolkit.** Trigger: APOE-excluded n_IV < 5 OR mean F < 15. Deliverable: power atlas + analytic minimum-detectable-OR per (phecode × instrument set) + recommendations for amyloid-PET GWAS sample-size targets.

2. **Pivot B — APOE-pleiotropy atlas against the FinnGen R12 phenome.** Trigger: H1 ≥ 5 phecodes pass but H2 returns 0–1 phecodes. Deliverable: per-phecode APOE ε4 vs other-APOE-region vs APOE-excluded decomposition; multivariable MR with ε4 dosage; mapping of which non-AD phecodes are APOE-pleiotropy-explained.

3. **Pivot C — Sex-stratified amyloid-PET MR atlas.** Trigger: H1 + H6 both supported, with ≥ 3 phecodes showing sex-heterogeneous effects (Q p < 0.05 after FDR). Deliverable: sex-stratified atlas; sex-difference mechanistic discussion; mediation through sex-specific intermediate phenotypes.

Operational triggers and full deliverable definitions are in §27.

---

## 4. Prior Art and Positioning

### 4.1 Mechanistic ancestor literature

| Study | Approach | What it does | What it does not |
|---|---|---|---|
| Ali et al. 2023 *Acta Neuropath Comm* (DOI 10.1186/s40478-023-01563-4) | Amyloid-PET GWAS (N=13,409, multi-ethnic + ancestry-stratified + sex-stratified) | Identifies novel non-APOE amyloid-PET loci | No downstream MR application |
| Kim et al. 2025 *Nat Commun* (DOI 10.1038/s41467-025-57751-4) | Cross-ancestry amyloid-PET GWAS (N=15,701; EUR + Korean EAS) | Identifies SORL1 as a novel cross-ancestry locus | No PheWAS-MR application; sumstats deposition pending |
| Bellenguez et al. 2022 *Nat Genet* (DOI 10.1038/s41588-022-01024-z) | AD GWAS (N=487,511 EUR) | Identifies 75+ AD risk loci | Not amyloid-PET-specific |
| Heo et al. 2025 *Dement Neurocog Disord* | MR: Aβ-PET → tau-PET | Single-outcome MR for tau accumulation | No PheWAS sweep |
| Ji et al. 2026 *Brain Imaging Behav* | MR: AD → brain imaging traits | Reverse direction | Amyloid-PET not the exposure |
| Lindbohm et al. 2022 *Nat Aging* | MR: immune markers → dementia | Immune-marker exposures | Amyloid-PET not the exposure |
| Smith et al. 2021 BIG40 | UKB brain imaging MR-adjacent | UKB brain phenotypes | Not PET; not amyloid-specific |
| Jansen et al. 2022 EADB | CSF Aβ42 GWAS | CSF biomarker exposure | Not PET |
| Brumpton et al. 2020 *Nat Comm*; Howe et al. 2022 *Nat Genet* | Within-family MR critiques | Demonstrate that conventional MR is vulnerable to assortative-mating, population stratification, and dynastic effects | No within-family amyloid-PET GWAS exists; this study uses conventional population-level MR and acknowledges this as a limitation (§34) |

### 4.2 Gap statement

To our knowledge, no published phenome-wide MR uses amyloid-PET cortical deposition as the exposure across the full FinnGen R12 phecode universe. No published amyloid-PET MR uses APOE-stratified instrument panels as a pre-registered decomposition design. No published amyloid-PET MR uses ancestry-stratified or sex-stratified Ali 2023 sub-cohorts as instrument sources. No published amyloid-PET MR uses GWAS-by-subtraction to decompose amyloid-direct from AD-mediated effects. This study addresses all four gaps in a single pre-registered design.

### 4.3 What this study adds

1. First PheWAS-MR with amyloid-PET as exposure across ~1,500 phecodes.
2. Pre-registered APOE-stratified decomposition (inclusive + excluded + ε4 MVMR + ε4 dosage-response).
3. Pre-registered ancestry-stratified instrument analysis using Ali 2023's NHW/ASN/AFR arms, formalised via MR-MEGA.
4. Pre-registered sex-stratified sub-analysis using Ali 2023's F/M arms — novel for amyloid-PET MR.
5. Pre-registered reverse MR for AD-spectrum phecodes using Bellenguez 2022 AD as the reverse exposure.
6. Pre-registered mediation analysis through AD / cognitive function as a candidate mediator.
7. Pre-registered amyloid-direct vs AD-mediated decomposition via MVMR with AD as a co-exposure (§14).
8. Pre-registered equivalence testing (TOST framework, §22) for phecodes that map to anti-amyloid therapeutic adverse events.
9. Pre-registered LDSC genetic-correlation cross-tabulation (§18) as complementary causal-inference signal.
10. Pre-registered polygenic-score-MR sensitivity as an alternative instrument strategy for the weak-instrument regime (§6.10).
11. Pre-registered colocalisation at top non-APOE loci using GTEx v8 brain tissues + eQTLGen blood.
12. Pre-registered power atlas as a standalone secondary deliverable if the primary sweep returns < 5 hits.
13. **Pre-registered APOE ε4 dose-response MR** (§11.4) — first amyloid-PET MR to stratify by ε4 dosage tier (0/1/2 alleles).
14. **Pre-registered pre-symptomatic vs symptomatic amyloid decomposition** (§14.4) — directly addresses anti-amyloid therapeutic timing decisions.
15. **Pre-registered tau-PET cascade two-step MR** (§13.5) — decomposes amyloid → tau-mediated vs amyloid-tau-independent effects per phecode.
16. **Pre-registered non-linear MR with fractional-polynomial dose-response curves** (§10.5) — tests amyloid → outcome linearity, directly relevant to partial-therapeutic-clearance reasoning.
17. **Pre-registered multi-ancestry SuSiE fine-mapping** (§17.7) — for FDR-significant cross-ancestry hits, push from causal-effect estimation to causal-variant identification across ancestries.
18. **Pre-registered post-marketing surveillance triangulation** (§23.4) — forward-looking commitment to cross-reference significant MR hits against FAERS/DAEN signals for lecanemab/donanemab at manuscript submission time.
19. **Pre-registered Amyloid-PET Causal Atlas** (§38.5) — full per-phecode results published as a downloadable parquet + static HTML/PDF report with DOI, available as a community resource regardless of the primary-paper outcome.

### 4.4 Clinical motivation

Anti-amyloid monoclonal antibody therapeutics (lecanemab, donanemab) entered routine clinical use in 2024–2025. Off-target causal effects of amyloid burden on non-AD phenotypes — both protective and harmful — directly inform the therapeutic label, screening indications, and the clinical decision of whether to extend anti-amyloid therapy to pre-symptomatic amyloid-PET-positive individuals. Existing MR on dementia uses clinical AD diagnosis or CSF biomarkers as exposures, missing the pre-clinical amyloid-positive window that PET captures and that anti-amyloid therapeutics target.

---

**Part B — Materials**

---

## 5. Exposure Definition

### 5.1 Primary exposure

**Phenotype:** Cortical amyloid-β deposition, measured by PET imaging.

**Exposure GWAS:** Ali et al. 2023 *Acta Neuropath Comm* multi-ethnic meta-analysis (N=13,409 across 14 cohorts). Genome build: GRCh38.

**Effect scaling:** All MR estimates reported per **1 standard deviation increase in genetically-proxied cortical amyloid-PET signal**. Effect-size translation to Centiloid units and to lecanemab-equivalent clearance magnitudes is reported per §23.

**Rationale.** Ali 2023 is the largest available amyloid-PET GWAS as of 2026-05-19, with multi-ethnic ancestry stratification (NHW, ASN, AFR sub-cohorts) and sex stratification (F, M sub-cohorts) available as separate sub-files. This stratification supports the H5 (cross-ancestry) and H6 (sex-stratified) hypotheses pre-registered in §3.

**Snapshot provenance.** SHA256 hashes for all 7 Ali 2023 sub-cohort files are logged in the data-snapshot companion document at the time of OSF deposition.

### 5.2 Contingent additional primary exposure (Kim 2025)

If the Kim et al. 2025 *Nat Commun* cross-ancestry amyloid-PET GWAS summary statistics (N=15,701 = 11,816 EUR + 3,885 EAS Korean; GWAS Catalog accession GCST90483382) become available from the corresponding author at SKKU (request sent 2026-05-19) by Phase 3 closure (Week 10), Kim 2025 is added as a co-primary exposure under the protocol in the exposure-choice companion document.

If Kim 2025 sumstats are NOT obtained by Phase 3 closure, the Ali 2023 ASN arm (N=1,494) serves as the EAS-ancestry instrument source for H5, with explicit acknowledgement in the manuscript that Kim 2025 would have provided substantially greater EAS power.

### 5.3 Replication exposure (CSF Aβ42)

For Phase 5 sensitivity, top non-APOE hits are re-tested using the Jansen et al. 2022 EADB CSF Aβ42 GWAS as an orthogonal exposure capturing soluble CSF amyloid rather than aggregated cortical amyloid. Concordant direction across CSF Aβ42 and amyloid-PET exposures strengthens the causal interpretation; discordance flags the possibility that PET-specific aggregated-amyloid effects differ from CSF-Aβ42 effects.

### 5.4 Secondary sensitivity exposure (Nho 2024 tau-PET)

For phecodes with significant amyloid-PET MR effects in the AD-spectrum or psychiatric chapters, Nho et al. 2024 ADNI tau-PET GWAS is used as a secondary exposure to test whether tau-PET shows concordant or discordant effects. Tau-PET is downstream of amyloid in the AD pathological pathway; concordance suggests a shared mechanism, discordance suggests amyloid-specific causal architecture.

---

## 6. Instrument Construction

### 6.1 Three-tier instrument strategy

Given the relatively modest exposure GWAS sample size (Ali 2023 N=13,409) and the expected sparsity of non-APOE instruments at genome-wide significance, three pre-registered instrument tiers are used:

| Tier | Threshold | Primary estimator | Role |
|---|---|---|---|
| **1 — Strict** | p < 5 × 10⁻⁸; F ≥ 10; r² < 0.001 | IVW random-effects | Primary inference |
| **2 — Relaxed** | p < 1 × 10⁻⁶; F ≥ 10; r² < 0.001 | MR-RAPS (weak-instrument robust) | Sensitivity widening if Tier 1 is sparse outside APOE; reported separately |
| **3 — PGS-MR (fallback)** | PRS-CS polygenic score as single composite instrument | Two-sample MR per Burgess et al. 2017 framework | Activated as primary if LDSC h² gate triggers Pivot A (§27) |

The tier-used-as-primary depends on the APOE-excluded gate (§6.4) and the LDSC heritability QC (§8.5). Tier 1 is the default; Tier 2 is reported as a transparency mechanism when Tier 1 has n_IV < 10 in the APOE-excluded panel; Tier 3 becomes primary if Pivot A is activated.

Pre-registering all three tiers prevents post-hoc threshold switching: the operative tier is determined by the gate outcomes, not by which tier "produces results".

### 6.2 LD clumping

- LD reference panel: 1000 Genomes Phase 3 EUR for primary multi-ethnic and NHW arms; 1000G EAS for ASN arm; 1000G AFR for AFR arm.
- Clumping r² < 0.001.
- Clumping window 10 Mb (consistent with conventional cis-MR practice).
- Implementation: `ieugwasr::ld_clump()` or local PLINK 1.9 binary against the reference `.bed/.bim/.fam` files; choice fixed before instrument extraction.

### 6.3 Instrument strength (F-statistic)

- Per-SNP F-statistic computed as F = β² / SE².
- Inclusion threshold F ≥ 10 (the conventional weak-instrument cutoff, Stock & Yogo 2005).
- Mean F-statistic reported per instrument panel; instrument-panel-level F < 15 flagged in Limitations as a weak-instrument risk.
- Conditional F-statistics (Sanderson et al. 2019) computed for multivariable MR settings (§11.2).

### 6.4 Two pre-specified instrument panels

- **APOE-inclusive panel:** all SNPs passing p < 5 × 10⁻⁸, F ≥ 10, post-clumping.
- **APOE-excluded panel:** as above, with all SNPs in chr19:44,906,000–45,916,000 (GRCh38) removed.

### 6.5 Pre-registered ancestry-stratified panels

In addition to the multi-ethnic primary, instruments are extracted independently from:

- Ali 2023 NHW arm (N=11,556) → LD clump against 1000G EUR.
- Ali 2023 ASN arm (N=1,494) → LD clump against 1000G EAS.
- Ali 2023 AFR arm (N=359) → LD clump against 1000G AFR.

Expected instrument count is low for ASN and AFR given sample size; these panels are pre-registered as ancestry-specific replication anchors, not as primary discovery panels.

### 6.6 Pre-registered sex-stratified panels

- Ali 2023 female arm (N=5,195) → multi-ethnic LD clump.
- Ali 2023 male arm (N=4,625) → multi-ethnic LD clump.

Sex-stratified IVW reported per phecode alongside the pooled estimate, with Cochran's Q across-sex heterogeneity test.

### 6.7 SNP type restrictions

- Bi-allelic SNVs only. Indels excluded.
- Palindromic SNPs (A/T or G/C) with MAF 0.42–0.58 excluded as ambiguous; otherwise aligned by effect-allele frequency.
- SNPs with MAF < 0.01 excluded.
- SNPs failing genome-build harmonisation (no GRCh38 position) excluded.

### 6.8 Proxy SNPs

If an instrument SNP is absent from an outcome GWAS, an LD proxy is sought at r² ≥ 0.8 from an ancestry- and population-matched LD reference panel:

- **For FinnGen R12 outcomes:** SISu Finnish-population LD reference panel is used preferentially (FinnGen is a founder population with distinct LD architecture due to historical bottlenecks; 1000G EUR LD does not adequately represent Finnish haplotype structure for some loci). 1000G EUR is the fallback if a variant is not in the SISu panel.
- **For BBJ outcomes:** 1000G EAS LD panel.
- **For Pan-UKB AFR outcomes:** 1000G AFR LD panel.
- **For Bellenguez 2022 AD (reverse-MR exposure):** 1000G EUR LD panel.

If no proxy at r² ≥ 0.8 exists in the appropriate panel, the SNP is dropped for that outcome and the per-(instrument × outcome) panel is documented. Proxy-source panel is recorded per proxy in `data/processed/proxy_sources.parquet`.

### 6.9 Minimum instrument count rules

| n instruments | Analysis permitted |
|---|---|
| 0 | Outcome excluded for that panel |
| 1 | Wald ratio only; flagged hypothesis-generating |
| 2 | IVW only (no Egger or heterogeneity tests) |
| 3 | IVW + Egger + weighted median |
| 4 | All primary + sensitivity methods including MR-PRESSO |
| ≥ 5 | All methods including CAUSE and MR-RAPS |

### 6.10 Polygenic-score sensitivity instruments (pre-registered)

For the weak-instrument regime (APOE-excluded panel where n_IV < 10 OR mean F < 20), pre-register a polygenic-score (PGS) MR sensitivity analysis:

1. Compute an amyloid-PET PGS using PRS-CS (Ge et al. 2019 *Nat Comm*; DOI 10.1038/s41467-019-09718-5) on the Ali 2023 multi-ethnic GWAS with 1000G EUR LD panel.
2. Use the PGS as a single composite instrument in two-sample MR per the Burgess et al. 2017 framework for genetic risk score MR.
3. Pre-register the inclusion threshold for PGS-MR primary: the APOE-excluded panel's polygenic score must exceed an LDSC-derived heritability floor (h² > 0.05) on the Ali 2023 multi-ethnic data.

The PGS-MR analysis sacrifices per-locus interpretability for greater statistical power in the weak-instrument regime. Reported as a sensitivity for the APOE-excluded primary analysis. If the clumped-IV primary fails (Pivot A activated) but PGS-MR detects effects, the manuscript reports both — with explicit framing that PGS-MR does not localise the causal signal to specific loci.

---

## 7. Outcome Definition

### 7.1 Primary outcome universe

**FinnGen R12 phecode summary statistics**, public release (Google Cloud Storage bucket `gs://finngen-public-data-r12/summary_stats/release/`; manifest URL: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`).

Total endpoints in R12: 2,469.

**Pre-registered filter** (companion CSV, deposited with this registration):
- Case count ≥ 500.
- Yields **1,574 phecodes** retained across 19 ICD-10–aligned chapters and 22 FinnGen-specific category groupings.

Parent/child phecode duplication is NOT auto-collapsed at this stage. Within-chapter correlation structure is handled at the multiple-testing layer via the Benjamini-Yekutieli (2008) hierarchical FDR procedure as a sensitivity (§24.2).

### 7.2 Outcome harmonisation

For each (instrument-panel × phecode) combination:
- TwoSampleMR `harmonise_data(action = 2)` for palindromic handling.
- Outcome effect direction harmonised to match the exposure effect allele.
- Outcome SNPs failing build harmonisation (GRCh38 alignment) excluded.

### 7.3 Replication outcomes

For each phecode showing a significant MR effect (FDR q < 0.05) in the FinnGen primary analysis:

- **Biobank Japan** matched phenotype sumstats (where available; phenotype-matching best-effort).
- **Pan-UKB AFR** matched phenotype sumstats (where available; small case counts expected).
- **Bellenguez 2022 non-FinnGen sub-meta-analysis** where AD-spectrum phecodes are involved (to avoid Finnish-cohort dependence in the AD replication).

Replication is reported as nominal p < 0.05 with concordant direction. No FDR is applied at the replication layer because the testing is hit-restricted.

### 7.4 Binary-outcome interpretation caveat

FinnGen R12 phecode summary statistics are produced from logistic-regression GWAS on a binary case/control definition (any time during follow-up). The MR estimates correspond to a binary log-odds-ratio interpretation. Time-to-event MR (Cox proportional-hazards MR; MR-CIST methods) would require individual-level data and is out of scope. The binary-OR interpretation represents the genetic association with lifetime risk of the outcome rather than instantaneous hazard, and is reported as such in Limitations.

---

## 8. Harmonisation Protocol

### 8.1 Genome build

All sumstats harmonised to **GRCh38**. Ali 2023 is native GRCh38. FinnGen R12 is native GRCh38. CSF Aβ42 (Jansen 2022) is native GRCh38. Nho 2024 tau-PET is GRCh37 → lifted to GRCh38 using UCSC liftOver chain file. BBJ is GRCh37 → lifted. Pan-UKB AFR is GRCh37 → lifted.

SNPs failing liftover are dropped and the failure rate reported.

### 8.2 rsID re-annotation

Ali 2023 uses `chr:pos:ref:alt` as the SNP ID. For `ieugwasr::ld_clump()` and TwoSampleMR compatibility, instruments are re-annotated to rsIDs via the `MungeSumstats` package against dbSNP build 155 (or the most recent at analysis time, fixed before instrument extraction begins). SNPs without a canonical rsID at the matched position+alleles are dropped at the re-annotation step and the dropped count reported.

### 8.3 Effect-allele coding

All exposure and outcome sumstats harmonised so that the effect allele is the alternate (ALT) allele. Effect-allele frequency (EAF) cross-checked between exposure and outcome; SNPs with |EAF_exposure − EAF_outcome| > 0.15 flagged as potential strand mismatches and excluded as sensitivity.

### 8.4 Inflation factor and population-stratification check

For every outcome sumstats file ingested, the genomic-control inflation factor (λ_GC) is computed on the full GWAS. Outcomes with λ_GC > 1.10 are flagged in the manuscript as having potential residual population stratification; MR estimates against such outcomes are reported with explicit caveat. λ_GC > 1.20 triggers exclusion from primary analysis (sensitivity-only).

### 8.5 LD score regression heritability + intercept QC

LD score regression (LDSC; Bulik-Sullivan et al. 2015 *Nat Genet*; DOI 10.1038/ng.3211) is run on every outcome sumstats and on the Ali 2023 exposure sumstats as a QC step. Pre-registered thresholds:

| Metric | Threshold | Action |
|---|---|---|
| Exposure heritability h² (Ali 2023) | h² > 0.05 | Required for primary analysis; if h² ≤ 0.05, switch to PGS-MR sensitivity (§6.10) as primary |
| Outcome heritability h² (per phecode) | h² > 0.01 | Required for inclusion in MR-meaningful analysis; below threshold flagged as "MR uninformative" |
| LDSC intercept | < 1.10 | If ≥ 1.10, possible population stratification or sample-overlap; flagged in Limitations |
| Exposure × outcome rg | reported | LDSC genetic correlation as complementary causal-inference signal (§18) |

LDSC reference panel: 1000G EUR HapMap3 SNPs.

---

**Part C — Primary analysis**

---

## 9. Primary MR Analysis

### 9.1 Primary estimator

**Inverse-variance weighted (IVW) random-effects** per (instrument-panel × phecode) combination.

Implementation: `TwoSampleMR::mr_ivw_re()`.

### 9.2 Output per primary run

| Field | Description |
|---|---|
| `phecode` | FinnGen endpoint identifier |
| `instrument_set` | APOE-inclusive \| APOE-excluded |
| `method` | IVW (primary); other methods in sensitivity |
| `nsnp` | Number of harmonised instruments |
| `beta` | IVW point estimate (per 1-SD amyloid-PET) |
| `se` | Standard error |
| `pval` | Primary p-value |
| `pval_fdr_chapter` | BH-adjusted within phecode chapter |
| `pval_fdr_global` | BH-adjusted across full 1,574-phecode sweep |
| `pval_bonferroni_global` | Bonferroni (n = 1,574 × 2) sensitivity |
| `pval_by_chapter` | Benjamini-Yekutieli within-chapter (correlation-robust FDR) |
| `steiger_pval` | Steiger directionality test p-value |
| `q_pval` | Cochran's Q heterogeneity p-value |
| `i_squared` | I² heterogeneity statistic |
| `instrument_f_mean` | Mean F-statistic of instrument panel |
| `case_count` | FinnGen R12 case count for the phecode |
| `control_count` | FinnGen R12 control count for the phecode |
| `ldsc_h2_outcome` | Outcome heritability from LDSC |
| `ldsc_intercept` | Outcome LDSC intercept |
| `ldsc_rg_exposure_outcome` | Genetic correlation between exposure and outcome |

### 9.3 Output deposit

Saved to `data/processed/mr_primary_results.parquet`; deposited to OSF at submission.

---

## 10. Sensitivity Suite

### 10.1 Required sensitivity methods for any FDR-significant hit

For every (instrument-panel × phecode) hit at FDR q < 0.05, the following sensitivity panel is mandatory:

| Method | Role | Reference |
|---|---|---|
| **MR-Egger regression** | Directional pleiotropy intercept test | Bowden et al. 2015 |
| **Robust MR-Egger** | Egger variant robust to invalid instruments | Bowden et al. 2018 |
| **Weighted median** | Robust to ≤ 50% invalid instruments | Bowden et al. 2016 |
| **Weighted mode** | Robust to plurality-invalid instruments | Hartwig et al. 2017 |
| **MR-PRESSO** (1000 distributions; outlier + global test) | Horizontal pleiotropy + outlier detection | Verbanck et al. 2018 |
| **MVMR-PRESSO** | MR-PRESSO extension for multivariable MR (used in §11.2 APOE ε4 MVMR) | Sanderson et al. 2020 |
| **RadialMR** | Modified Cochran's Q outlier detection | Bowden et al. 2018 |
| **MR-RAPS** (Robust Adjusted Profile Score) | Robust to weak-instrument bias AND mild horizontal pleiotropy | Zhao et al. 2020 |
| **CAUSE** (Causal Analysis Using Summary Effects) | Robust to correlated horizontal pleiotropy (critical given APOE) | Morrison et al. 2020 |
| **MR-ConMix** (Contamination Mixture) | Robust to outlier instruments using mixture model | Burgess et al. 2020 |
| **MR-Lasso** | Penalised regression-style invalid-instrument detection | Slob & Burgess 2020 |
| **Leave-one-out IVW** | Single-SNP influence on the IVW estimate | Standard |
| **Steiger directionality** (per-SNP) | SNPs with R²_outcome > R²_exposure removed; IVW re-estimated | Hemani et al. 2017; Burgess & Thompson 2017 corrected R² scaling |

### 10.2 Concordance threshold

A primary IVW hit is declared **robust** if:
- IVW, weighted median, weighted mode, and MR-PRESSO outlier-corrected estimates all share the same direction.
- MR-Egger intercept p > 0.05 (no significant directional pleiotropy) OR the MR-Egger slope direction concords with IVW.
- CAUSE posterior probability of a causal effect (γ ≠ 0) > 0.5.
- MR-RAPS estimate direction concords with IVW.
- MR-Lasso instrument-set agrees with PRESSO outlier-detection set (no contradiction).

A hit failing 2+ of these is downgraded to "suggestive, requires replication" in the manuscript.

### 10.3 Steiger filter implementation

Per-SNP Steiger filtering with corrected R² scaling on both sides of the comparison:

- **Exposure side (continuous amyloid-PET):** Burgess & Thompson (2017) corrected R² scaling. The default `TwoSampleMR::directionality_test()` uses a binary-trait approximation that is biased for continuous exposures.
- **Outcome side (binary phecode):** Lee et al. (2012) *Genet Epidemiol* (DOI 10.1002/gepi.21614) prevalence-adjusted R² formula. Without this correction, low-prevalence phecodes (most of the 1,574 are at <1% case prevalence in FinnGen R12) systematically receive inflated outcome-R² estimates, biasing the Steiger test toward false-positive reverse-causality flags.

The combined approach — Burgess & Thompson 2017 for the continuous exposure plus Lee 2012 for the binary outcome — is the methodologically correct two-sided Steiger implementation for a continuous-exposure × binary-outcome MR.

Steiger-filter multiple-testing correction (Bonferroni across applied SNPs) reported as sensitivity.

### 10.4 Alternative LD-clumping sensitivity

Instruments re-extracted at r² < 0.01 and r² < 0.05; IVW re-estimated. Concordance with the primary r² < 0.001 panel reported.

### 10.5 Non-linear MR and dose-response curve construction

The conventional linear-MR assumption (constant β across the exposure range) is questionable for amyloid burden — biological intuition and therapeutic-trial observations suggest possible threshold effects (e.g., a Centiloid threshold above which AD risk accelerates) and non-linear pharmacodynamics under therapeutic clearance.

#### 10.5.1 Polygenic-score stratified non-linear MR

For FDR-significant hits, the amyloid-PET PGS (§6.10) is used to stratify the FinnGen R12 outcome population into amyloid-burden quintiles (Q1 lowest predicted amyloid, Q5 highest). Within each quintile, the genotype-exposure association is estimated from external cohorts (Ali 2023) and within-quintile MR estimates are constructed using the doubly-ranked stratification method (Tian et al. 2023 *Stat Med*) to avoid the assumption of homogeneous genetic effects across the exposure distribution.

Doubly-ranked stratification (Tian 2023) addresses the principal critique of earlier residual-stratification methods (Hamilton et al. 2024) — the previous residual-method NLMR results were shown to be biased under realistic genetic effect heterogeneity. The doubly-ranked method is the current methodologically defensible non-linear MR estimator for binary outcomes.

#### 10.5.2 Fractional-polynomial dose-response curves

For each FDR-significant non-AD-spectrum hit and for the AD positive-control phecode, construct the full dose-response curve using fractional-polynomial MR per Sulc et al. 2022 (*Genet Epidemiol*; DOI 10.1002/gepi.22507). Outputs:

- Best-fitting fractional-polynomial shape (linear, J-shape, U-shape, threshold).
- Likelihood-ratio test of non-linearity vs linear model.
- Per-Centiloid risk curve with 95% bootstrap CI ribbon.
- Identification of any inflection / threshold point (if non-linear).

Pre-registered interpretation:
- **Linear:** therapeutic amyloid reduction produces proportional risk reduction across the burden range.
- **Threshold (e.g., risk accelerates above ~50 Centiloid):** therapeutic amyloid reduction below the threshold may not produce risk reduction; clinical-trial inclusion criteria should target above-threshold populations.
- **U-shape or J-shape:** complex relationship requiring further investigation.

The dose-response curve is reported per phecode of clinical interest (§21) and per AD positive control as a primary deliverable. This is the first dose-response amyloid-PET MR at the phecode-wide scale.

#### 10.5.3 Pre-registered honest power statement

Non-linear MR has substantially higher data requirements than linear MR. With Ali 2023 N=13,409 and even the largest FinnGen phecodes, the dose-response curve estimates will be noisy. Curves are reported with full uncertainty quantification (bootstrap CIs) and the linearity-test verdict is reported as the primary inferential output, with the curve itself reported descriptively.

---

**Part D — Decomposition analyses** (per FDR-significant hit)

---

## 11. APOE Stratification and ε4 Dosage Decomposition

### 11.1 Three-tier APOE analysis

For every phecode passing FDR q < 0.05 in either H1 or H2:

| Tier | Instrument set | Interpretation |
|---|---|---|
| 1 — Inclusive | All SNPs (incl. APOE region) | Total amyloid-burden effect; can be APOE-confounded |
| 2 — Excluded | APOE region (chr19:44.9–45.9 Mb GRCh38) removed | Non-APOE-driven amyloid effect |
| 3 — APOE-ε4 only | Single instrument: rs429358 (APOE ε4 tag SNP) | Pure APOE ε4 dosage effect (Wald ratio) |

### 11.2 APOE ε4 dosage multivariable MR

For phecodes with significant Tier 1 effect, MVMR is performed with two exposures:
1. Genetically-proxied amyloid-PET burden (instrument panel: APOE-excluded SNPs from Ali 2023).
2. APOE ε4 dosage (instrument: rs429358 as a single SNP).

MVMR implementation: `MendelianRandomization::mr_mvivw()`. Conditional F-statistics reported per exposure (Sanderson et al. 2019). MVMR-PRESSO (Sanderson et al. 2020) applied as outlier-detection sensitivity within the MVMR framework.

Interpretation:
- **Both effects significant in MVMR:** independent contributions of non-APOE amyloid AND APOE ε4 dosage.
- **Only APOE ε4 significant after conditioning:** the apparent total effect is APOE-pleiotropy-driven.
- **Only non-APOE amyloid significant:** APOE-independent amyloid causal effect.

### 11.3 APOE pleiotropy classification

Each FDR-significant phecode classified into one of four categories:

| Category | Tier 1 | Tier 2 | Tier 3 | Interpretation |
|---|---|---|---|---|
| A | sig | sig | sig | Independent APOE ε4 + non-APOE amyloid effects |
| B | sig | sig | ns | Primarily non-APOE-amyloid–driven |
| C | sig | ns | sig | Primarily APOE-ε4–driven (APOE pleiotropy) |
| D | sig | ns | ns | APOE-region pleiotropy outside ε4 |

This classification frames the manuscript's per-hit interpretation and is reported in the Results.

### 11.4 APOE ε4 dose-response stratified MR

Beyond the binary inclusive/excluded panels, the amyloid → outcome MR effect is estimated stratified by APOE ε4 dosage:

- **ε4-0 stratum** (non-carriers): instruments derived from Ali 2023 with rs429358 weighted as 0 dose.
- **ε4-1 stratum** (heterozygous carriers): 1 dose.
- **ε4-2 stratum** (homozygous carriers): 2 dose.

#### 11.4.1 Implementation

Two complementary approaches:

1. **Outcome-stratified MR.** If FinnGen R12 provides APOE-ε4-stratified GWAS summary statistics for the outcomes (per-stratum case/control counts and per-stratum SNP-outcome betas), run separate per-stratum IVW for each phecode. If FinnGen does NOT release stratified outcomes, fall back to (2).
2. **Interaction-term MVMR.** Within an MVMR framework, include the (amyloid-PET × ε4-dosage) interaction term as a third exposure. The interaction-term coefficient is the dose-response modifier estimate.

The chosen route is fixed before analysis (route 1 if stratified outcomes are publicly available; route 2 otherwise). Pre-registered decision logged in the amendment log at instrument-extraction time.

#### 11.4.2 Decision rule

A phecode is declared **ε4-dose-modified** if:
- The IVW estimate in the ε4-2 stratum differs from the ε4-0 stratum by ≥ 1 effect-size tier (§26.3) AND
- Cochran's Q test across the three dosage tiers returns p < 0.05.

#### 11.4.3 Novelty statement

No prior amyloid-PET MR has stratified by ε4 dosage tier, despite the well-established dosage-response modifier role of APOE ε4 on amyloid accumulation kinetics (Reiman et al. 2009 *PNAS*; Jansen et al. 2015 *JAMA*). The dose-response analysis directly informs whether anti-amyloid therapeutic benefit would be expected to scale with ε4 carrier status — a question of immediate clinical relevance for therapeutic-trial enrichment strategies.

---

## 12. Reverse MR (Directionality of Causal Architecture)

### 12.1 Steiger filter

Per-SNP Steiger filter applied as in §10.3. SNPs where R²_outcome > R²_exposure are flagged as potentially reverse-causal and re-tested with the per-SNP filter.

### 12.2 Formal reverse MR for AD-spectrum phecodes

For phecodes F00, F03, G30, G31, run a formal reverse MR with:
- **Reverse exposure:** Bellenguez et al. 2022 AD GWAS (N=487,511 EUR; DOI 10.1038/s41588-022-01024-z).
- **Reverse outcome:** Ali 2023 multi-ethnic amyloid-PET.

Same instrument-construction protocol (§6) on AD as the reverse exposure. Reverse-MR IVW estimate compared to forward-MR IVW: H3 supported if forward effect magnitude > reverse effect magnitude AND Steiger directionality favours forward.

---

## 13. Mediation Analysis (MVMR via Cognitive Function / AD)

### 13.1 Trigger

For any non-AD-spectrum FDR-significant phecode in H1 or H2, run a mediation analysis to estimate how much of the amyloid → outcome effect is direct vs mediated through AD or cognitive function.

### 13.2 Method

Multivariable MR (Sanderson et al. 2019) with two or three exposures:
1. Genetically-proxied amyloid-PET (Ali 2023 instruments).
2. AD risk (Bellenguez 2022 instruments) OR general cognitive function (Davies et al. 2018 *Nat Commun*; DOI 10.1038/s41467-018-04362-x).
3. (If both AD and cognitive function show MVMR effects, sequential MVMR conditioning out the dominant mediator.)

Implementation: `MendelianRandomization::mr_mvivw()` with MVMR-PRESSO sensitivity.

### 13.3 Decomposition

- **Total effect** = univariate amyloid → outcome IVW.
- **Direct effect** = MVMR amyloid coefficient conditional on AD (or cognition).
- **Indirect (mediated) effect** = total − direct.

Bootstrap confidence intervals for the indirect-effect estimate via 1,000 bootstrap resamples of the instrument set.

Two-step MR (sequential univariate amyloid → mediator → outcome) reported as a methodologically distinct sensitivity for mediation.

### 13.4 Interpretation classifications

Each mediated phecode classified:
- **Direct-dominant:** direct effect > 70% of total effect.
- **Mediated-dominant:** indirect effect > 70% of total effect.
- **Mixed:** neither component > 70%.

### 13.5 Tau-PET cascade two-step MR

The canonical AD pathological cascade is amyloid → tau → clinical phenotype. Extending the §14 amyloid-direct framework to the cascade level, a pre-registered two-step MR decomposes per-phecode effects into (a) amyloid → tau → outcome (tau-mediated) and (b) amyloid → outcome direct of tau.

#### 13.5.1 Method

Three-exposure MVMR with:
1. Genetically-proxied amyloid-PET (Ali 2023 instruments).
2. Genetically-proxied tau-PET (Nho 2024 ADNI tau-PET GWAS instruments).
3. Genetically-proxied AD risk (Bellenguez 2022 instruments, retained from §13 to control for the residual AD signal).

Implementation: `MendelianRandomization::mr_mvivw()` with conditional F-statistics per exposure and MVMR-PRESSO outlier sensitivity.

#### 13.5.2 Decomposition

- **β_amyloid|tau,AD** = MVMR coefficient on amyloid-PET conditional on tau and AD = amyloid effect on outcome not mediated through tau or AD.
- **β_tau|amyloid,AD** = MVMR coefficient on tau-PET conditional on amyloid and AD = tau effect on outcome not mediated through amyloid or AD.
- **Indirect (amyloid → tau → outcome) effect** = β_amyloid × β_outcome_per_tau, computed via the product-of-coefficients method with bootstrap CI (1,000 resamples).

#### 13.5.3 Per-phecode cascade classification

| β_amyloid|tau,AD | β_tau|amyloid,AD | Indirect (amyloid → tau → outcome) | Cascade interpretation |
|---|---|---|---|
| significant | significant | significant | Amyloid acts via tau AND has direct amyloid-tau-independent effect |
| significant | null | null | Amyloid acts directly, independent of tau |
| null | significant | significant | Effect is tau-mediated; amyloid acts only via tau |
| null | null | null | No genetic signal on this outcome from either pathway |

#### 13.5.4 Honest power statement

Tau-PET GWAS (Nho 2024 ADNI; N≈3,000) is substantially smaller than the amyloid-PET exposure. Three-exposure MVMR with one weak instrument inflates standard errors substantially. The cascade decomposition is reported with explicit power caveats and the §10.5 non-linear MR results inform whether the cascade-classification verdicts are achievable at current GWAS sample sizes.

#### 13.5.5 Novelty statement

This is the first pre-registered application of the amyloid → tau → outcome cascade decomposition at the phenome-wide MR scale. The conventional MR literature treats amyloid and tau as separate exposures; this design integrates them into the canonical pathological cascade as a single MVMR framework.

---

## 14. Amyloid-Direct Effects via Multivariable MR (AD as Co-Exposure)

For each FinnGen R12 phecode, estimate the amyloid-direct effect using multivariable MR with two simultaneous exposures:
1. Genetically-proxied amyloid-PET (Ali 2023 instruments).
2. Genetically-proxied AD risk (Bellenguez 2022 instruments).

### 14.1 Primary statistic — MVMR amyloid-direct coefficient

The MVMR coefficient on amyloid-PET (β_amyloid|AD) is the amyloid-direct effect on the outcome, conditional on AD. This is the statistically principled decomposition: it accounts for the covariance between the two exposure estimates (which share cohorts and overlapping genetic architecture), unlike a simple difference of univariate estimates.

Implementation: `MendelianRandomization::mr_mvivw()` per (phecode × MVMR-exposure-set). Conditional F-statistics reported per exposure (Sanderson et al. 2019). MVMR-PRESSO (Sanderson et al. 2020) applied as outlier-detection sensitivity.

Note: §14 differs in scope from the Mediation Analysis in §13. §13 runs MVMR for non-AD-spectrum FDR-significant hits only, framed as a mediation decomposition (direct vs indirect). §14 runs the same MVMR estimator across the full 1,574-phecode universe to identify which outcomes show amyloid-direct effects independent of AD, regardless of FDR status in the univariate sweep.

### 14.2 Descriptive sanity check — β_residual

As a descriptive complement, compute β_residual = β_amyloid − β_AD using the univariate IVW estimates from each exposure. Bootstrap-resample both instrument sets jointly (1,000 resamples) for the bootstrap CI of β_residual.

β_residual is reported alongside β_amyloid|AD as a descriptive check. The two estimators should agree in direction for amyloid-direct effects; disagreement (β_amyloid|AD positive but β_residual null, or vice versa) is interpretable evidence that the simple subtraction is unreliable in this setting (typically due to instrument-sharing or sample overlap) and the MVMR estimate is the operative one.

### 14.3 Per-phecode classification

Classification is based on the MVMR primary statistic, with descriptive β_residual reported alongside:

| β_amyloid|AD (MVMR) | β_AD|amyloid (MVMR) | Interpretation |
|---|---|---|
| significant | null | Amyloid-direct effect distinct from AD pathway |
| significant | significant (concordant direction) | Both amyloid and AD contribute independently |
| null | significant | Effect runs via AD pathway; amyloid contributes nothing independent of AD |
| null | null | No genetic signal on this outcome from either exposure |

This decomposition directly answers "does amyloid have effects beyond what is captured by clinical AD?" — the central translational question for the anti-amyloid therapeutic label.

### 14.4 Pre-symptomatic vs symptomatic amyloid effects

The §14 MVMR amyloid-direct estimate β_amyloid|AD captures amyloid effects not mediated through clinical AD diagnosis — in other words, effects that operate in the pre-symptomatic amyloid-positive state. To make this interpretation explicit and to address the central anti-amyloid-therapeutic-timing question, two complementary analyses are pre-registered:

#### 14.4.1 Cognitively-normal-only amyloid exposure (contingent)

If Ali 2023 provides a cognitively-normal-only sub-meta-analysis (i.e., excluding AD-case-enriched cohorts and restricted to cognitively-unimpaired participants), instruments are re-extracted from this sub-cohort and a parallel MR sweep is run. The cognitively-normal-only exposure isolates the pre-symptomatic amyloid signal at the source.

**Contingency:** sub-cohort availability is verified from the Ali 2023 supplementary data at Phase 1 closure. If not available, this analysis is omitted and reported as "data unavailable" in the manuscript.

#### 14.4.2 CSF Aβ42–conditioned MVMR (always run)

If Ali 2023 cognitively-normal-only sumstats are unavailable, the pre-symptomatic decomposition is operationalised differently: a three-exposure MVMR with (Ali 2023 amyloid-PET, Jansen 2022 CSF Aβ42, Bellenguez 2022 AD). CSF Aβ42 declines earlier in the disease cascade than amyloid-PET positivity, so the CSF-Aβ42-conditioned amyloid-PET coefficient (β_amyloid|CSF-Aβ42,AD) approximates the symptomatic / late-stage amyloid effect, whereas the CSF-Aβ42 coefficient conditional on amyloid-PET approximates the pre-symptomatic effect.

This is a methodological approximation, not a clean separation; it is reported with the limitation explicitly noted.

#### 14.4.3 Decision-context interpretation

For each phecode, the pre-symptomatic vs symptomatic decomposition maps to one of three anti-amyloid-therapeutic decision contexts:

| Pre-symptomatic | Symptomatic | Decision-context implication |
|---|---|---|
| significant | significant | Effect operates across the full amyloid-positive trajectory; therapeutic intervention at any stage may benefit |
| null | significant | Effect emerges only in the symptomatic stage; pre-symptomatic intervention may not affect this outcome |
| significant | null | Effect operates pre-symptomatically and stabilises; later intervention may be too late for this outcome |
| null | null | No genetic signal on this outcome at either stage |

#### 14.4.4 Novelty statement

The pre-symptomatic vs symptomatic decomposition operationalises, for the first time at the phenome-wide MR scale, the central anti-amyloid-therapeutic-timing question. The decision-context table provides a direct mapping from genetic findings to therapeutic-timing recommendations per outcome.

---

## 15. Cross-Ancestry Analysis

### 15.1 Trans-ancestry meta-regression (MR-MEGA)

For each FDR-significant phecode, run MR-MEGA (Mägi et al. 2017; DOI 10.1093/hmg/ddx280) across all available ancestry-stratified instrument arms:

- Ali 2023 NHW (N=11,556)
- Ali 2023 ASN (N=1,494)
- Ali 2023 AFR (N=359)
- Kim 2025 EUR+EAS combined (if obtained)

Outputs reported:
- Pooled trans-ancestry effect estimate.
- Residual heterogeneity p-value (p_anc) testing whether effect varies along the ancestry axis.
- Ancestry-axis correlated component (PC1, PC2 of population genetic structure).

### 15.2 Outcome-side cross-ancestry replication

For FDR-significant phecodes, replicate using:
- BBJ matched phenotype sumstats (where available; phenotype-matching best-effort).
- Pan-UKB AFR matched phenotype sumstats (where available).

H4 supported if ≥ 1 non-APOE hit replicates at nominal p < 0.05 with concordant direction in at least one cross-ancestry outcome panel.

### 15.3 Honest power statement and primary-verdict reframing for ancestry analyses

Ali 2023 ASN N=1,494 and AFR N=359 are small. Pre-registered expectation: ancestry-stratified instrument panels will frequently have n_IV < 5 outside APOE, and formal heterogeneity tests will be underpowered.

For the H5 primary verdict, the operative test is **directional consistency**, not statistical heterogeneity:

- For each FDR-significant non-APOE hit in the EUR-heavy primary, the sign of the IVW point estimate in each ancestry arm (where n_IV ≥ 2) is recorded.
- Directional-consistency verdict per hit: "concordant" (≥ 2 of 3 arms agree in sign), "discordant" (≥ 2 arms disagree), "inconclusive" (insufficient n_IV in ≥ 2 arms).
- A non-APOE hit with concordant direction across NHW + ASN + AFR is the most defensible cross-ancestry evidence currently obtainable.

The formal MR-MEGA trans-ancestry meta-regression (§15.1) is reported as the secondary, heterogeneity-focused inference, with explicit acknowledgement that p_anc estimates are power-limited at these sample sizes.

This reframing — directional concordance primary, heterogeneity secondary — matches the achievable inference at the available cohort sizes and avoids treating an underpowered heterogeneity-null result as evidence of cross-ancestry agreement.

---

## 16. Sex-Stratified Sub-Analysis

### 16.1 Sex-stratified instrument extraction

- Ali 2023 female arm (N=5,195) and male arm (N=4,625) clumped independently.
- Sex-stratified IVW per phecode using the pooled multi-ethnic outcome.

### 16.2 Between-sex heterogeneity

Per phecode, Cochran's Q across the female and male IVW estimates with corresponding p-value (Q_sex). BH-FDR within the FDR-significant H1 hit-restricted sub-family.

### 16.3 Pre-specified direction expectation

For AD-spectrum phecodes (F00, F03, G30, G31) and psychiatric phecodes (F00–F69 ICD-10 chapters), the documented sex-difference literature predicts greater magnitude of effect in females. The sex-stratified analysis is pre-registered as a directional test (one-sided after the across-sex Q-test) for these chapters; the direction expectation is reversed (or absent) for cardiovascular and traumatic-injury chapters.

---

## 17. Colocalisation

### 17.1 Method

Bayesian colocalisation via `coloc::coloc.abf()` and `coloc::coloc.susie()` for each (exposure locus × outcome locus) pair, for the top 3 non-APOE exposure loci per FDR-significant phecode.

### 17.2 Loci definition

Per FDR-significant phecode, identify the top 3 non-APOE exposure SNPs by exposure-GWAS p-value within the APOE-excluded instrument panel. For each, extract full summary statistics within ±500 kb (1 Mb total window) for both exposure and outcome GWASs.

### 17.3 Priors

| Analysis | p1 (exposure-only) | p2 (outcome-only) | p12 (shared) |
|---|---|---|---|
| Primary | 1 × 10⁻⁴ | 1 × 10⁻⁴ | 1 × 10⁻⁵ |
| Sensitivity 1 (conservative) | 1 × 10⁻⁴ | 1 × 10⁻⁴ | 1 × 10⁻⁶ |
| Sensitivity 2 (liberal) | 1 × 10⁻⁴ | 1 × 10⁻⁴ | 1 × 10⁻⁴ |

### 17.4 SuSiE multi-causal-variant fine-mapping

For loci where standard coloc.abf returns intermediate H4 posterior (0.5 < H4 < 0.8), SuSiE multi-causal-variant fine-mapping is applied as a sensitivity to distinguish "no shared signal" (H3) from "multiple distinct causal variants at the locus".

### 17.4.1 Prior-sensitivity sweep

The choice of p12 (shared causal-variant prior) is methodologically contested. In addition to the three pre-registered point priors in §17.3, `coloc::sensitivity()` is invoked for every locus with H4 < 0.8 to produce a prior-sensitivity sweep showing how H4 varies across the full prior-probability space (p12 ranging from 10⁻⁷ to 10⁻³). The reported colocalisation result is considered robust if H4 > 0.5 across at least 80% of the prior-space tested. Sensitivity-sweep heatmaps are deposited as supplementary material per FDR-significant locus.

### 17.5 Tissue-specific eQTL colocalisation

For each (exposure SNP, outcome) significant coloc, additionally colocalise the exposure signal against:
- **GTEx v8** brain tissues (Cortex; Frontal Cortex BA9; Hippocampus; Hypothalamus; Substantia nigra; Cerebellum) and whole-blood.
- **eQTLGen** consortium whole-blood eQTL (N≈31,684).

Reported as a heatmap of H4 per (exposure SNP × tissue × candidate gene).

### 17.6 Interpretation

| H4 posterior | Interpretation |
|---|---|
| H4 > 0.8 | Strong evidence of shared causal variant; supports causal interpretation |
| H3 > 0.8 | Evidence of distinct causal variants at the locus; MR estimate may be LD-confounded |
| 0.5 < H4 < 0.8 | Suggestive; trigger SuSiE multi-causal fine-mapping |
| H4 < 0.5 | No evidence of colocalisation |

### 17.7 Multi-ancestry SuSiE fine-mapping for cross-ancestry hits

For phecodes that pass H5 directional-consistency (i.e., concordant non-APOE effect direction across ≥ 2 of NHW / ASN / AFR arms; §15.3) AND show colocalisation H4 > 0.5 at a non-APOE locus, pre-register multi-ancestry SuSiE fine-mapping at that locus.

#### 17.7.1 Method

Multi-ancestry SuSiE (Yuan et al. 2023 *PLOS Genet*; or MsCAVIAR-style multi-ancestry fine-mapping) jointly fine-maps the locus across ancestry-stratified GWASs using ancestry-specific LD reference panels (1000G EUR, EAS, AFR). The output is a posterior credible set of putative causal variants, with per-variant ancestry-specific posterior inclusion probabilities (PIPs).

#### 17.7.2 Interpretation

- **Shared causal variant across ancestries:** the credible set contains the same variant(s) with high PIP across ≥ 2 ancestries. Supports universal causal architecture; the genetic effect transports across populations.
- **Ancestry-specific causal variants:** the credible set differs by ancestry. Supports ancestry-specific causal architecture; therapeutic targeting strategies may need ancestry-specific adaptation.
- **Ambiguous:** credible sets overlap but with low PIP. Reported as inconclusive.

#### 17.7.3 Novelty statement

Pushes the analysis from causal-effect estimation (the H5 cross-ancestry MR) to causal-variant identification across ancestries. For any cross-ancestry-consistent amyloid effect identified in this study, the multi-ancestry fine-mapping output is the highest-resolution causal-architecture statement the available data permit.

#### 17.7.4 Honest power statement

Ali 2023 ASN N=1,494 and AFR N=359 are small; multi-ancestry fine-mapping is heavily power-limited at these sample sizes. The fine-mapping result is reported with explicit per-ancestry credible-set size and PIP uncertainty.

---

## 18. LDSC Genetic-Correlation Cross-Tabulation

For every phecode in the 1,574-phecode universe, the LDSC genetic correlation (rg) between Ali 2023 amyloid-PET and the phecode is reported alongside the MR IVW estimate. The cross-tabulation of (MR significant × rg significant) provides complementary causal-inference signal beyond MR alone:

| MR | rg (\|rg\| > 0.1 at p_FDR < 0.05) | Interpretation |
|---|---|---|
| significant | significant | Strong shared causal architecture; MR causal interpretation supported by genetic correlation |
| significant | null | Possibility of horizontal pleiotropy via instruments OR weak rg power; cautious interpretation |
| null | significant | Shared genetic architecture without detectable causal signal; could be limited MR power or genuinely-shared-but-non-causal genetic component |
| null | null | No genetic relationship detected at either level |

Reported in Supplementary Material as a cross-tab matrix and as a scatter of MR β vs rg per phecode.

---

**Part E — Quality controls**

---

## 19. Sample-Overlap Audit

### 19.1 Pre-registered audit

Sample overlap between exposure and outcome GWASs biases MR estimates toward the observational confounded estimate. Pre-registered audit:

| Comparison | Expected overlap | Audit method |
|---|---|---|
| Ali 2023 ↔ FinnGen R12 | Near-zero (Ali 2023 cohorts are US/EU/AU; FinnGen is Finnish) | Cohort-level overlap check via Ali 2023 supplementary Table 1 and FinnGen R12 cohort documentation |
| Ali 2023 ↔ BBJ | Near-zero | Cohort-level cross-reference |
| Ali 2023 ↔ Pan-UKB AFR | Possibly non-zero | Cohort-level cross-reference |
| Ali 2023 ↔ Bellenguez 2022 AD | Non-zero (ADNI, ROS/MAP appear in both) | Cohort-level cross-reference + sensitivity analysis using Bellenguez 2022 non-FinnGen non-ADNI sub-meta-analyses if available |
| Kim 2025 ↔ FinnGen R12 | Near-zero | Cohort-level cross-reference |

### 19.2 Overlap correction — MRlap as primary tool

For any non-zero overlap identified, the **MRlap** package (Mounier & Kutalik 2023 *PLOS Genet*) is the pre-registered primary correction tool. MRlap jointly models exposure and outcome GWAS summary statistics under known or estimated sample overlap and produces overlap-corrected MR estimates with appropriate standard errors.

For settings where MRlap is inapplicable (e.g., reverse-MR with Bellenguez 2022 as the reverse exposure, where ADNI / ROS-MAP overlap with Ali 2023 is non-trivial), sub-cohort sensitivity using Bellenguez 2022 non-ADNI / non-ROS-MAP sub-meta-analyses (where available) is the alternative correction route.

The MRlap overlap-corrected estimate is reported alongside the uncorrected IVW for every (exposure, outcome) pair where overlap > 0%. Bias direction and magnitude documented in Limitations.

---

## 20. Falsification Checks (Pre-specified Negative and Positive Controls)

### 20.1 Positive controls

| Phecode | Endpoint | Expected | If null |
|---|---|---|---|
| `G30_ALZHEIMER` | Alzheimer's disease | Strong positive signal (OR > 1.30) | Suggests instrument failure or harmonisation error; full audit before any reporting |
| `F00_DEMENTIA_VASCULAR` | Vascular dementia | Positive signal | Less stringent; not invalidating |
| `F03_DEMENTIA_UNSPECIFIED` | Unspecified dementia | Positive signal | Less stringent |

### 20.2 Within-phenome negative controls

| Phecode (or category) | Endpoint | Expected |
|---|---|---|
| Traumatic fracture phecodes | Trauma | Null (no causal link to amyloid) |
| Acute trauma chapter (S00–S99) | Traumatic injuries | Null |
| Pregnancy / puerperium chapter (O15) | O-chapter endpoints | Null |
| Refractive eye errors (H7) | Eye/refractive | Null |

**Excluded from within-phenome negative-control panel:** Down syndrome / chromosomal disorders (Q-chapter) — Down syndrome IS associated with early-onset amyloid and is excluded for this biological reason.

### 20.2.1 High-power external negative controls (population-stratification diagnostic)

In addition to the within-FinnGen negative controls above, the amyloid-PET instruments are tested against external high-power GWASs of traits with massive genetic signal but near-zero biological relationship to amyloid burden. A "significant causal effect" on any of these traits flags systematic confounding (population stratification, instrument-side pleiotropy, or harmonisation error) that within-FinnGen negative controls may not detect.

| External outcome GWAS | N | Use |
|---|---|---|
| Yengo et al. 2022 Height GWAS (*Nature*; DOI 10.1038/s41586-022-05275-y) | ~5,400,000 (multi-ancestry; EUR primary) | Population-stratification diagnostic |
| Hysi et al. 2018 Eye colour GWAS (*Nat Genet*; DOI 10.1038/s41588-018-0049-4) — if accessible at deposition time | ~195,000 | Secondary stratification diagnostic |

These high-power outcomes have orders of magnitude greater statistical power than typical FinnGen phecodes; a spurious "causal effect" of amyloid on Height would be detectable at very small effect sizes and is therefore a sensitive diagnostic for cryptic confounding. Pre-registered threshold: any p < 0.05 effect of amyloid-PET on Height or Eye Colour triggers a methodological audit before reporting the primary findings.

### 20.2.2 Negative-control failure protocol

If any within-phenome OR external negative-control passes the pre-registered failure threshold, the testing methodology is audited for systematic confounding before the primary findings are reported. Audit checklist:

1. Re-verify instrument harmonisation (allele coding, build, palindromic handling).
2. Re-check LDSC intercept on the offending negative-control outcome.
3. Re-extract instruments with stricter palindrome exclusion and re-run.
4. If audit fails to resolve the negative-control signal, primary findings are reported with explicit acknowledgement of the negative-control failure and the systematic-confounding implication.

### 20.3 Effect-magnitude bound check

The strongest reported MR effect across the 1,574-phecode sweep should not exceed the strongest reported observational amyloid-PET ↔ outcome effect by more than 3×. A reported MR OR > 3× the strongest observational effect triggers a methodological audit before publication.

---

**Part F — Safety and clinical translation**

---

## 21. Phecodes of Clinical Interest (anti-amyloid therapeutic safety panel)

In addition to the primary 1,574-phecode sweep, the following phecodes are pre-specified as **phecodes of clinical interest** because they map to clinically-relevant safety endpoints for anti-amyloid therapeutic decisions. The actual FinnGen R12 phenocodes for these conditions are documented at analysis time in `data/processed/clinical_interest_phecodes.csv`; the table below specifies the clinical-relevance rationale:

| Condition (ICD-10) | Clinical relevance |
|---|---|
| Intracerebral haemorrhage (I61) | Direct ARIA-H concern; relevant to anti-amyloid bleeding risk |
| Ischaemic cerebrovascular accident (I63) | Vascular event in amyloid-positive |
| Cerebral oedema (G93.6) | Direct ARIA-E concern |
| Subarachnoid haemorrhage (I60) | Vascular safety endpoint |
| Cerebral infarction with concurrent imaging abnormality | ARIA-related concern |
| Anaphylactic reactions (T78.0–T78.2) | Anti-amyloid administration risk |
| Depression (F32) | Common AD comorbidity; off-target signal worth checking |
| Bipolar disorder (F31) | Psychiatric off-target signal |
| Schizophrenia (F20) | Psychiatric off-target signal |

For phecodes of clinical interest:
1. The equivalence test (§22) is reported alongside the conventional null-hypothesis test.
2. The GWAS-by-subtraction (§14) decomposition is reported.
3. Effect-size translation to lecanemab-equivalent clearance (§23) is reported.

Reporting these specifically — rather than waiting for them to emerge from the 1,574-phecode sweep — anchors the clinical-translation framing to a pre-registered safety-relevance commitment.

---

## 22. Equivalence Testing for Safety-Relevant Phecodes

For the phecodes of clinical interest pre-specified in §21, pre-register a two-one-sided test (TOST) equivalence framework alongside the conventional null-hypothesis significance test.

### 22.1 Equivalence margin

- Δ_eq = 0.05 (absolute log-OR magnitude). This corresponds to OR ≈ 1.05 per 1-SD genetically-proxied amyloid-PET increase.
- Rationale: an effect smaller than OR 1.05 is below the magnitude of risk-factor effects typically considered clinically actionable in genetic epidemiology, and below the noise floor of observational risk-stratification practice.

### 22.2 Test

For each safety-relevant phecode, simultaneously test:
1. H_low: IVW point estimate is significantly less than +Δ_eq (one-sided).
2. H_high: IVW point estimate is significantly greater than −Δ_eq (one-sided).

Both at α = 0.025 (two one-sided tests; family-wise α = 0.05). If both reject, the phecode is declared **safety-equivalent to no effect** at the pre-registered margin.

### 22.3 Interpretation

- **Equivalence-positive:** "No evidence of clinically-meaningful effect at the pre-registered safety margin." This is a distinct verdict from "no statistically-significant evidence" — it actively rules out the existence of a meaningful effect.
- **Equivalence-negative + conventional null:** "Cannot rule out effect of clinically-meaningful magnitude; cannot establish presence either." This is the underpowered-cannot-conclude verdict in the safety framing.
- **Equivalence-negative + conventional positive:** standard "evidence of effect" verdict.

The equivalence-test verdicts are reported as a separate column in the safety-phecode results table.

---

## 23. Effect-Size Translation to Therapeutic Context

The MR effect estimate is per 1-SD genetically-proxied increase in cortical amyloid-PET signal. For clinical translation, the manuscript reports the effect on a Centiloid scale using the Centiloid Working Group rescaling formula (Klunk et al. 2015 *Alzheimers & Dementia*; DOI 10.1016/j.jalz.2014.07.003).

### 23.1 Foundational caveat — lifetime exposure vs therapeutic intervention

**This caveat applies to every translation statement below and must be carried through in any manuscript reporting.**

The lifetime cumulative genetic exposure captured by MR instruments is NOT directly comparable to an 18-month therapeutic amyloid clearance:

- The genetic effect reflects decades of cumulative amyloid burden modulation from conception onward.
- A therapeutic intervention reduces amyloid burden at a single time-point in adulthood, with potentially non-linear pharmacodynamics, reversible binding kinetics, and downstream pathological cascades already in progress.
- The two are not interchangeable, and a 25-Centiloid genetic increase does NOT equate to "what would happen if lecanemab were withheld" in a therapeutic-effect sense.

The genetic effect is an **etiological causal effect**, not a **therapeutic intervention effect**. Every reported translation statement carries this distinction explicitly.

### 23.2 Centiloid conversion

Conversion formula (pre-registered): 1 SD of Ali 2023's pooled amyloid-PET phenotype is computed from the per-cohort SDs documented in Ali 2023 supplementary materials. The resulting "1-SD in Centiloid units" is reported for the manuscript audience. The Centiloid scale is a standardised quantification of cortical amyloid burden harmonised across PET tracers (Klunk et al. 2015).

### 23.3 Lecanemab-equivalent magnitude (descriptive only, under §23.1 caveat)

Lecanemab CLARITY-AD (van Dyck et al. 2023 *NEJM*; DOI 10.1056/NEJMoa2212948) reported a mean ~25-Centiloid amyloid reduction at 18 months of treatment. For each FDR-significant phecode, the effect of a 25-Centiloid genetically-proxied increase is reported alongside the per-1-SD effect — purely as a magnitude reference for clinical readers, NOT as a claim that a 25-Centiloid therapeutic reduction would produce a phenotypic change of the same magnitude. The §23.1 caveat is repeated in every such translation statement in the manuscript.

### 23.4 Post-marketing surveillance triangulation (forward-looking commitment)

Lecanemab received FDA traditional approval in July 2023; donanemab in July 2024. Real-world adverse-event reporting is accumulating in public pharmacovigilance databases (FAERS in the US; DAEN in Australia; EudraVigilance in the EU). For any FDR-significant amyloid → outcome MR hit identified in this study, the manuscript will pre-specify a triangulation step:

#### 23.4.1 Method

At the time of manuscript submission:
1. Query FAERS, DAEN, and EudraVigilance for lecanemab and donanemab adverse-event reports.
2. For each FDR-significant MR hit (and for each phecode of clinical interest in §21), compute the disproportionality statistic (reporting odds ratio, ROR; or proportional reporting ratio, PRR) for the corresponding adverse-event term against background reporting.
3. Concordance check: does the post-marketing signal (ROR direction and effect size) agree with the genetic prediction (MR direction and effect size)?

#### 23.4.2 Pre-registered interpretation

| MR direction | Post-marketing signal | Triangulation verdict |
|---|---|---|
| Positive (amyloid → ↑ outcome) | ROR > 1 + 95% CI excludes 1 | Triangulation supports causal interpretation |
| Positive | ROR null or negative | Genetic effect not yet emerging in post-marketing data; possibly underpowered or possibly not therapeutically operative |
| Negative (amyloid → ↓ outcome) | ROR < 1 + 95% CI excludes 1 | Triangulation supports protective causal effect; therapeutic clearance may reduce outcome |
| Null MR | ROR > 1 | Therapeutic adverse event not predicted by lifetime-burden genetics; could be pharmacological off-target rather than amyloid-mechanism-mediated |

#### 23.4.3 Honest power caveat

Lecanemab and donanemab uptake is still modest as of 2026-05-19; FAERS / DAEN power for most outcomes is limited. The triangulation is framed as a forward-looking commitment: the manuscript will be drafted at a time when at least 2 years of accumulated post-marketing data are available. If FAERS / DAEN remain underpowered at manuscript submission, the triangulation section is reported as "data accumulating; signal-to-noise inadequate for triangulation at this time" with the analysis re-runnable as data accumulate.

#### 23.4.4 Novelty statement

Pre-registering a forward-looking triangulation between genetic causal inference and post-marketing therapeutic surveillance is, to our knowledge, unprecedented in the MR-PheWAS literature. The commitment makes the registration's conclusions falsifiable not only by future genetic data but also by accumulating real-world therapeutic evidence.

---

**Part G — Inference framework**

---

## 24. Multiple Testing

### 24.1 Primary correction

Benjamini-Hochberg (BH) FDR within each of 19 phecode chapters (per the companion phecode-filter CSV `chapter_prefix` column) at family-wise q < 0.05. Within-chapter correction is more powerful than global BH because prior probability of association varies by chapter (e.g., AD/neurology chapter has higher prior than musculoskeletal chapter).

### 24.2 Sensitivity correction

- **Global BH** across the full 1,574-phecode sweep (more conservative; reported alongside).
- **Bonferroni** at α = 0.05 / (1,574 × 2 instrument panels) = 1.59 × 10⁻⁵ (most conservative; reported alongside).
- **Benjamini-Yekutieli (2008) hierarchical FDR** to control FDR under within-chapter correlation (handles parent/child phecode dependence formally).
- **Monte Carlo–calibrated FDR threshold** per chapter (1,000 permutations of phecode-outcome assignments) reported as a tertiary correction-aware comparison.

All correction outputs reported in the results parquet; the BH within-chapter is the primary inference framework.

### 24.3 No post-hoc threshold adjustment

Post-hoc α adjustment based on the observed p-distribution is explicitly forbidden. The four correction methods above are the only thresholds applied.

---

## 25. Power Calculation and Atlas

### 25.1 Method

Analytic two-sample MR power calculation per **Burgess (2014) *Int J Epidemiol* 43(3):922–929** (DOI 10.1093/ije/dyu005), computed for every (phecode × instrument-panel) combination:

- α = 0.05 / (1,574 × 2) = 1.59 × 10⁻⁵ (Bonferroni-conservative) for primary power.
- α = 0.05 for the more liberal report.
- Power target: 80%.
- Variance explained by instruments (R²_exposure) computed from observed F-statistics.
- Case/control counts from FinnGen R12 manifest.

### 25.2 Reported per phecode

- Minimum detectable OR (MDOR) at 80% power at α = 1.59 × 10⁻⁵.
- Minimum detectable OR at 80% power at α = 0.05.
- Required outcome case count to detect OR = 0.85 at 80% power.
- Power status classification: **adequately powered** (MDOR ≤ 1.20) / **moderately powered** (1.20 < MDOR ≤ 1.50) / **underpowered** (MDOR > 1.50).

### 25.3 Standalone deliverable

The power atlas is the primary deliverable of Pivot A (§3.2 / §27). Even in the primary-analysis-supported scenario, it is reported as Supplementary Material.

### 25.4 Bootstrap-based empirical power

Beyond the analytic Burgess formula, for FDR-significant H1 hits a bootstrap-based empirical power calculation is performed: 1,000 bootstrap resamples of the instrument set; per-resample IVW estimated; fraction of resamples returning p < 1.59 × 10⁻⁵ gives the empirical power conditional on the observed instrument F-distribution.

---

## 26. Decision Rules and Verdict Tiers

### 26.1 Rationale

Pure p-value cutoffs are insufficient for confirmatory inference. Following recent best-practice in MR (Burgess et al. 2023 *eLife*), each hit's verdict is determined by a combination of:
1. Direction of effect across primary + sensitivity methods.
2. Confidence interval exclusion of the null.
3. Effect-size magnitude (effect-size tier).

### 26.2 Verdict definitions

A phecode passes FDR q < 0.05 AND meets the §10.2 concordance threshold AND meets the §26.3 effect-size tier requirement to be declared a "**robust causal evidence**" finding.

A phecode meeting FDR q < 0.05 but failing concordance or effect-size tiers is reported as "**suggestive; requires replication**".

A phecode that fails FDR but meets nominal p < 0.05 with concordance across ≥ 3 of 4 (IVW + WM + WMode + MR-PRESSO) is reported as "**directionally suggestive; hypothesis-generating**".

### 26.3 Effect-size tiers (per 1-SD amyloid-PET increase)

| Tier | OR range (binary outcome) | Beta range (continuous) | Manuscript label |
|---|---|---|---|
| **Strong** | OR ≥ 1.50 or ≤ 0.67 | \|β\| ≥ 0.30 SD | "Strong evidence" |
| **Moderate** | OR 1.20–1.50 or 0.67–0.83 | \|β\| 0.15–0.30 SD | "Moderate evidence" |
| **Weak support** | OR 1.10–1.20 or 0.83–0.91 | \|β\| 0.07–0.15 SD | "Weak support" |
| **Negligible** | OR 1.00–1.10 or 0.91–1.00 | \|β\| < 0.07 SD | Not reported as a finding |

Effect-size tier "Weak support" is the minimum threshold for any positive verdict.

### 26.4 Honest power-tier cross-reference

For every reported verdict, the power tier from §25.2 is cross-referenced. A phecode where the verdict is "negative" (no FDR hit) AND the power tier is "underpowered" is explicitly framed as **"underpowered; cannot conclude"** rather than as a negative finding.

### 26.5 "What would change my mind" — pre-registered retraction triggers

Following Burgess et al. 2023 *eLife* recommendations, the following pre-specified observations would lead the investigator to substantially revise or retract the primary findings:

1. **Positive control (G30 AD) returns null in the APOE-inclusive panel** — indicates instrument-extraction or harmonisation failure; primary findings cannot be trusted; methodological audit precedes any reporting.
2. **Any negative-control phecode (§20.2) OR external negative-control trait (§20.2.1 Height / Eye colour) passes the pre-registered failure threshold** — indicates systematic confounding; audit precedes any reporting.
3. **Effect-magnitude bound check fails (MR estimate > 3× the strongest observational estimate for the same outcome)** — indicates instrument over-weighting or pleiotropy; specific hit requires methodological audit.
4. **MR-Egger intercept p < 0.05 for ≥ 50% of significant hits** — indicates widespread directional pleiotropy; primary IVW estimates are systematically biased; CAUSE / MR-RAPS estimates replace IVW as primary in the manuscript.
5. **Reverse MR (§12) returns AD → amyloid effect larger in magnitude than amyloid → AD** — indicates reverse causation; AD-spectrum primary findings downgraded.
6. **Sample-overlap audit (§19) reveals > 10% Ali 2023 ↔ FinnGen overlap** — primary findings re-estimated using overlap-correction methods (MRlap).
7. **Cross-ancestry replication (§15) shows opposite direction in ≥ 2 of 3 ancestry arms** — primary findings are EUR-bounded only; ancestry-divergent finding becomes the central novel claim.
8. **LDSC heritability of exposure h² < 0.05** — MR is fundamentally uninformative; PGS-MR sensitivity becomes the primary analytic frame.
9. **LDSC intercept of any positive-control phecode > 1.20** — population stratification cannot be ruled out; primary findings re-examined.
10. **Robust MR-Egger (Bowden 2018) returns null for ≥ 30% of IVW-positive hits** — indicates that the IVW positive findings are not robust to invalid-instrument inclusion.

Each retraction-trigger observation is documented in the amendment log with timestamp.

---

## 27. Pre-Specified Pivots (operational triggers)

### 27.1 Trigger conditions and deliverables

| Pivot | Trigger | Deliverable |
|---|---|---|
| A — Weak-instrument atlas | APOE-excluded n_IV < 5 OR mean F < 15 OR LDSC exposure h² < 0.05 | Per-phecode minimum-detectable-OR; future-GWAS sample-size recommendations; PGS-MR sensitivity becomes primary |
| B — APOE-pleiotropy atlas | H1 ≥ 5 hits BUT H2 returns 0–1 | Per-phecode APOE decomposition; ε4 MVMR; APOE pleiotropy classification (A/B/C/D) |
| C — Sex-difference atlas | H1 + H6 both supported with ≥ 3 sex-heterogeneous phecodes | Sex-stratified atlas; sex-difference mechanistic discussion |

### 27.2 No silent re-framing

The pivot framing is locked at this pre-registration. Any deviation — including running a "fourth pivot" or combining two — is recorded in the amendment log with the trigger condition, the deviation, and whether the deviation was identified before or after viewing FDR-significance results.

---

## 28. Pre-Specified Manuscript Framing Matrix

To prevent post-hoc narrative drift, the manuscript's Results-section framing is pre-specified for the most likely verdict combinations across H1, H2, H3, H4, H5, H6. Notation: ✓ = supported, ✗ = not supported, * = any.

| H1 | H2 | H3 | H4 | H5 | H6 | Manuscript framing |
|---|---|---|---|---|---|---|
| ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | Headline finding: APOE-independent causal effects on multiple non-AD phenotypes, replicated cross-ancestry and with sex-specific architecture. Strongest causal interpretation. |
| ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | Headline finding without sex difference. Sex-stratified analysis null; cross-ancestry positive. Sex-difference framing dropped. |
| ✓ | ✓ | ✓ | ✓ | ✗ | * | EUR-anchored finding. Cross-ancestry replication weak; discussion explicitly bounds generalisability to EUR. |
| ✓ | ✓ | ✗ | * | * | * | Forward-MR confirmed, reverse-MR ambiguous. Steiger fails for ≥ 1 AD phecode; reverse causation cannot be ruled out for AD-spectrum hits. Non-AD hits retained as primary findings. |
| ✓ | ✗ | * | * | * | * | APOE-pleiotropy paper (Pivot B activated). Primary findings are APOE-driven; APOE-pleiotropy atlas is the manuscript framing. |
| ✗ | * | * | * | * | * | Weak-instrument atlas (Pivot A activated). Power atlas + PGS-MR sensitivity is the primary deliverable. |
| any | any | any | any | any | ✓ (sex-difference primary) | Sex-difference paper (Pivot C activated). Sex-stratified findings dominate. |

### 28.1 Clinical-translation framing

Every verdict combination ends with a "Clinical implications" subsection in the Discussion that translates the finding into one of three pre-specified decision contexts:

1. **Anti-amyloid therapeutic label.** Does the finding inform off-target safety or screening indications for lecanemab/donanemab?
2. **Pre-symptomatic amyloid-positive intervention timing.** Does the finding shift the calculus for treating amyloid-positive cognitively-unimpaired individuals?
3. **Routine clinical screening for amyloid-PET.** Does the finding strengthen or weaken the case for population-level amyloid-PET screening?

The clinical-translation commitment is binding: even a null verdict (Pivot A activated) yields a concrete recommendation rather than just a methodological null.

---

**Part H — Reporting**

---

## 29. Reporting Standards

- **STROBE-MR** (Skrivankova et al. 2021 *JAMA*; DOI 10.1001/jama.2021.18236) — primary reporting checklist; completed and submitted with manuscript.
- **PRISMA-MR** where applicable for any meta-analytic component.
- All instrument extraction, harmonisation, and analytic decisions are logged with timestamps in `data/audit/` (see §36).

---

## 30. Pre-Registered Visualizations

The following manuscript figures are pre-registered as structural commitments (specific content depends on results):

1. **Figure 1.** Study schematic (instrument extraction → MR → sensitivity → colocalisation → replication).
2. **Figure 2.** Manhattan-style plot by phecode chapter showing -log10(p) for APOE-inclusive (top half) and APOE-excluded (bottom half).
3. **Figure 3.** Forest plot of all FDR-significant phecodes with IVW + sensitivity methods (Egger, WM, WMode, MR-PRESSO outlier-corrected, MR-RAPS, CAUSE).
4. **Figure 4.** APOE pleiotropy classification heatmap (Categories A/B/C/D) per FDR-significant phecode.
5. **Figure 5.** Cross-ancestry effect-estimate comparison per phecode of clinical interest.
6. **Figure 6.** Sex-stratified effect comparison per FDR-significant phecode.
7. **Figure 7.** Colocalisation H4 heatmap per (hit × top non-APOE locus × tissue from GTEx + eQTLGen).
8. **Figure 8.** Power atlas — per-phecode minimum-detectable-OR heatmap with FDR-significant hits overlaid.
9. **Figure 9.** β_residual (amyloid − AD) vs β_amyloid scatter, identifying phecodes with direct vs AD-mediated effects.
10. **Figure 10** (supplementary). LDSC rg vs MR estimate scatter with the 2×2 cross-tabulation flagged.
11. **Figure 11** (supplementary). Equivalence-test (TOST) results for the phecodes-of-clinical-interest panel.

Specific figure axis ranges, colour palettes, and statistical annotations are determined at manuscript drafting; the structural commitment above prevents fishing for the most-photogenic visualization.

---

**Part I — Limitations**

---

## 31. Ancestry and Generalisability Limits

### 31.1 Primary ancestry framing

The primary IVW analysis uses the multi-ethnic Ali 2023 cohort (predominantly European). FinnGen R12 outcomes are predominantly European Finnish. The primary findings are bounded as European-applicable; non-European generalisability rests on the H4 cross-ancestry replication (BBJ + Pan-UKB AFR + Kim 2025 if obtained) and the H5 ancestry-stratified analysis.

### 31.2 Limitations documentation

The manuscript will explicitly document:
1. Ali 2023 ancestry composition (NHW 86%, ASN 11%, AFR 3% of total N=13,409).
2. FinnGen Finnish-population-specific founder effects on phecode prevalence.
3. The expected weakness of ancestry-stratified instrument panels (ASN N=1,494, AFR N=359).
4. Absence of any large-scale East Asian or African amyloid-PET GWAS as of 2026-05-19.

---

## 32. Tracer and Cohort Heterogeneity

### 32.1 Tracer-pooling acknowledgement

Ali 2023 pools across amyloid-PET tracers including florbetapir, florbetaben, flutemetamol, and Pittsburgh Compound B (PiB). Different tracers have different cortical retention profiles and different relative affinities for diffuse vs neuritic amyloid plaques. The pooled effect estimate represents an averaged tracer-agnostic amyloid signal. This is acknowledged in the Limitations.

### 32.2 Cohort-meta-analysis sensitivity

If per-cohort sub-meta-analyses become available, a cohort-restricted sensitivity is performed. Pre-registered: this is conditional on data availability, not a primary commitment.

### 32.3 APOE ε4 enrichment in source cohorts

Many amyloid-PET cohorts (DIAN, ADNI) are enriched for APOE ε4 carriers because of recruitment for AD-spectrum research. This enrichment biases the APOE ε4 effect-size estimate in Ali 2023. The APOE ε4 dosage MVMR (§11.2) is reported with explicit acknowledgement that the ε4 effect-size estimate may overstate the population-level ε4 effect.

---

## 33. Survival Bias

Older participants who develop high amyloid burden may die before reaching outcome data collection in FinnGen (which has a defined registry follow-up window). This creates a competing-risks selection bias: outcomes correlated with mortality (e.g., cardiovascular events, cancer) may be selectively under-represented among high-amyloid individuals. Pre-register acknowledgement in Limitations. Where FinnGen provides per-endpoint mean age-at-diagnosis, the bias is described qualitatively; quantitative correction requires individual-level data and is out of scope.

---

## 34. Within-Family MR Limitation

Brumpton et al. 2020 *Nat Commun* (DOI 10.1038/s41467-020-19147-4) and Howe et al. 2022 *Nat Genet* (DOI 10.1038/s41588-022-01062-7) demonstrated that conventional population-level MR is vulnerable to population stratification, assortative-mating effects, and dynastic genetic effects. The gold-standard mitigation — within-family MR using sibling-pair genotype data — is unavailable for amyloid-PET (no within-family amyloid-PET GWAS exists). Reported as a fundamental methodological limitation in the manuscript Discussion.

---

**Part J — Process and governance**

---

## 35. Software and Reproducibility

| Software | Version | Role |
|---|---|---|
| R | ≥ 4.3.0 | Primary analysis language |
| `TwoSampleMR` | Latest CRAN/GitHub at analysis start; locked in `renv.lock` | Primary MR estimation |
| `MendelianRandomization` | Latest CRAN | MVMR; secondary MR estimators |
| `MRPRESSO` | GitHub `rondolab/MR-PRESSO` | Outlier detection |
| `MVMR` | GitHub for MVMR-PRESSO | MVMR outlier detection |
| `RadialMR` | Latest CRAN | Outlier detection |
| `mr.raps` | GitHub `qingyuanzhao/mr.raps` | Weak-instrument robustness |
| `cause` | GitHub `jean997/cause` | Correlated pleiotropy robustness |
| `MR-MEGA` | Software from Mägi et al. 2017 | Trans-ancestry meta-regression |
| `coloc` | Latest CRAN | Colocalisation |
| `susieR` | Latest CRAN | Multi-causal fine-mapping |
| `ieugwasr` | Latest GitHub | LD clumping; proxy SNPs |
| `MungeSumstats` | Latest Bioconductor | rsID re-annotation; build harmonisation |
| `PRScs` | GitHub `getian107/PRScs` | Polygenic score computation (§6.10) |
| `ldsc` | GitHub `bulik/ldsc` | LD score regression (§8.5) |
| PLINK | 1.9 | Local LD clumping fallback |

A `renv.lock` file is committed at instrument-extraction start to lock all package versions. A `seed.txt` file at project root records all random seeds used (analysis seed: `set.seed(81)`; bootstrap seed: `set.seed(8101)`; permutation seed: `set.seed(8102)`).

All code is published on GitHub at manuscript submission. A Dockerfile is provided for the full analytic environment, with image versions pinned. A reproducibility capsule (Docker image + OSF data deposit + GitHub code release) is the standard reproducibility commitment.

---

## 36. Audit Trail and Pre-Registered Analysis Sequence

The following analytic execution sequence is pre-registered and committed to the project repository as `R/00_pipeline.R`:

1. Ingest and SHA256-log all input sumstats.
2. Compute FinnGen R12 manifest case-count histogram; freeze phecode filter.
3. LDSC heritability + intercept QC on exposure and outcomes (§8.5).
4. rsID re-annotation and harmonisation.
5. Extract instruments per the APOE-inclusive panel and APOE-excluded panel; compute F-statistics. Tier 1 (strict) primary; Tier 2 (relaxed) sensitivity (§6.1).
6. **Decision gate (§6.4):** if APOE-excluded F-statistic or LDSC h² gate fails, activate Pivot A; Tier 3 (PGS-MR) becomes primary.
7. Compute amyloid-PET PGS for Tier 3 sensitivity / fallback (§6.10).
8. Harmonise instruments × phecodes.
9. Run primary IVW per (instrument-panel × phecode).
10. Run sensitivity suite per FDR-significant hit, including non-linear MR + fractional-polynomial dose-response (§10.5).
11. Run APOE three-tier analysis + ε4 MVMR + ε4 dose-response stratified MR (§11).
12. Run reverse MR for AD-spectrum phecodes (§12).
13. Run mediation MVMR through AD / cognitive function (§13.1–§13.4).
14. Run tau-PET cascade two-step MR (§13.5).
15. Run amyloid-direct MVMR (β_amyloid|AD) per phecode (§14).
16. Run pre-symptomatic vs symptomatic decomposition (§14.4).
17. Run cross-ancestry analysis: directional-consistency primary + MR-MEGA secondary + BBJ + Pan-UKB AFR (§15).
18. Run sex-stratified sub-analysis (§16).
19. Run colocalisation at top non-APOE loci per hit; coloc::sensitivity() prior-sweep for intermediate-H4 loci (§17.1–§17.4.1).
20. Run multi-ancestry SuSiE fine-mapping for cross-ancestry-consistent hits (§17.7).
21. Compute LDSC rg cross-tabulation (§18).
22. Run equivalence testing (TOST) for safety-relevant phecodes (§22).
23. Compute power atlas (§25).
24. Generate manuscript figures and tables (§30).
25. Assemble Amyloid-PET Causal Atlas (parquet + static report; §38.5).
26. At manuscript submission: query FAERS / DAEN / EudraVigilance for post-marketing triangulation against FDR-significant hits (§23.4).

Every script logs its execution timestamp, git SHA at execution, input file SHA256s, output file SHA256s, and full R session info (`sessionInfo()`) to `data/audit/`. This audit trail is committed alongside the OSF data deposit at manuscript submission.

Deviations from the execution sequence (e.g., re-running step 9 after step 10) are logged but not automatically problematic. What IS problematic, and what would be logged as an amendment, is any deviation from the analytic specification — not from execution order.

---

## 37. Code Review and Pre-Analysis Audit Checklist

Pre-register that the analysis code under `R/` will undergo a self-audit checklist review before execution on the FDR-discovery step. The checklist:

1. Random seeds explicitly set in every script that uses stochastic methods (MR-PRESSO permutation, bootstrap, CAUSE MCMC, PRS-CS sampling).
2. Per-script `sessionInfo()` written to audit log at execution.
3. Per-script input/output file SHA256 logged.
4. No interactive runs on the primary analysis (everything via `Rscript` or batched executor).
5. Function-level unit tests written for critical functions (instrument extraction, harmonisation, IVW, TOST, GWAS-by-subtraction) before being executed on production data.
6. A dry run of every script on synthetic data (random uniform p-values; simulated allele frequencies) completed without errors before the production run.
7. Git working tree must be clean (no uncommitted changes) at every script-execution timestamp.

The checklist completion is logged in `data/audit/code_review_audit.md` at the time of primary analysis.

---

## 38. Data Management Plan

### 38.1 Raw data

All raw GWAS sumstats hosted at their original public sources. Per-source SHA256 hashes and access dates logged in the data-snapshot companion document.

### 38.2 Derived data

All processed parquets (instrument tables, harmonised outcome tables, MR results, sensitivity outputs, coloc outputs, power atlas, equivalence tests, GWAS-by-subtraction tables, LDSC rg matrix) deposited on OSF as `data/processed/*.parquet` at manuscript submission.

### 38.3 Code

All R scripts under `R/` deposited on GitHub at submission with a Zenodo DOI for permanent archive. Tagged release at submission.

### 38.4 Retention and access

OSF deposit + Zenodo deposit + GitHub repo all maintained ≥ 5 years post-publication.

### 38.5 Amyloid-PET Causal Atlas (pre-registered deliverable)

The full per-phecode results are pre-registered as a standalone community-resource deliverable: the **Amyloid-PET Causal Atlas**. This deliverable is binding regardless of the primary-paper outcome (Pivot A activation does not omit the Atlas).

#### 38.5.1 Atlas contents

A queryable parquet file with one row per (phecode × instrument-panel) combination, containing every statistic reported in the §9.2 schema PLUS:

- MVMR amyloid-direct coefficient β_amyloid|AD (from §14).
- MVMR amyloid-direct conditional on tau coefficient β_amyloid|tau,AD (from §13.5).
- Pre-symptomatic vs symptomatic decomposition coefficients (from §14.4).
- APOE ε4-stratified IVW estimates (from §11.4).
- Equivalence-test (TOST) verdict where applicable (from §22).
- LDSC genetic correlation (rg) with the phecode (from §18).
- Colocalisation H4 per top non-APOE locus (from §17).
- Power-tier classification (from §25.2).
- Cross-ancestry directional-consistency verdict (from §15.3).
- Sex-stratified IVW estimates and Cochran's Q (from §16).
- Non-linear MR fractional-polynomial shape verdict (from §10.5).
- Centiloid-rescaled effect and lecanemab-equivalent magnitude (from §23).

#### 38.5.2 Static report

In addition to the parquet, a static HTML/PDF report is generated with per-phecode pages, each summarising:
- The MR estimate plus sensitivity-method forest plot.
- The decomposition decision (amyloid-direct vs AD-mediated vs tau-mediated).
- The clinical-interpretation paragraph if the phecode is in the §21 clinical-interest panel.
- The power-tier and "underpowered; cannot conclude" flag where applicable.

The static report is generated programmatically from the parquet using a templated Quarto / R Markdown pipeline. No live server is required.

#### 38.5.3 Deposit and DOI

Both the parquet and the static report are deposited at:
- OSF (under the project node, alongside this pre-registration).
- Zenodo (with an independent DOI for citation as a standalone resource).

The Atlas DOI is reported in the manuscript Abstract as a separate citable artefact.

#### 38.5.4 Novelty statement

No public Amyloid-PET Causal Atlas exists. The pre-registered commitment shifts the deliverable from "one paper" to "one paper + one community resource", and the Atlas is publishable independently of whether the primary-paper inference produces a positive headline finding. Pivot A activation (weak-instrument result) does NOT omit the Atlas — under Pivot A the Atlas IS the primary deliverable.

---

## 39. Amendment Log Policy

This registration is **frozen upon OSF deposition**. Deviations during analysis are documented in the amendment-log companion document with:

- Date of deviation.
- Section of the registration affected.
- Specific change and rationale.
- Whether the deviation was identified **before or after viewing FDR-significance results**.

Deviations identified **after** viewing results are clearly labelled as **post-hoc** in the manuscript. Deviations identified **before** viewing results (e.g., a package version change, a software bug fix, a methodological correction) are labelled as **pre-hoc protocol clarifications**.

---

**Part K — Manifest and sign-off**

---

## 40. Pre-Registration Deposit Manifest

The following files are deposited together at OSF as the frozen pre-registration:

| File | Role |
|---|---|
| `registration.md` | This document — master pre-registration |
| `ethics_determination.md` | HREC-exemption determination |
| `exposure_choice_prespec.md` | Exposure-selection rationale and contingency rules |
| `phecode_filter_prespec.csv` | Frozen filtered phecode list (n=1,574; case-count ≥ 500) |
| `power_atlas_prespec.md` | Power-calculation methodology |
| `data_snapshot_log.md` | SHA256 hashes + access dates for input GWAS files |
| `amendment_log.md` | Live during analysis; populated with deviations |

The deposit is uploaded to OSF using the Open-Ended Registration v3 schema, with the `summary` field populated with the abstract from §1 of this document.

---

## 41. Investigator Sign-Off

**Investigator:** Hayden Farquhar MBBS MPHTM
**ORCID:** 0009-0002-6226-440X
**Date:** 2026-05-19
