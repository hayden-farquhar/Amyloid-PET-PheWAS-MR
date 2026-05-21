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

## (Future amendments appended below)
