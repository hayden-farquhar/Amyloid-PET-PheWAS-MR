# Deposit MANIFEST (human-readable)

Total: ~108 files, ~12 MB (chr1q32.2 Brain_Cortex full-pairs extract is ~9 MB; the rest is small). SHA256 hashes are in `MANIFEST.sha256` for machine-readable verification.

## Top-level

| File | Purpose |
|---|---|
| `README.md` | Public-facing description of the deposit; how to reproduce; how to cite |
| `LICENSE` | MIT (code) + CC-BY 4.0 (data) dual licence |
| `CITATION.cff` | Citation File Format metadata for Zenodo / GitHub ingestion |
| `MANIFEST.md` | This file |
| `MANIFEST.sha256` | SHA256 hashes for every file in the deposit |

## `prereg/` — Locked pre-registration + amendments

| File | Purpose |
|---|---|
| `registration.md` | The locked 41-section pre-registration document (frozen 2026-05-19) |
| `amendment_log.md` | 4 amendments (all pre-results, all publicly versioned) |
| `data_snapshot_log.md` | Snapshot log of all upstream data acquisitions with dates and hashes |
| `ethics_determination.md` | HREC-exempt determination (independent-researcher use of public summary-level data) |
| `exposure_choice_prespec.md` | Pre-specified rationale for using Ali 2023 multi-ethnic primary |
| `phecode_filter_prespec.csv` | Pre-specified filter list of 1,574 FinnGen R12 phecodes (case-count ≥ 500) |
| `power_atlas_prespec.md` | Pre-specified power calculation across the 1,574-phecode atlas |

## `code/R/` — R analytical pipeline

| File | Purpose |
|---|---|
| `00_hash_raw.R` | SHA256 audit of all upstream sumstats; builds the audit log |
| `05_phase2_locus_mr.R` | Phase 2: per-phecode × per-panel locus-level MR sweep (5 panels × 4 estimators) |
| `06_pgs_mr.R` | PGS-MR (narrow + wide); Pivot A primary inferential layer |
| `07_sensitivity_mr.R` | Phase 3: MR-PRESSO (1,000 bootstraps) + MR-RAPS |
| `08_coloc.R` | Phase 4: coloc.abf against GTEx v8 brain + eQTLGen |
| `09_build_causal_atlas.R` | Master Atlas assembly: joins all evidence layers; 7,870 × 46 |
| `10_atlas_figures.R` | Atlas-visualisation helper figures |
| `11_csf_cascade_mr.R` | Amendment 2: CSF biomarker cascade against Timsina 2026 |
| `12_build_atlas_browser.R` | Single-file HTML browser using DataTables.js |
| `13_replication_ali_asn.R` | Amendment 3: within-study cross-ancestry (Ali ASN N=1,494) |
| `13b_replication_kim_korean.R` | Amendment 4 Component A: between-study cross-ancestry (Kim Korean N=3,387) |
| `14_steiger_directionality.R` | Steiger reverse-causation guard per (panel × phecode × IV) |
| `15_kim_korean_secondary_primary.R` | Amendment 4 Component B: Kim Korean as secondary primary discovery arm |
| `16_manuscript_figures.R` | Generates Fig 2 (CR1 locus zoom), Fig 3 (cross-ancestry concordance), Fig 4 (forest plot) |
| `17_coloc_susie_chr1q322.R` | Per-gene SuSiE fine-mapping at chr1q32.2 on Brain_Cortex full-pairs eQTL data (CR1, CR1L, CD46) + amyloid GWAS lookup |
| `18_figure_chr1q322_regional.R` | chr1q32.2 regional plot: amyloid GWAS Manhattan track + per-gene SuSiE credible sets |
| `19_coloc_susie_chr1q322_eur_formal.R` | Formal amyloid-side coloc.susie attempt with SuSiE-RSS kriging diagnostic + L-progression + outlier removal (documents the multi-SNP-tight-LD limitation; transparency artefact) |
| `19b_coloc_abf_nhw_check.R` | Cross-cohort coloc.abf cross-check using Ali NHW EUR-only sumstats; replicates PP.H4 = 0.977 at CR1 |
| `20_adni_cascade_pgs.R` | ADNI tau-PET PGS cascade biological-validation arm (OSF Amendment 5; post-hoc with model frozen pre-results). Reads from a user-provided local ADNI checkout (DUA-restricted; not redistributed). Computes the pre-registered primary endpoint (META_TEMPORAL_SUVR ~ amyloid_PGS + covariates) + six pre-specified sensitivity analyses (tau-PVC, cross-modality amyloid centiloid, APOE-stratified, EUR-only, longitudinal LMM). |

## `code/scripts/` — Auxiliary scripts

| File | Purpose |
|---|---|
| `build_unified_iv_list.py` | Builds the 16-SNP unified IV list across the 5 panels |
| `build_wide_pgs_noapoe.py` | Lifts top-500 PRS-CS-weighted non-APOE SNPs from GRCh37 → GRCh38 |
| `build_wide_pgs_noapoe.R` | R-side companion to the wide-PGS lifter |
| `make_prscs_sumstats.py` | Converts Ali 2023 GRCh38 sumstats to PRS-CS input format (15.6M → 10.4M SNPs) |
| `run_prscs_amyloid.sh` | PRS-CS MCMC driver (phi=1e-2, seed=81, default 1000 iter / 500 burnin) |
| `sweep_finngen_r12.sh` | Hardened tabix sweep over FinnGen R12 sumstats (~300× faster than full download); includes macOS-specific portability traps |
| `build_references_bib.py` | CrossRef API-based references.bib builder; 32 entries with DOI resolution; manual-overrides dict for ambiguous lookups |

## `instruments/` — Instrument-variable construction outputs

| File | Purpose |
|---|---|
| `unified_iv_list.tsv` | 16 SNPs × panel-membership (the master IV table; source of Supp Table S1) |
| `phase1_instrument_construction.md` | Phase 1 instrument-construction record: per-panel n_IV, mean F, Pivot A activation rationale |
| `ali2023_significant_annotated.tsv` | 313 pre-clump SNPs (p < 1e-6) after liftover + rsID annotation |
| `pgs_wide_iv_list.tsv` | Full top-500 PRS-CS-weighted SNP list (wide-PGS source) |
| `pgs_wide_iv_list_noAPOE.tsv` | Top-50 non-APOE PRS-CS-weighted SNPs (wide-PGS arm) |
| `per_panel/tier1_apoe_inclusive*` | Tier 1 APOE-inclusive (n=6 IVs): `.clumped` (PLINK output) + `_index_snps.tsv` + `.assoc` |
| `per_panel/tier1_apoe_narrow_excl*` | Tier 1 narrow-excluded (n=5): same 3 files |
| `per_panel/tier1_apoe_wide_excl*` | Tier 1 wide-excluded (n=3; Pivot A primary): same 3 files |
| `per_panel/tier2_apoe_narrow_excl*` | Tier 2 narrow-excluded (n=13): same 3 files |
| `per_panel/tier2_apoe_wide_excl*` | Tier 2 wide-excluded (n=11): same 3 files |
| `per_panel/process_logs/*.log` | PLINK execution logs for the 5 per-panel clumping runs (audit value; not analytically informative) |

## `results/` — Derived data products

| File | Purpose |
|---|---|
| `amyloid_pet_causal_atlas.parquet` | **The §38.5 binding deliverable**: 7,870 rows × 46 cols (phecode × panel × all evidence layers joined) |
| `amyloid_pet_causal_atlas.tsv` | TSV mirror of the Atlas parquet (5.1 MB; human-readable) |
| `atlas_browser.html` | Single-file interactive Atlas browser (DataTables.js bundled; no server required) |
| `causal_atlas_top_hits.parquet` | Filtered top-hits view of the Atlas |
| `locus_mr_results.parquet` | Phase 2: 31,480 raw MR estimates (1,574 phecodes × 5 panels × 4 estimators) |
| `csf_cascade_mr_results.parquet` | Amendment 2: CSF cascade against Timsina 2026 (5 panels × 3 biomarkers = 15 estimates) |
| `coloc_results.parquet` | Phase 4: coloc.abf posteriors at lead loci × GTEx brain tissues + eQTLGen |
| `steiger_results.parquet` | Per (panel × phecode × IV) Steiger directionality test outcome |
| `steiger_summary.parquet` | Per (panel × phecode) Steiger pass-fraction summary |
| `sensitivity_mr_results.parquet` | Phase 3: MR-PRESSO + MR-RAPS on screened phecodes |
| `pgs_mr_results_narrow.parquet` | PGS-MR narrow scope (16 unified IVs × PRS-CS shrinkage) |
| `pgs_mr_results_wide.parquet` | PGS-MR wide scope (top-50 non-APOE PRS-CS-weighted SNPs) |
| `kim_korean_locus_mr_results.parquet` | Amendment 4 Component B: full Kim-exposure locus-MR sweep |
| `kim_korean_iv_lookup.parquet` | Amendment 4 Component A: per-IV Kim Korean lookup table |
| `replication_kim_korean_per_iv.parquet` | Per-IV Kim Korean replication outcome (concordance + p) |
| `replication_kim_korean_per_panel.parquet` | Per-panel Kim Korean replication summary |
| `replication_ali_asn_per_iv.parquet` | Per-IV Ali ASN within-study replication |
| `replication_ali_asn_per_panel.parquet` | Per-panel Ali ASN within-study summary |
| `figures/fig1A_volcano_per_panel.{pdf,png}` | Fig 1A: volcano plot per panel |
| `figures/fig1B_concordance_t1wide_t2narrow.{pdf,png}` | Fig 1B: concordance between Tier 1 wide-excl and Tier 2 narrow-excl |
| `figures/fig2_CR1_locuszoom.{pdf,png}` | Fig 2: CR1 locus-zoom |
| `figures/fig3_cross_ancestry_concordance.{pdf,png}` | Fig 3: cross-ancestry concordance |
| `figures/fig4_forest_top_hits.{pdf,png}` | Fig 4: forest plot of top hits |
| `figures/fig5_chr1q322_regional.{pdf,png}` | Fig 5: chr1q32.2 regional plot — amyloid GWAS Manhattan track + per-gene SuSiE credible sets (CR1, CD46; CR1L has no credible set) |
| `susie_chr1q322_brain_cortex_gene_summary.{parquet,tsv}` | Per-gene SuSiE summary at chr1q32.2 (CR1, CR1L, CD46): credible-set lead + amyloid GWAS lookup + distance-from-amyloid-lead |
| `susie_chr1q322_brain_cortex_cs_members.{parquet,tsv}` | Per-credible-set SNP membership (25 rows) with PIPs + amyloid GWAS z/p + eQTL z |
| `coloc_susie_chr1q322_eur_formal.{parquet,tsv}` | Formal amyloid-side coloc.susie outputs (NHW + 1KG EUR LD) — transparency artefact; PP.H4 ≈ 0 / PP.H3 ≈ 1 attributable to known multi-SNP-tight-LD limitation; the replicated coloc.abf is the substantive evidence |
| `coloc_susie_diagnostic_nhw.tsv` | SuSiE-RSS kriging diagnostic (1,758 SNPs; 23 LD-inconsistent outliers flagged) |
| `coloc_susie_intermediates/ali_chr1q322.tsv` | Ali 2023 multi-ethnic amyloid-PET sumstats for chr1q32.2 ± 500 kb (4,938 SNPs); input to R/17 |
| `coloc_susie_intermediates/ali_nhw_chr1q322.tsv` | Ali 2023 NHW sub-meta sumstats for chr1q32.2 ± 500 kb (4,416 SNPs); input to R/19 + R/19b |
| `coloc_susie_intermediates/chr1q32.2_brain_cortex_raw.tsv` | GTEx_V8 Brain_Cortex full-pairs nominal eQTL stats for chr1q32.2 ± 500 kb (59,594 rows = 3,255 SNPs × ~18 genes); fetched via tabix from EBI eQTL Catalogue mirror |
| `coloc_susie_intermediates/chr1q32.2_ld.bim` | PLINK BIM for the chr1q32.2 1KG EUR Phase 3 LD reference subset (2,292 SNPs) |

## `results/adni_cascade/` — ADNI tau-PET PGS cascade arm (Amendment 5; group-level summary only)

| File | Purpose |
|---|---|
| `adni_cascade_lm_results.tsv` | Primary cross-sectional cascade + 5 sensitivity OLS analyses (6 rows × 9 cols: analysis, n, beta, SE, 95% CI, p, R², adj R²) |
| `adni_cascade_lmm_results.tsv` | Longitudinal LMM cascade (tau + amyloid centiloid; 2 rows × 8 cols) |

**DUA-compliance note:** Individual-level ADNI data are NOT redistributed. The `adni_cascade/` subdirectory contains only group-level summary statistics. The full `R/20` analysis script is shareable; the raw ADNI checkout must be obtained from adni.loni.usc.edu under their DUA.

## `supplementary/`

| File | Purpose |
|---|---|
| `supplementary_table_S1.md` | Per-panel SNP inventory (16 SNPs × full annotation × panel membership) |
| `supplementary_table_S2.md` | Per-gene SuSiE credible-set membership at chr1q32.2 in GTEx_V8 Brain_Cortex (CR1 + CD46 credible-set members with PIPs + amyloid GWAS z/p + eQTL z; CR1L noted as having no credible set) |
| `supplementary_table_S3.md` | ADNI cascade arm (Amendment 5): cohort details + per-analysis full coefficients (mirrors the Table S3 in the manuscript submission supplementary file) |

## `docs/`

| File | Purpose |
|---|---|
| `software_environment.txt` | R / Python / system tool versions with macOS-portability notes |
| `data_provenance.md` | Upstream-source URLs, accession IDs, retrieval dates, SHA256 verifications |
