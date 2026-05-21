# Deposit MANIFEST (human-readable)

Total: 89 files, ~2.5 MB. SHA256 hashes are in `MANIFEST.sha256` for machine-readable verification.

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

## `code/scripts/` — Auxiliary scripts

| File | Purpose |
|---|---|
| `build_unified_iv_list.py` | Builds the 16-SNP unified IV list across the 5 panels |
| `build_wide_pgs_noapoe.py` | Lifts top-500 PRS-CS-weighted non-APOE SNPs from GRCh37 → GRCh38 |
| `build_wide_pgs_noapoe.R` | R-side companion to the wide-PGS lifter |
| `make_prscs_sumstats.py` | Converts Ali 2023 GRCh38 sumstats to PRS-CS input format (15.6M → 10.4M SNPs) |
| `run_prscs_amyloid.sh` | PRS-CS MCMC driver (phi=1e-2, seed=81, default 1000 iter / 500 burnin) |
| `sweep_finngen_r12.sh` | Hardened tabix sweep over FinnGen R12 sumstats (~300× faster than full download); includes macOS-specific portability traps |

## `instruments/` — Instrument-variable construction outputs

| File | Purpose |
|---|---|
| `unified_iv_list.tsv` | 16 SNPs × panel-membership (the master IV table; source of Supp Table S1) |
| `PHASE1_VERDICT.md` | Phase 1 verdict document: per-panel n_IV, mean F, and Pivot A activation rationale |
| `ali2023_significant_annotated.tsv` | 313 pre-clump SNPs (p < 1e-6) after liftover + rsID annotation |
| `pgs_wide_iv_list.tsv` | Full top-500 PRS-CS-weighted SNP list (wide-PGS source) |
| `pgs_wide_iv_list_noAPOE.tsv` | Top-50 non-APOE PRS-CS-weighted SNPs (wide-PGS arm) |
| `per_panel/tier1_apoe_inclusive*` | Tier 1 APOE-inclusive (n=6 IVs): .clumped, .log, .nosex, _index_snps.tsv, .assoc |
| `per_panel/tier1_apoe_narrow_excl*` | Tier 1 narrow-excluded (n=5): same 5 files |
| `per_panel/tier1_apoe_wide_excl*` | Tier 1 wide-excluded (n=3; Pivot A primary): same 5 files |
| `per_panel/tier2_apoe_narrow_excl*` | Tier 2 narrow-excluded (n=13): same 5 files |
| `per_panel/tier2_apoe_wide_excl*` | Tier 2 wide-excluded (n=11): same 5 files |

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

## `supplementary/`

| File | Purpose |
|---|---|
| `supplementary_table_S1.md` | Per-panel SNP inventory (16 SNPs × full annotation × panel membership) |

## `docs/`

| File | Purpose |
|---|---|
| `software_environment.txt` | R / Python / system tool versions with macOS-portability notes |
| `data_provenance.md` | Upstream-source URLs, accession IDs, retrieval dates, SHA256 verifications |
