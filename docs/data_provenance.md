# Data provenance manifest

This document records the canonical sources, accession identifiers, retrieval dates, and SHA256 hashes of every upstream data resource used in this analysis. The deposit does NOT redistribute these upstream resources; re-download via the URLs below and verify against the hashes given here. The derived data products in `results/`, `instruments/`, and `prereg/` of this deposit are the outputs of the pipeline operating on these upstream resources.

## Upstream summary-level GWAS sources

### Ali et al. 2023 — Multi-ethnic amyloid-PET cortical SUVR GWAS

- **Primary publication:** Ali M, Archer DB, Gorijala P, Western D, Timsina J, Fernández-Manjón B, Western D et al. *Acta Neuropathologica Communications* (2023). DOI: [10.1186/s40478-023-01563-4](https://doi.org/10.1186/s40478-023-01563-4)
- **N:** 13,409 multi-ethnic (11,556 NHW European + 1,494 East Asian + 359 African American; 5,195 female + 4,625 male sub-meta arms also available)
- **Genome build:** GRCh38 native
- **Acquisition:** WashU NeuroGenomics and Informatics Center, [https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l](https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l)
- **Retrieved (primary acquisition):** 2026-05-19
- **SHA256 hashes** (from `data/raw/ali_2023/SHA256SUMS` in source project):
  - `795a4f78970f4513d90902f867d5f5825645d3c8c6d04c9fd07d874e4e3b5f4e` — `AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt` (primary exposure; ~960 MB)
  - `626bab2b560b364cd2fd612a8a75f76665778d6a199f4d6000709cef3b810c5e` — `AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt` (Ali ASN, within-study cross-ancestry reference)
  - `f489b6f38594ba983ae88fe7ab24a36257beaf9580d265858f0e278fb4b89788` — `AmyloidPET_GWAS_AFR_joint_cohorts8_n359_hg38_joint_WashU.txt`
  - `d20cb72a369e429ca407a516c0e59a72f3da3dcb7a886126fabb4fdc57772330` — `AmyloidPET_GWAS_Female_multiEthnic_metaAnalysis_cohorts9_n5195_hg38_WashU.txt`
  - `df1a0766dac854bce23632609eae281e1f73fd77b957247992b0cbd469b23078` — `AmyloidPET_GWAS_Male_multiEthnic_metaAnalysis_cohorts9_n4625_hg38_WashU.txt`
- **Note on local copy:** Only the multi-ethnic primary arm (~960 MB) is retained on the development workstation. The 5 secondary arms (NHW, ASN, AFR, Female, Male; ~3.3 GB total) were deleted post-analysis for disk-space reasons; the hashes above suffice for re-download verification.

### FinnGen R12 — Outcome phecode summary statistics

- **Source:** [FinnGen results portal](https://www.finngen.fi/en/access_results)
- **N:** 500,348 Finnish-ancestry individuals; 1,574 phecodes after case-count ≥ 500 filter (full release contains ~2,400 endpoints)
- **Genome build:** GRCh38 (FinnGen-imputed against the global LD reference panel)
- **Acquisition pattern:** bgzip-compressed, tabix-indexed sumstats on the FinnGen public Google Cloud bucket. Per-IV slices pulled via remote byte-range queries (`tabix -h <URL> -R <regions.bed>`) — see `code/scripts/sweep_finngen_r12.sh` for the hardened pattern (~300× faster than full-file download).
- **Retrieved (FinnGen sweeps):** 2026-05-19 through 2026-05-21
- **Slices retained:** Per-phecode per-IV-position outcome statistics for 16 unified IVs × 1,574 phecodes = ~25,000 sliced rows total. Original development-host storage is on Google Drive at `_shared_mirror/finngen_R12/sliced/`; not redistributed in this deposit (FinnGen consortium release-policy terms).
- **Verification:** FinnGen-published per-phecode checksums on the FinnGen results portal. Re-download verification against the FinnGen-bgzip-and-tbi files is sufficient.

### Timsina et al. 2026 — CSF biomarker GWAS (Aβ42, t-tau, p-tau181)

- **Primary publication:** Timsina J et al. *Nat Commun* (2026). PMID 42014397 (verify when bib is built — this is a future-dated PMID; confirm the actual accession at submission time)
- **N:** 18,948 EUR for each of three biomarkers (Aβ42, total tau, phosphorylated-tau 181)
- **GWAS Catalog accessions:** [GCST90726396](https://www.ebi.ac.uk/gwas/studies/GCST90726396) (Aβ42), [GCST90726397](https://www.ebi.ac.uk/gwas/studies/GCST90726397) (t-tau), [GCST90726398](https://www.ebi.ac.uk/gwas/studies/GCST90726398) (p-tau181)
- **Genome build:** GRCh38 (confirm at re-download)
- **Retrieved (cascade acquisition):** 2026-05-20
- **Verification:** EBI GWAS Catalog publishes checksums alongside each summary statistic file.

### Kim et al. 2025 — Korean amyloid-PET GWAS

- **Primary publication:** Kim JP, Won H-H et al. *Nat Commun* (2025), 16:3133. DOI: [10.1038/s41467-025-57751-4](https://doi.org/10.1038/s41467-025-57751-4)
- **N:** 3,387 East Asian individuals (Korean cohort; the publication also reports a cross-ancestry meta-analysis at larger N that is not separately deposited)
- **GWAS Catalog accession:** [GCST90483382](https://www.ebi.ac.uk/gwas/studies/GCST90483382)
- **Genome build:** GRCh38 (Kim uses TOPMed-R2 imputation, which differs from Ali 2023's reference panel; this is the source of the 7-of-16 IV unavailability that affects the cross-ancestry analysis — see Results)
- **Acquisition:** EBI GWAS Catalog FTP per correspondence with the Won lab (Sungkyunkwan University) confirming public accessibility of the Korean cohort component (see acknowledgements)
- **Retrieved:** 2026-05-21

## Upstream reference and eQTL sources

### GTEx v8 — Brain cis-eQTL data (Phase 4 colocalisation + chr1q32.2 per-gene SuSiE)

- **Tissues used:** 13 brain regions: Amygdala, Anterior cingulate cortex BA24, Caudate basal ganglia, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex BA9, Hippocampus, Hypothalamus, Nucleus accumbens basal ganglia, Putamen basal ganglia, Spinal cord cervical c-1, Substantia nigra
- **Genome build:** GRCh38
- **Source (genome-wide screen — significant-pairs):** [GTEx Portal Google Cloud bucket](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/) (significant-pairs only; FDR<0.05 per tissue; ~154 MB total across 13 tissues). Used for the Phase 4 coloc.abf screen across all lead loci.
- **Source (chr1q32.2 full-pairs for per-gene SuSiE):** [EBI eQTL Catalogue tabix mirror](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cortex.tsv.gz) — bgzip + tabix-indexed; full-pairs nominal eQTL data for the chr1q32.2 ± 500 kb window fetched via tabix HTTP byte-range query (~9 MB region pull rather than the 3.4 GB tissue file). Used for `code/R/17_coloc_susie_chr1q322.R` (per-gene SuSiE on CR1 / CR1L / CD46) and `code/R/19_coloc_susie_chr1q322_eur_formal.R` (formal amyloid-side coloc.susie attempt).
- **Citation:** GTEx Consortium *Science* 2020; EBI eQTL Catalogue (Kerimov et al. 2021 *Nat Genet*)

### eQTLGen Phase 1 — Whole-blood cis-eQTL significant-pairs

- **N:** 31,684 (Võsa et al. 2021 *Nat Genet*)
- **Source:** [eQTLGen download portal](https://www.eqtlgen.org/)
- **Genome build:** GRCh37 (used as comparison reference to GTEx v8 brain tissues)

### 1000 Genomes Phase 3 EUR LD reference

- **Source:** [International Genome Sample Resource](https://www.internationalgenome.org/)
- **Use in this analysis:** PLINK 1.9 clumping (r²<0.001, kb<10,000) for IV construction in all 5 panels
- **Local symlink:** Reused from a companion analysis to avoid duplicating the ~5 GB reference; this deposit does not redistribute the reference panel.
- **Verification:** 1000G-published checksums

### PRS-CS 1KG EUR LD reference (Broad ldblk_1kg_eur)

- **Source:** Broad-hosted archive linked from [Ge et al. 2019 PRS-CS GitHub](https://github.com/getian107/PRScs)
- **Size:** ~4.3 GB unpacked (22 HDF5 LD blocks)
- **Use in this analysis:** PRS-CS MCMC for Bayesian shrinkage of amyloid-PET GWAS effect sizes → polygenic-score weighted MR
- **Citation:** Ge T, Chen CY, Ni Y, Feng YA, Smoller JW (2019) *Nat Commun*

## Pre-registration and amendment trail

- **OSF DOI:** [10.17605/OSF.IO/HVEDJ](https://doi.org/10.17605/OSF.IO/HVEDJ)
- **Frozen date:** 2026-05-19
- **Amendments:** 4 (all pre-results, all publicly versioned at the OSF project node):
  - Amendment 1 (wide-APOE-exclusion sensitivity arm; 2026-05-19)
  - Amendment 2 (Timsina 2026 CSF biomarker GWAS as primary cascade exposure; 2026-05-20T04:04:19 UTC)
  - Amendment 3 (Ali ASN substitution for BBJ cross-ancestry; 2026-05-20T04:09:13 UTC)
  - Amendment 4 (Kim Korean as primary cross-ancestry + secondary primary discovery arm; 2026-05-20/21)
- **Full text:** `prereg/registration.md` (locked) + `prereg/amendment_log.md` (4 amendments)
