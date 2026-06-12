# Amyloid-PET PheWAS-MR — public reproducibility package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20334183.svg)](https://doi.org/10.5281/zenodo.20334183)
[![Pre-registration](https://img.shields.io/badge/OSF-10.17605%2FOSF.IO%2FHVEDJ-blue)](https://doi.org/10.17605/OSF.IO/HVEDJ)
[![Licence: MIT + CC-BY 4.0](https://img.shields.io/badge/licence-MIT%20%2B%20CC--BY%204.0-lightgrey)](LICENSE)

**Manuscript title:** Pre-registered phenome-wide Mendelian randomization of amyloid-PET: instrument limits and a colocalisation-input artefact

**Author:** Hayden Farquhar, MBBS MPHTM (independent researcher; Finley, New South Wales, Australia)
**ORCID:** [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)
**Contact:** hayden.farquhar@icloud.com

**Pre-registration:** OSF DOI [10.17605/OSF.IO/HVEDJ](https://doi.org/10.17605/OSF.IO/HVEDJ) (frozen 2026-05-19; 4 amendments publicly versioned at the OSF project node)

**Zenodo conceptDOI (all versions):** [10.5281/zenodo.20334183](https://doi.org/10.5281/zenodo.20334183) — stable identifier that always resolves to the latest version of this deposit.

**Latest version (v1.1.0, 2026-06-12):** corrects a colocalisation-input artefact and supersedes the colocalisation numbers of earlier releases (see "What this deposit contains" and the v1.1.0 changelog below). The conceptDOI above always resolves to the latest version; the version-specific DOI is assigned at release.

**License:** MIT for code (R + Python + bash) | CC-BY 4.0 for derived data products (Causal Atlas, instrument tables, MR result parquets, figures) — see `LICENSE`

**Reporting standards:** STROBE-MR (Skrivankova et al. 2021); PRISMA-MR where applicable

---

## What this deposit contains

This is the reproducibility package for a pre-registered phenome-wide Mendelian randomization (PheWAS-MR) of amyloid-PET cortical deposition against 1,574 FinnGen R12 disease endpoints. The scan is instrument-limited: outside APOE the conservative panel resolves to three instruments (effectively the CR1 variant rs6656401), and a per-phecode statistical-power floor (R/21) shows the scan detects only large effects (median minimum detectable odds ratio ≈ 2.8 at the within-panel Bonferroni threshold). The locus that survives is the previously-known CR1 region, where genetically-predicted amyloid, AD risk, and brain-cortical CR1 expression share a causal variant for late-onset AD-spectrum dementia (Brain_Cortex coloc.abf PP.H4 = 0.98) — a shared-variant result that is method- and prior-dependent (coloc.susie disagrees) and that does not establish amyloid mediation (reverse-direction MR, R/25, is itself significant).

**v1.1.0 changelog (2026-06-12) — supersedes earlier colocalisation numbers.** Earlier releases ran coloc.abf on GTEx significant-pairs eQTL data and summarised it as a maximum posterior across tissues; at the tight-LD CR1 locus this inflated false-positive colocalisation (all 16 coloc-tested phecodes, including a forearm-fracture negative control, returned ≈ 0.98). The colocalisation layer is re-run with full-pairs eQTL data in the pre-specified Brain_Cortex tissue (R/24, R/29): the colocalisation is then specific to the late-onset AD-spectrum phecodes and the negative control passes. A pre-specified 24-phenotype negative-control panel (R/23, R/27) and a simulation on the real CR1 LD matrix (R/28) bound the artefact, which is a known colocalisation failure mode (Giambartolomei 2014; Wallace 2021). The earlier "evidence-symmetry / tissue-of-action" framing is retired. New analyses also add the power floor (R/21), colocalisation-prior robustness (R/22), and reverse-direction MR (R/25). Do not cite the pre-v1.1.0 colocalisation numbers.

The deposit includes:

| Directory | Contents | Purpose |
|---|---|---|
| `prereg/` | Locked 41-section pre-registration, 4 amendment logs, ethics determination, exposure-choice pre-spec, power-atlas pre-spec, phecode-filter pre-spec CSV, data-snapshot log | Audit trail of the locked protocol; every analytical decision is referenced to a §-numbered registration section |
| `code/R/` | 18 R scripts (00–19, including R/13b) implementing Phase 1–5 + atlas build + manuscript figures + the chr1q32.2 per-gene SuSiE decomposition (R/17), the regional plot (R/18), and the formal coloc.susie attempt with kriging diagnostic (R/19) + the coloc.abf cross-cohort cross-check (R/19b) | Re-runnable analytical pipeline |
| `code/scripts/` | Python + shell utilities (tabix sweep, IV-list builder, PRS-CS input formatter, PRS-CS driver, wide-PGS lifter, CrossRef-based references.bib builder) | Auxiliary code |
| `instruments/` | 16-SNP unified instrument list, Phase 1 verdict (Pivot A activation rationale), 5 per-panel PLINK clumping outputs, wide-PGS IV lists | Reproducible IV construction |
| `results/` | Amyloid-PET Causal Atlas (parquet + TSV mirror + interactive HTML browser), top-hits parquet, per-stage result parquets (locus MR, CSF cascade, coloc, sensitivity MR, PGS-MR narrow + wide, cross-ancestry per-IV + per-panel, Steiger), per-gene SuSiE outputs at chr1q32.2 (per-gene credible-set membership + amyloid GWAS lookup; formal coloc.susie outputs + kriging diagnostic), regional figure for the chr1q32.2 cascade, and 5 + 1 manuscript figures (PDF + PNG) | All evidence layers + Atlas browser |
| `supplementary/` | Supplementary Table S1 (per-panel SNP inventory) + Supplementary Table S2 (per-gene SuSiE credible-set membership at chr1q32.2) | Submission-ready supplementary material |
| `docs/` | Software environment specification, data-provenance manifest with upstream SHA256 hashes | Reproducibility metadata |

**Total deposit size:** ~12 MB (~110 files; the chr1q32.2 Brain_Cortex full-pairs eQTL extract is ~9 MB). Well within Zenodo single-record limits (50 GB; 100 MB per file).

---

## What this deposit does NOT contain (and where to get it)

This package is the *derived* reproducibility layer. The upstream summary-level GWAS data and reference panels are not redistributed because (a) they are large (tens of GB to hundreds of GB), (b) they are already deposited at canonical sources, and (c) they carry their own licences and data-use agreements that must be honoured at the point of download. The data-provenance manifest (`docs/data_provenance.md`) lists exact source URLs, accession IDs, and SHA256 hashes for verification at re-download.

| Required input | Source | SHA256 verification |
|---|---|---|
| Ali 2023 multi-ethnic amyloid-PET GWAS (N=13,409) | [WashU Box](https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l) | See `docs/data_provenance.md` |
| FinnGen R12 phecode summary statistics (N=500,348; 1,574 endpoints) | [FinnGen results portal](https://www.finngen.fi/en/access_results); tabix-indexed bgzip on the FinnGen public bucket | bgzip + .tbi indexed; verify via FinnGen-published checksums |
| Timsina 2026 CSF biomarker GWAS (Aβ42, t-tau, p-tau181; N=18,948 EUR) | GWAS Catalog [GCST90726396](https://www.ebi.ac.uk/gwas/studies/GCST90726396), [GCST90726397](https://www.ebi.ac.uk/gwas/studies/GCST90726397), [GCST90726398](https://www.ebi.ac.uk/gwas/studies/GCST90726398) | EBI-published checksums |
| Kim 2025 Korean amyloid-PET GWAS (N=3,387 EAS) | GWAS Catalog [GCST90483382](https://www.ebi.ac.uk/gwas/studies/GCST90483382) | EBI-published checksums |
| GTEx v8 brain cis-eQTL significant-pairs (13 brain tissues) | [GTEx Portal Google Cloud bucket](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/) | GTEx-published checksums |
| eQTLGen Phase 1 whole-blood cis-eQTLs (Võsa 2021) | [eQTLGen download portal](https://www.eqtlgen.org/) | eQTLGen-published checksums |
| 1000 Genomes Phase 3 EUR LD reference | [International Genome Sample Resource](https://www.internationalgenome.org/) | 1KG-published checksums |
| PRS-CS 1KG EUR LD reference | [Ge et al. 2019 PRS-CS GitHub](https://github.com/getian107/PRScs) (Broad ldblk_1kg_eur archive) | Broad-published checksums |

---

## How to reproduce the analysis from scratch

1. **Pull the upstream data** per the table above. Place into the directory structure expected by the pipeline (see `code/R/00_hash_raw.R` for the canonical layout). The total disk footprint is ~10-15 GB once GTEx v8 brain cis-eQTLs are pulled.

2. **Install software** per `docs/software_environment.txt`. The critical pins are R 4.5.2 with `TwoSampleMR` 0.7.6 / `MRPRESSO` 1.0 / `coloc` 5.2.3 / `arrow` 24.0.0 / `data.table` 1.18.2.1; Python 3.14 with `pyliftover`, `numpy` 2.x (with two patches to upstream PRS-CS — see `code/scripts/` for the patched files); PLINK 1.9; htslib 1.23.1 (tabix).

   **Project-root resolution:** All R and Python scripts use a portable pattern to find the project root. Either (a) export `AMYLOID_PHEWASMR_ROOT=/path/to/your/working/copy` in your shell before running, or (b) ensure your R working directory is set so that `here::here()` resolves to the project root, or (c) simply run each script from the project root (the scripts fall back to `getwd()` / `here::here()` when the env var is not set). The R scripts depend on the `here` package; install with `install.packages("here")` if not already present.

3. **Reproduce Phase 1 instrument construction**: run `code/scripts/build_unified_iv_list.py` then `code/R/00_hash_raw.R`. This produces `instruments/unified_iv_list.tsv` and the per-panel PLINK clumping outputs.

4. **Reproduce Phase 2 locus-level MR**: run `code/R/05_phase2_locus_mr.R`. This pulls FinnGen R12 outcome statistics per IV-position via `tabix` byte-range queries (see `code/scripts/sweep_finngen_r12.sh` for the hardened tabix pattern; ~300× faster than full-file downloads).

5. **Reproduce Phase 3-5**: run `code/R/07_sensitivity_mr.R` (Phase 3), `code/R/08_coloc.R` (Phase 4 colocalisation), `code/R/11_csf_cascade_mr.R` (CSF cascade, Amendment 2), `code/R/13_replication_ali_asn.R` + `code/R/13b_replication_kim_korean.R` + `code/R/15_kim_korean_secondary_primary.R` (Amendments 3 + 4 cross-ancestry).

6. **Reproduce PGS-MR** (Pivot A primary inferential layer): run `code/scripts/run_prscs_amyloid.sh` (uses patched PRS-CS source for numpy 2.x compatibility), then `code/R/06_pgs_mr.R` (narrow + wide modes).

7. **Reproduce Steiger directionality**: run `code/R/14_steiger_directionality.R`.

8. **Build the Causal Atlas + browser**: run `code/R/09_build_causal_atlas.R` (produces `results/amyloid_pet_causal_atlas.parquet`) and `code/R/12_build_atlas_browser.R` (produces `results/atlas_browser.html`).

9. **Reproduce manuscript figures**: run `code/R/16_manuscript_figures.R`.

10. **Reproduce the chr1q32.2 per-gene SuSiE decomposition + cross-cohort coloc.abf cross-check** (post-hoc; addresses the multi-causal architecture that coloc.abf alone cannot resolve): run `code/R/17_coloc_susie_chr1q322.R` (per-gene SuSiE on full-pairs Brain_Cortex eQTL data fetched byte-range from the EBI eQTL Catalogue mirror; ~9 MB region pull), `code/R/18_figure_chr1q322_regional.R` (regional plot), `code/R/19_coloc_susie_chr1q322_eur_formal.R` (formal amyloid-side coloc.susie attempt with SuSiE-RSS kriging diagnostic + L-progression + outlier removal; documents the multi-SNP-tight-LD limitation), and `code/R/19b_coloc_abf_nhw_check.R` (cross-cohort coloc.abf cross-check using Ali NHW EUR-only sumstats; replicates PP.H4 = 0.977 at CR1, ruling out EAS-arm-driven inflation).

**Wall-clock estimate**: Phase 1 + Phase 2 + Phases 3-5 = ~2-4 hours on a modern laptop. PRS-CS MCMC for the polygenic-score layer is the longest single step (~3-6 hours single-threaded on CPU). Total end-to-end ~6-10 hours assuming all upstream data are local.

---

## How to use the Atlas without re-running anything

The `results/atlas_browser.html` is a single-file deployment (~5-10 MB) that uses DataTables.js to render the 7,870-row Causal Atlas as a sortable, filterable, exportable browser. Open it in any modern web browser. No web server required, no JavaScript dependencies beyond what is bundled inline.

For programmatic access, load `results/amyloid_pet_causal_atlas.parquet` in Python with `pandas.read_parquet()` or in R with `arrow::read_parquet()`. The 46-column schema is documented in `code/R/09_build_causal_atlas.R`.

For human-readable inspection, `results/amyloid_pet_causal_atlas.tsv` is a tab-separated mirror (5.1 MB).

---

## Citation

If you use this Atlas or the analytical pipeline in derivative work, please cite both the manuscript and this deposit:

> Farquhar H. Pre-registered phenome-wide Mendelian randomization of amyloid-PET: instrument limits and a colocalisation-input artefact. *[Journal]* (year). DOI: [manuscript DOI]
>
> Farquhar H. Amyloid-PET PheWAS-MR — public reproducibility package (Version 1.0.1) [Data set]. Zenodo. DOI: [10.5281/zenodo.20334254](https://doi.org/10.5281/zenodo.20334254). Latest version: [10.5281/zenodo.20334183](https://doi.org/10.5281/zenodo.20334183).

Pre-registration:

> Farquhar H. Phenome-wide Mendelian randomization of genetically-proxied cortical amyloid-β deposition against 1,574 FinnGen R12 disease endpoints (pre-registration). OSF (2026). DOI: [10.17605/OSF.IO/HVEDJ](https://doi.org/10.17605/OSF.IO/HVEDJ)

A machine-readable `CITATION.cff` is included at the root of this deposit for direct ingestion by Zenodo and GitHub.

---

## Funding, ethics, and conflicts of interest

- **Funding:** No external funding. This work was conducted as part of an independent research portfolio.
- **Ethics:** HREC-exempt; uses only publicly-available summary-level data. See `prereg/ethics_determination.md`.
- **Competing interests:** None declared.
- **AI use disclosure:** Anthropic's Claude (claude-opus-4-7 with 1M-context) was used during analytical pipeline development, manuscript drafting, and reproducibility-package assembly. All scientific judgments, results, and editorial decisions are the author's. The use of Claude is consistent with this deposit's transparency commitments; logs and prompts can be made available on request.

---

## Contact

For questions about the analysis, the data, or the Atlas, email hayden.farquhar@icloud.com. The OSF project at [10.17605/OSF.IO/HVEDJ](https://doi.org/10.17605/OSF.IO/HVEDJ) is monitored for follow-up.
