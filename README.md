# Amyloid-PET PheWAS-MR — public reproducibility package

**Manuscript title:** Genetically-predicted amyloid-PET burden localises Alzheimer's disease causation to the chr1q32.2 CR1 region: a pre-registered phenome-wide Mendelian randomization study with cascade triangulation and cross-ancestry replication

**Author:** Hayden Farquhar, MBBS MPHTM (independent researcher; Finley, New South Wales, Australia)
**ORCID:** [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)
**Contact:** hayden.farquhar@icloud.com

**Pre-registration:** OSF DOI [10.17605/OSF.IO/HVEDJ](https://doi.org/10.17605/OSF.IO/HVEDJ) (frozen 2026-05-19; 4 amendments publicly versioned at the OSF project node)

**Zenodo DOI:** (to be assigned at deposit minting; will be inserted here)

**License:** MIT for code (R + Python + bash) | CC-BY 4.0 for derived data products (Causal Atlas, instrument tables, MR result parquets, figures) — see `LICENSE`

**Reporting standards:** STROBE-MR (Skrivankova et al. 2021); PRISMA-MR where applicable

---

## What this deposit contains

This is the reproducibility package for the first pre-registered phenome-wide Mendelian randomization (PheWAS-MR) of amyloid-PET cortical deposition against 1,574 FinnGen R12 disease endpoints. The headline finding is that the conservative wide-APOE-excluded instrument panel localises the AD causal signal to the chr1q32.2 CR1 region (PP.H4 = 0.983 with brain-cortical CR1 expression), with a pre-registered negative-control diagnostic that bounds the strength of the CR1-as-mediator claim.

The deposit includes:

| Directory | Contents | Purpose |
|---|---|---|
| `prereg/` | Locked 41-section pre-registration, 4 amendment logs, ethics determination, exposure-choice pre-spec, power-atlas pre-spec, phecode-filter pre-spec CSV, data-snapshot log | Audit trail of the locked protocol; every analytical decision is referenced to a §-numbered registration section |
| `code/R/` | 14 R scripts (00-16) implementing Phase 1-5 + atlas build + manuscript figures | Re-runnable analytical pipeline |
| `code/scripts/` | 4 Python + 2 shell utilities (tabix sweep, IV-list builder, PRS-CS input formatter, PRS-CS driver, wide-PGS lifter) | Auxiliary code |
| `instruments/` | 16-SNP unified instrument list, Phase 1 verdict (Pivot A activation rationale), 5 per-panel PLINK clumping outputs, wide-PGS IV lists | Reproducible IV construction |
| `results/` | Amyloid-PET Causal Atlas (parquet + TSV mirror + interactive HTML browser), top-hits parquet, per-stage result parquets (locus MR, CSF cascade, coloc, sensitivity MR, PGS-MR narrow + wide, cross-ancestry per-IV + per-panel, Steiger), and 5 manuscript figures (PDF + PNG) | All evidence layers + Atlas browser |
| `supplementary/` | Supplementary Table S1 (per-panel SNP inventory) | Submission-ready supplementary material |
| `docs/` | Software environment specification, data-provenance manifest with upstream SHA256 hashes | Reproducibility metadata |

**Total deposit size:** ~2.5 MB (60 files). Well within Zenodo single-record limits (50 GB; 100 MB per file).

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

**Wall-clock estimate**: Phase 1 + Phase 2 + Phases 3-5 = ~2-4 hours on a modern laptop. PRS-CS MCMC for the polygenic-score layer is the longest single step (~3-6 hours single-threaded on CPU). Total end-to-end ~6-10 hours assuming all upstream data are local.

---

## How to use the Atlas without re-running anything

The `results/atlas_browser.html` is a single-file deployment (~5-10 MB) that uses DataTables.js to render the 7,870-row Causal Atlas as a sortable, filterable, exportable browser. Open it in any modern web browser. No web server required, no JavaScript dependencies beyond what is bundled inline.

For programmatic access, load `results/amyloid_pet_causal_atlas.parquet` in Python with `pandas.read_parquet()` or in R with `arrow::read_parquet()`. The 46-column schema is documented in `code/R/09_build_causal_atlas.R`.

For human-readable inspection, `results/amyloid_pet_causal_atlas.tsv` is a tab-separated mirror (5.1 MB).

---

## Citation

If you use this Atlas or the analytical pipeline in derivative work, please cite both the manuscript and this deposit:

> Farquhar H. Genetically-predicted amyloid-PET burden localises Alzheimer's disease causation to the chr1q32.2 CR1 region: a pre-registered phenome-wide Mendelian randomization study with cascade triangulation and cross-ancestry replication. *[Journal]* (year). DOI: [manuscript DOI]
>
> Farquhar H. Amyloid-PET PheWAS-MR — public reproducibility package (Version 1.0) [Data set]. Zenodo. DOI: [Zenodo DOI to be assigned]

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
