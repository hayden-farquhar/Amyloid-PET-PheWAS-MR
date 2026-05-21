# Data Snapshot Log

**Project:** Amyloid-PET PheWAS-MR vs FinnGen R12
**Investigator:** Hayden Farquhar MBBS MPHTM (ORCID 0009-0002-6226-440X)
**Pre-registration date:** 2026-05-19

This log records SHA256 hashes, file sizes, source URLs, and access dates for all GWAS summary-statistics files used in this study. Hashes are the GNU coreutils `sha256sum` format. Per-source download logs are maintained as `data/raw/<source>/download_log.tsv` and SHA256SUMS files in the same directory.

---

## Ali et al. 2023 — Multi-ethnic amyloid-PET GWAS

**Publication:** Ali M, et al. *Acta Neuropath Comm* (2023) 11:68. DOI 10.1186/s40478-023-01563-4
**Host:** Washington University NeuroGenomics + Informatics Center (public Box folder share)
**Source URL:** https://wustl.box.com/s/ahiy99e3rkmlrsln28yy4dck8i9hq64l
**Access date:** 2026-05-19
**Genome build:** GRCh38 (native; no liftover applied)

## SHA256 Hashes

```
f489b6f38594ba983ae88fe7ab24a36257beaf9580d265858f0e278fb4b89788  AmyloidPET_GWAS_AFR_joint_cohorts8_n359_hg38_joint_WashU.txt
626bab2b560b364cd2fd612a8a75f76665778d6a199f4d6000709cef3b810c5e  AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt
d20cb72a369e429ca407a516c0e59a72f3da3dcb7a886126fabb4fdc57772330  AmyloidPET_GWAS_Female_multiEthnic_metaAnalysis_cohorts9_n5195_hg38_WashU.txt
df1a0766dac854bce23632609eae281e1f73fd77b957247992b0cbd469b23078  AmyloidPET_GWAS_Male_multiEthnic_metaAnalysis_cohorts9_n4625_hg38_WashU.txt
795a4f78970f4513d90902f867d5f5825645d3c8c6d04c9fd07d874e4e3b5f4e  AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt
5bcc261c41e5269eb1dbd2b9fbf650b906a9e3fcbe1c7776f275900bfd7c12eb  AmyloidPET_GWAS_NHW_metaAnalysis_cohorts13_n11556_hg38_WashU.txt
deb72e2ca2d72eee1288ef8d5eddd344a2cd6fb6fba5588c4beeeea5f93c2a71  README.txt
```

## Per-file metadata

| filename | bytes | use |
|---|---|---|
| AmyloidPET_GWAS_AFR_joint_cohorts8_n359_hg38_joint_WashU.txt | 683795584 | African arm (N=359) — H5 AFR instrument source |
| AmyloidPET_GWAS_ASN_metaAnalysis_cohorts10_n1494_hg38_WashU.txt | 748096352 | Asian arm (N=1,494) — H5 EAS instrument source |
| AmyloidPET_GWAS_Female_multiEthnic_metaAnalysis_cohorts9_n5195_hg38_WashU.txt | 616449940 | Female arm (N=5,195) — H6 sex-stratified |
| AmyloidPET_GWAS_Male_multiEthnic_metaAnalysis_cohorts9_n4625_hg38_WashU.txt | 591929837 | Male arm (N=4,625) — H6 sex-stratified |
| AmyloidPET_GWAS_multiEthnic_metaAnalysis_cohorts14_n13409_hg38_WashU.txt | 1006881564 | Primary multi-ethnic exposure (N=13,409) |
| AmyloidPET_GWAS_NHW_metaAnalysis_cohorts13_n11556_hg38_WashU.txt | 876386514 | Non-Hispanic White / EUR arm (N=11,556) — H5 cross-ancestry |
| README.txt | 1396 | Dataset documentation |

---

## (Other sources — populated when files are downloaded)

The following sources will be logged here when files are downloaded during Phases 1–5:

- **Kim et al. 2025** *Nat Commun* (DOI 10.1038/s41467-025-57751-4) — GWAS Catalog GCST90483382 (deposit pending; corresponding-author request sent 2026-05-19)
- **Bellenguez et al. 2022** *Nat Genet* AD GWAS (DOI 10.1038/s41588-022-01024-z)
- **Davies et al. 2018** *Nat Commun* general cognitive function GWAS (DOI 10.1038/s41467-018-04362-x)
- **Jansen et al. 2022** EADB CSF Aβ42 GWAS
- **Nho et al. 2024** ADNI tau-PET GWAS
- **FinnGen R12** phecode summary statistics (per `osf/phecode_filter_prespec.csv` — 1,574 endpoints)
- **Biobank Japan** matched phenotype summary statistics (per-hit only)
- **Pan-UKB AFR** matched phenotype summary statistics (per-hit only)
- **GTEx v8** brain + whole-blood eQTL
- **eQTLGen** consortium whole-blood eQTL
- **1000 Genomes Phase 3** LD reference panels (EUR, EAS, AFR)

Per-source SHA256SUMS and download_log.tsv files are populated by `R/00_hash_raw.R` at the time of download.
