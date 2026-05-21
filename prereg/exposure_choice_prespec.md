# Exposure-Choice Pre-Specification

**Project:** Amyloid-PET PheWAS-MR vs FinnGen R12
**Investigator:** Hayden Farquhar MBBS MPHTM (ORCID 0009-0002-6226-440X)
**Date:** 2026-05-19
**Status:** FROZEN at OSF deposition

This document pre-specifies the choice of exposure GWAS for the Mendelian randomization analysis and the contingency rules that govern any switch between candidate exposure datasets. It exists to prevent post-hoc exposure-selection bias.

---

## 1. Primary exposure (locked)

**Ali et al. 2023, *Acta Neuropath Comm* 11:68** — multi-ethnic amyloid-PET GWAS (N=13,409).
DOI: 10.1186/s40478-023-01563-4.
Sumstats access: WashU NeuroGenomics + Informatics Center public folder, accessed 2026-05-19.
SHA256s logged in `osf/data_snapshot_log.md`.

**Genome build:** GRCh38 (native; no liftover required).
**Effect-allele coding:** EFFECT_ALLELE = alternate (ALT); harmonised consistently across instrument extraction.
**Sample composition:** 14 cohorts; ancestry composition NHW 86%, ASN 11%, AFR 3%.

The Ali 2023 file inventory includes 7 sub-cohort sumstats:
- Multi-ethnic meta-analysis (N=13,409) — **primary**
- NHW arm (N=11,556) — used in H5 cross-ancestry analysis
- ASN arm (N=1,494) — used in H5; backup EAS replication if Kim 2025 unavailable
- AFR arm (N=359) — used in H5; small-N ancestry anchor
- Female arm (N=5,195) — used in H6 sex-stratified analysis
- Male arm (N=4,625) — used in H6 sex-stratified analysis
- README.txt — file format documentation

## 2. Contingent additional primary exposure (Kim 2025)

**Kim JP, Won H-H, et al. 2025, *Nat Commun* 16:3133** — cross-ancestry amyloid-PET GWAS (N=15,701 = 11,816 EUR + 3,885 EAS Korean).
DOI: 10.1038/s41467-025-57751-4.
Sumstats: GWAS Catalog accession GCST90483382 (**deposit registered, sumstats file "Not available" as of 2026-05-19**).

### Contingency rule (locked)

Corresponding-author contact (Won H-H, SKKU SAIHST; Seo S-W) initiated 2026-05-19 requesting cross-ancestry meta-analysis summary statistics.

**If Kim 2025 sumstats are received by Phase 3 closure (Week 10):**
- Kim 2025 is added as a co-primary exposure.
- Identical instrument-extraction protocol applied (§6 of `osf/registration.md`).
- Per-phecode IVW reported from both Ali 2023 and Kim 2025 instruments.
- MR-MEGA trans-ancestry meta-regression integrates Ali 2023 (multi-ethnic + 3 ancestry arms) + Kim 2025 (EUR+EAS combined) as five exposure panels.

**If Kim 2025 sumstats are NOT received by Phase 3 closure:**
- Ali 2023 ASN arm (N=1,494) serves as the EAS-ancestry instrument source for H5 cross-ancestry analysis.
- The manuscript explicitly notes that Kim 2025 would have provided substantially greater EAS power if obtained.
- No retroactive switch to Kim 2025 is permitted after Phase 3 closure — this commitment prevents post-hoc exposure-selection bias if Kim 2025 sumstats arrive late.

### Trigger for amendment log

The arrival or non-arrival of Kim 2025 sumstats by Phase 3 closure is logged in `osf/amendment_log.md` with the date and decision.

## 3. Replication exposure (locked)

**Jansen IE et al. 2022, EADB CSF Aβ42 GWAS.**
Role: orthogonal exposure for replication of top non-APOE hits (§5.3 of `osf/registration.md`).
Used in Phase 5 only; not used for primary discovery.

## 4. Secondary sensitivity exposure (locked)

**Nho K et al. 2024, ADNI tau-PET GWAS.**
Role: downstream-pathway sensitivity test for AD-spectrum and psychiatric phecode hits, since tau-PET lies downstream of amyloid in the AD biological pathway.
Used in Phase 5 only.

## 5. Excluded exposure candidates and rationale

The following amyloid-PET-related GWASs were considered and explicitly excluded as primary exposures:

| Excluded source | Reason |
|---|---|
| Apostolova et al. 2018 PiB-PET GWAS | Smaller N; tracer-restricted (PiB only); superseded by Ali 2023 |
| Korean Aβ-PET GWAS (Kim 2021 *Alzheimers Res Ther*) | Subsumed by Kim 2025 (Kim 2025 is the cross-ancestry extension) |
| Single-cohort sub-meta-analyses (ADNI-only, DIAN-only) | Not pre-registered as primary; used only as sensitivity if available |

## 6. Locked decision diagram

```
Phase 0 (Week 1-2):    Primary exposure = Ali 2023 multi-ethnic (FROZEN)
                       ├─ Ancestry arms (NHW, ASN, AFR) freeze with primary
                       └─ Sex arms (F, M) freeze with primary

Phase 0-3 (Weeks 1-10): Won lab contact in flight; Kim 2025 contingent

Phase 3 closure (Week 10):
   ├─ Kim 2025 received → add as co-primary; MR-MEGA across 5 panels
   └─ Kim 2025 NOT received → Ali 2023 ASN arm = EAS instrument source
                              No retroactive switch permitted

Phase 5 (Week 12-14):  Replication = Jansen 2022 CSF Aβ42 (LOCKED)
                       Secondary sensitivity = Nho 2024 tau-PET (LOCKED)
```

## 7. End

This document is signed off by the investigator on the date below and frozen at OSF deposit.

**Investigator:** Hayden Farquhar MBBS MPHTM
**Date:** 2026-05-19
