# Supplementary Table S1 — Per-panel instrument-variable inventory

**Source:** `data/interim/instruments/unified_iv_list.tsv` (Phase 1, locked 2026-05-19; sourced from Ali et al. 2023 multi-ethnic amyloid-PET cortical SUVR GWAS, N=13,409, GRCh38 native build).

**Construction:** PLINK 1.9 clumping against 1000 Genomes Phase 3 EUR reference (r²<0.001, kb<10,000), per-SNP F-statistic = β²/SE², F-statistic threshold ≥10 enforced at panel level. Panel definitions per Methods §"Per-panel definitions". GRCh38 positions are the native Ali 2023 build; hg19 positions are derived via the UCSC `hg38ToHg19.over.chain.gz` lift via pyliftover and are reported for cross-reference with legacy GRCh37 resources (PRS-CS LD reference, 1KG legacy panels).

## Full instrument list (16 SNPs across 5 panels)

| # | rsID | Chr | hg38 pos | hg19 pos | EA | OA | EAF | β (Ali EUR) | SE | p | F | Panel membership |
|---:|:---|:---:|---:|---:|:---:|:---:|---:|---:|---:|---:|---:|:---|
| 1 | **rs6656401** | 1 | 207,518,704 | 207,692,049 | A | G | 0.175 | +0.098 | 0.016 | 2.4×10⁻¹⁰ | 40.2 | T1-incl, T1-narrow, **T1-wide**, T2-narrow, T2-wide — CR1 (canonical AD locus); Pivot A primary instrument |
| 2 | rs4971567 | 2 | 50,897,030 | 51,124,168 | T | G | NA | +1.054 | 0.207 | 3.3×10⁻⁷ | 26.0 | T2-narrow, T2-wide |
| 3 | rs12713280 | 2 | 54,854,239 | 55,081,376 | T | C | 0.253 | −0.060 | 0.012 | 4.5×10⁻⁷ | 25.3 | T2-narrow, T2-wide — chr2 NRXN1 region |
| 4 | rs10050940 | 5 | 15,819,160 | 15,819,269 | T | C | 0.193 | +0.051 | 0.010 | 6.3×10⁻⁷ | 24.6 | T2-narrow, T2-wide — chr5 FBXL7 region |
| 5 | rs7982 | 8 | 27,604,964 | 27,462,481 | A | G | 0.385 | −0.063 | 0.012 | 2.4×10⁻⁷ | 26.6 | T2-narrow, T2-wide |
| 6 | rs138650168 | 8 | 113,833,617 | 114,845,846 | A | G | 0.039 | +0.145 | 0.030 | 8.2×10⁻⁷ | 24.3 | T2-narrow, T2-wide |
| 7 | rs10748157 | 12 | 40,785,080 | 41,178,882 | A | G | NA | +0.744 | 0.144 | 2.2×10⁻⁷ | 26.9 | T2-narrow, T2-wide |
| 8 | **rs117834516** | 14 | 52,848,729 | 53,315,447 | T | C | 0.053 | +0.156 | 0.026 | 2.0×10⁻⁹ | 35.9 | T1-incl, T1-narrow, **T1-wide**, T2-narrow, T2-wide — chr14q21.3 (TXNDC16 / GPR137C region); Pivot A primary instrument |
| 9 | rs146124233 | 16 | 72,988,790 | 73,022,689 | G | A | 0.041 | +0.149 | 0.030 | 9.4×10⁻⁷ | 24.0 | T2-narrow, T2-wide — chr16 ZFHX3 region |
| 10 | **rs12151021** | 19 | 1,050,875 | 1,050,874 | A | G | 0.327 | +0.072 | 0.013 | 9.2×10⁻⁹ | 32.9 | T1-incl, T1-narrow, **T1-wide**, T2-narrow, T2-wide — chr19p13.3 distal (~44 Mb from APOE; outside APOE LD block); Pivot A primary instrument |
| 11 | rs4803748 | 19 | 44,743,791 | 45,247,048 | T | C | 0.397 | −0.092 | 0.012 | 1.7×10⁻¹⁴ | 59.3 | T1-narrow, T2-narrow — APOE LD-spillover SNP (~165 kb upstream of rs429358); inside Amendment 1 wide-APOE exclusion |
| 12 | rs2965169 | 19 | 44,747,899 | 45,251,156 | C | A | 0.403 | −0.101 | 0.012 | 3.6×10⁻¹⁷ | 71.5 | T1-incl only — APOE LD-spillover SNP (~161 kb upstream of rs429358) |
| 13 | rs56394238 | 19 | 44,832,419 | 45,335,676 | G | A | 0.460 | +0.081 | 0.012 | 6.0×10⁻¹² | 47.2 | T1-incl only — APOE LD-spillover SNP (~76 kb upstream of rs429358) |
| 14 | rs6857 | 19 | 44,888,997 | 45,392,254 | T | C | 0.205 | +0.335 | 0.010 | 1.4×10⁻²⁵⁸ | 1,170.0 | T1-narrow, T2-narrow — APOE LD-spillover (~20 kb upstream of rs429358) |
| 15 | rs429358 | 19 | 44,908,684 | 45,411,941 | C | T | 0.188 | +0.350 | 0.009 | 6.2×10⁻³¹¹ | 1,419.6 | T1-incl only — APOE ε4 defining variant |
| 16 | rs142855795 | 22 | 21,838,290 | 22,192,579 | A | G | 0.018 | +0.252 | 0.049 | 3.3×10⁻⁷ | 26.1 | T2-narrow, T2-wide |

**Bolded SNPs (rs6656401, rs117834516, rs12151021)** constitute the Tier-1 wide-APOE-excluded panel — the Pivot A primary instrument set (n=3; min F=32.9; mean F=36.3; max F=40.2).

**EAF = NA** for 3 SNPs (rs4971567, rs10748157, rs10748157 [duplicated entry — verify Ali sumstats]; rs4971567): Ali 2023 sumstats had missing allele-frequency values at these positions. Tier 2 panels carry these SNPs via their effect sizes only; harmonisation handles allele matching against FinnGen R12 outcome statistics via TwoSampleMR::harmonise_data action=2 (palindromic-SNP exclusion at MAF 0.42–0.58 — none of these three SNPs are palindromic at non-extreme MAF, so all 16 IVs survive harmonisation against FinnGen).

## Panel composition summary

| Panel | n_IV | Mean F | Min F | Max F | Member rsIDs |
|:---|---:|---:|---:|---:|:---|
| Tier 1 APOE-inclusive | 6 | 274.6 | 32.9 (rs12151021) | 1,419.6 (rs429358) | rs429358, rs2965169, rs56394238, rs117834516, rs6656401, rs12151021 |
| Tier 1 APOE-narrow excluded | 5 | 267.6 | 32.9 (rs12151021) | 1,170.0 (rs6857) | rs6857, rs4803748, rs117834516, rs6656401, rs12151021 |
| **Tier 1 APOE-wide excluded (Pivot A primary)** | **3** | **36.3** | **32.9 (rs12151021)** | **40.2 (rs6656401)** | **rs6656401, rs117834516, rs12151021** |
| Tier 2 APOE-narrow excluded | 13 | 118.6 | 24.0 (rs146124233) | 1,170.0 (rs6857) | rs4971567, rs12713280, rs10050940, rs7982, rs138650168, rs10748157, rs117834516, rs146124233, rs12151021, rs4803748, rs6857, rs142855795, rs6656401 |
| Tier 2 APOE-wide excluded | 11 | 28.4 | 24.0 (rs146124233) | 40.2 (rs6656401) | rs4971567, rs12713280, rs10050940, rs7982, rs138650168, rs10748157, rs117834516, rs146124233, rs12151021, rs142855795, rs6656401 |

## APOE LD-spillover annotation

The Amendment 1 wide-APOE exclusion boundary (chr19:44,400,000–46,500,000 GRCh38) was introduced after a pre-clump scan revealed an LD-spillover SNP cluster at chr19:44.88 Mb falling outside the §6.4 registered narrow boundary (chr19:44,906,000–45,916,000). The four LD-spillover SNPs are:

| rsID | hg38 pos | Distance from rs429358 (chr19:44,908,684) | In T1-narrow? | In T1-wide? |
|:---|---:|---:|:---:|:---:|
| rs2965169 | 44,747,899 | ~160.8 kb upstream | No | No |
| rs4803748 | 44,743,791 | ~164.9 kb upstream | Yes | No |
| rs56394238 | 44,832,419 | ~76.3 kb upstream | No | No |
| rs6857 | 44,888,997 | ~19.7 kb upstream | Yes | No |

Two of these (rs6857, rs4803748) survive the narrow-boundary exclusion but are dropped by the wide-boundary exclusion. They are the principal motivation for Amendment 1 and they are the SNPs whose effect direction flips between Ali EUR and Kim Korean — empirical confirmation that the narrow boundary admits APOE LD-spillover instruments (see Results §"Negative-control diagnostic and coloc-evidence symmetry" for the discussion).

## Reproducibility

All clumping outputs (`.clumped`, `.log`, `.nosex`) are deposited in the OSF project (DOI 10.17605/OSF.IO/HVEDJ) under `data/interim/instruments/`. The clumping command is recorded in `osf/registration.md` §6.4. SHA256 verification of the source Ali 2023 sumstats is in `data/raw/ali_2023/SHA256SUMS`.

---

*Footer:* Format suitable for submission as Supplementary Table S1 in PDF or Word at submission time. Re-render as Excel or CSV for journals that prefer machine-readable supplementary tables.
