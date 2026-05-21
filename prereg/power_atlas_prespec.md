# Power Atlas Pre-Specification

**Project:** Amyloid-PET PheWAS-MR vs FinnGen R12
**Investigator:** Hayden Farquhar MBBS MPHTM (ORCID 0009-0002-6226-440X)
**Date:** 2026-05-19
**Status:** FROZEN at OSF deposition

This document pre-specifies the methodology used to compute the per-phecode minimum-detectable-OR power atlas. The atlas is the primary deliverable of Pivot A (weak-instrument atlas) and a mandatory supplementary deliverable of the primary analysis.

---

## 1. Method

Analytic two-sample MR power calculation per **Burgess (2014) *Int J Epidemiol* 43(3):922–929** (DOI 10.1093/ije/dyu005). The formula expresses power as a function of:
- α (significance threshold)
- N_outcome (total participants in outcome GWAS)
- K (proportion of cases in outcome GWAS)
- R²_exposure (variance in the exposure explained by the instruments)
- True causal OR (the effect-size to detect)

The minimum-detectable-OR (MDOR) at a fixed power target is the value of the true causal OR that yields 80% (or other specified) power at the specified α.

## 2. Per-(phecode × instrument-panel) computation

For each combination of FinnGen R12 phecode (n=1,574 from `osf/phecode_filter_prespec.csv`) and instrument panel (APOE-inclusive; APOE-excluded), compute:

| Field | Computation | Source |
|---|---|---|
| `N_outcome` | num_cases + num_controls | FinnGen R12 manifest |
| `K` | num_cases / N_outcome | FinnGen R12 manifest |
| `R2_exposure` | Sum across instruments of 2 × β² × MAF × (1−MAF) / Var(exposure); per-instrument F-statistic available in `data/interim/instruments_*.csv` | Ali 2023 (instrument panel post-clumping, F≥10) |
| `mean_F` | Mean F-statistic across instruments | Ali 2023 |
| `n_IV` | Count of instruments | Ali 2023 |
| `MDOR_alpha_005` | Burgess 2014 formula at α=0.05, power=80% | Computed |
| `MDOR_alpha_bonferroni` | At α = 0.05 / (1,574 × 2) = 1.59 × 10⁻⁵, power=80% | Computed |
| `MDOR_alpha_fdr_chapter` | At chapter-specific BH-FDR threshold (computed post-hoc; approximated pre-hoc) | Computed |
| `required_N_for_OR085` | Outcome N required to detect OR=0.85 at 80% power | Computed |
| `power_status` | adequately_powered (MDOR ≤ 1.20) / moderately_powered (1.20 < MDOR ≤ 1.50) / underpowered (MDOR > 1.50) | Computed |

## 3. Reported deliverables

- `data/processed/power_atlas.parquet` — full atlas table.
- `results/figures/power_atlas_heatmap.png` — heatmap of MDOR per (phecode chapter × instrument panel).
- `results/tables/power_classification_summary.csv` — counts of phecodes by power status × chapter.

## 4. Pivot A (weak-instrument atlas) trigger and deliverables

If APOE-excluded instrument panel mean F < 15 OR n_IV < 5 (the gate from §3.2 / §21.1 of `osf/registration.md`), Pivot A is activated and the power atlas becomes the manuscript's primary deliverable.

In that case, the manuscript additionally reports:
- A "future amyloid-PET GWAS sample-size recommendation" curve: for each adequately-powered OR (1.10, 1.15, 1.20, 1.25, 1.30), what outcome-side N (case + control combined) is required to achieve 80% power?
- A "future exposure-side GWAS recommendation": what exposure-side N would be required to bring R²_exposure (APOE-excluded) to the level that enables MDOR ≤ 1.20 against the 25th-percentile FinnGen phecode case count?

These two curves frame the recommendations for future amyloid-PET GWAS efforts and guide the manuscript's Discussion under Pivot A.

## 5. Bootstrap-based empirical power (supplementary)

For phecodes that pass FDR q < 0.05 in the primary analysis, supplement the analytic Burgess 2014 power calculation with a bootstrap-based empirical power estimate:

1. Bootstrap-resample the instrument panel (1,000 resamples; with replacement).
2. Per resample, compute IVW estimate.
3. Empirical power at α = 1.59 × 10⁻⁵ = fraction of resamples returning p < 1.59 × 10⁻⁵.
4. Report empirical power alongside analytic power; concordance check (analytic vs empirical) reported.

Bootstrap implementation seed: `set.seed(8101)`.

## 6. Honest interpretation guidance

For every reported FDR-negative phecode, the power tier is cross-referenced (§20.4 of `osf/registration.md`). A phecode with FDR-negative AND underpowered status (MDOR > 1.50) is explicitly framed in the manuscript as "underpowered; cannot conclude", not as a "negative finding".

This commitment is binding: the manuscript Limitations section reports the proportion of phecodes that fall in each power-status × verdict category, so reviewers can independently assess the power-conditional interpretation.
