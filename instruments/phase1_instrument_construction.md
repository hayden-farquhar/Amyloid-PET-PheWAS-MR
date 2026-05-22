# Phase 1 Instrument-Extraction Verdict

**Date:** 2026-05-19
**Source:** Ali et al. 2023 *Acta Neuropath Comm* multi-ethnic amyloid-PET GWAS (N=13,409)
**Pre-registration:** OSF DOI 10.17605/OSF.IO/HVEDJ
**Amendment 1:** APOE boundary sensitivity (narrow primary + wide sensitivity)
**Clumping:** PLINK 1.9, 1000 Genomes Phase 3 EUR reference panel (GRCh37 with positional liftover applied to Ali 2023 GRCh38 sumstats); `--clump-r2 0.001 --clump-kb 10000`

## Post-clump instrument count per panel

| Panel | Pre-clump SNPs | Post-clump n_IV | Mean F | Decision |
|---|---|---|---|---|
| Tier 1 APOE-inclusive | 249 | 6 | 274.6 | Anchored by rs429358 (F=1,420) |
| Tier 1 APOE-narrow-excl (§6.4 primary) | 206 | 5 | 267.6 | At the n_IV ≥ 5 boundary |
| **Tier 1 APOE-wide-excl (Amendment 1 sensitivity)** | **12** | **3** | **36.3** | **Pivot A TRIGGER FIRES (n_IV < 5)** |
| Tier 2 APOE-narrow-excl (relaxed p<1e-6) | 260 | 13 | 118.6 | Adequate for MR-RAPS primary |
| Tier 2 APOE-wide-excl (relaxed + wide) | 43 | 11 | 28.4 | Adequate for MR-RAPS primary |

## Tier 1 APOE-wide-excluded index SNPs (3 IVs — the conservative non-APOE panel)

| rsID | chr | bp (GRCh38) | EAF | β | SE | p | F | Likely locus |
|---|---|---|---|---|---|---|---|---|
| rs6656401 | 1 | 207,518,704 | 0.175 | +0.0983 | 0.0155 | 2.40e-10 | 40.2 | **CR1** — canonical AD locus |
| rs117834516 | 14 | 52,848,729 | 0.053 | +0.1557 | 0.0260 | 1.98e-09 | 35.9 | chr14q21.3 — TXNDC16 / GPR137C region |
| rs12151021 | 19 | 1,050,875 | 0.327 | +0.0717 | 0.0125 | 9.25e-09 | 32.9 | chr19p13.3 distal — GAMT / NDUFS7 |

## LD-spillover SNPs that fall outside narrow APOE but inside wide APOE LD-block

These two SNPs are kept by the §6.4 narrow boundary but excluded by Amendment 1 wide boundary. They have extraordinarily strong genetic signal (F = 1,170 and 59) but sit within the broader APOE LD block:

| rsID | chr | bp (GRCh38) | EAF | β | SE | p | F | Distance from rs429358 (chr19:44,908,684) |
|---|---|---|---|---|---|---|---|---|
| rs6857 | 19 | 44,888,997 | 0.205 | +0.3352 | 0.0098 | 1.38e-258 | 1,169.9 | ~20 kb upstream |
| rs4803748 | 19 | 44,743,791 | 0.397 | −0.0924 | 0.0120 | 1.67e-14 | 59.3 | ~165 kb upstream |

These two SNPs are exactly the LD-spillover instruments that Amendment 1 was designed to exclude. Without the wide-boundary sensitivity, the narrow-panel n_IV = 5 would have appeared to pass the Pivot A gate, but two of those five are not independent of APOE.

## Pivot A activation verdict (per Amendment 1 conservative interpretation)

**Pivot A activates** at the Tier 1 stage under the wide-boundary interpretation:
- APOE-wide-excluded Tier 1 n_IV = 3 < 5 (n_IV threshold)
- APOE-wide-excluded Tier 1 mean F = 36.3 > 15 (F-statistic threshold passes)

The conservative interpretation per Amendment 1 ("Pivot A activates if EITHER boundary fails") triggers Pivot A at the n_IV criterion under the wide boundary.

## Operational consequence

Per §27 of the registration, Pivot A activation means:
1. **The Amyloid-PET Causal Atlas (§38.5) is now the primary deliverable.** The Atlas was always pre-registered as binding regardless of pivot status, but under Pivot A it becomes the central artefact of the manuscript.
2. **PGS-MR (Tier 3, §6.10) becomes primary** for non-APOE inference. The PRS-CS polygenic-score instrument from the full Ali 2023 GWAS provides a single composite IV with dramatically higher statistical power than the 3 independent loci.
3. **Tier 2 relaxed (p < 1e-6) with MR-RAPS as primary estimator** (§6.1) becomes the parallel inferential layer at the locus level. Both APOE-narrow (n_IV = 13) and APOE-wide (n_IV = 11) Tier 2 panels are adequately powered.
4. **The boundary-sensitivity contrast itself** becomes a methodological mini-finding for the manuscript — the chr19:44.88 Mb cluster characterisation, and the demonstration that the conventional narrow APOE boundary admits LD-spillover instruments.

## Genuine non-APOE-LD loci confirmed for downstream MR (under wide-boundary interpretation)

1. **CR1** (chr1, rs6656401) — canonical AD risk locus, replicates from Bellenguez 2022.
2. **chr14q21.3** (rs117834516) — less well-characterised for amyloid; warrants colocalisation at this locus.
3. **chr19p13.3 distal** (rs12151021) — distal chr19, well outside the APOE LD block (rs12151021 is at chr19:1.05 Mb, APOE at chr19:44.9 Mb; the distance precludes any LD relationship).

## Files written

- `data/interim/instruments/ali2023_significant_annotated.tsv` — 313 SNPs (p < 1e-6, post-liftover, post-rsID-annotation)
- `data/interim/instruments/tier{1,2}_apoe_{inclusive,narrow_excl,wide_excl}_index_snps.tsv` — clumped instrument tables per panel
- `data/interim/instruments/tier*_apoe_*_clumped.{clumped,log}` — PLINK clumping outputs

All SHA256-hashable for audit trail commitment to OSF post-Phase 5.
