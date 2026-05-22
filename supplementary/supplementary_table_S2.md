# Supplementary Table S3 — Per-gene SuSiE credible-set members at chr1q32.2 (GTEx_V8 Brain_Cortex)

Per-gene SuSiE fine-mapping of cis-eQTL signals at chr1q32.2 in GTEx_V8 Brain_Cortex (n=1,758 SNPs after rsID match to 1KG EUR Phase 3 LD reference). Credible-set members are listed with their PIP (posterior inclusion probability), the amyloid-PET GWAS z-score and p-value at the same SNP (from Ali 2023 multi-ethnic, sign-aligned to the eQTL alt allele), and the per-gene eQTL z-score. CR1L is included as a tested gene that returned no credible set (omitted from the table). Analysis script: `R/17_coloc_susie_chr1q322.R`. Source data: `data/processed/susie_chr1q322_brain_cortex_cs_members.parquet`.

| Gene | Gene ID | CS | rsID | Position (GRCh38) | PIP | Amyloid z | Amyloid p | eQTL z |
|---|---|---|---|---|---|---|---|---|
| CD46 | ENSG00000117335 | 1 | rs2724384 | chr1:207,756,858 | 0.187 | -1.26 | 2.09e-01 | +13.91 |
| CD46 | ENSG00000117335 | 1 | rs2724374 | chr1:207,767,846 | 0.143 | -1.42 | 1.57e-01 | +13.81 |
| CD46 | ENSG00000117335 | 1 | rs2796265 | chr1:207,747,443 | 0.143 | -1.34 | 1.81e-01 | +13.81 |
| CD46 | ENSG00000117335 | 1 | rs2466572 | chr1:207,761,142 | 0.132 | -1.21 | 2.26e-01 | +13.81 |
| CD46 | ENSG00000117335 | 1 | rs2761434 | chr1:207,744,154 | 0.132 | -1.24 | 2.16e-01 | +13.81 |
| CD46 | ENSG00000117335 | 1 | rs2724360 | chr1:207,769,813 | 0.132 | -1.43 | 1.53e-01 | +13.81 |
| CD46 | ENSG00000117335 | 1 | rs2761437 | chr1:207,749,736 | 0.129 | -1.19 | 2.34e-01 | +13.81 |
| CR1 | ENSG00000203710 | 1 | rs679515 | chr1:207,577,223 | 0.101 | -6.33 | 2.47e-10 | -6.53 |
| CR1 | ENSG00000203710 | 1 | rs12037841 | chr1:207,510,847 | 0.096 | -5.76 | 8.56e-09 | -6.52 |
| CR1 | ENSG00000203710 | 1 | rs7515905 | chr1:207,564,732 | 0.096 | -6.03 | 1.68e-09 | -6.52 |
| CR1 | ENSG00000203710 | 1 | rs1752684 | chr1:207,573,951 | 0.096 | -2.65 | 7.98e-03 | -6.52 |
| CR1 | ENSG00000203710 | 1 | rs2093760 | chr1:207,613,483 | 0.067 | -2.74 | 6.05e-03 | -6.45 |
| CR1 | ENSG00000203710 | 1 | rs3818361 | chr1:207,611,623 | 0.067 | -2.71 | 6.64e-03 | -6.45 |
| CR1 | ENSG00000203710 | 1 | rs6701713 | chr1:207,612,944 | 0.067 | -6.30 | 3.03e-10 | -6.45 |
| CR1 | ENSG00000203710 | 1 | rs6661489 | chr1:207,524,699 | 0.062 | -2.50 | 1.24e-02 | -6.42 |
| CR1 | ENSG00000203710 | 1 | rs4844610 | chr1:207,629,207 | 0.055 | -6.20 | 5.70e-10 | -6.39 |
| CR1 | ENSG00000203710 | 1 | rs10863417 | chr1:207,623,552 | 0.053 | -0.71 | 4.76e-01 | -6.39 |
| CR1 | ENSG00000203710 | 1 | rs4562624 | chr1:207,512,620 | 0.042 | -6.26 | 3.74e-10 | -6.33 |
| CR1 | ENSG00000203710 | 1 | rs6656401 | chr1:207,518,704 | 0.042 | -6.34 | 2.27e-10 | -6.33 |
| CR1 | ENSG00000203710 | 1 | rs4266886 | chr1:207,512,441 | 0.031 | -2.87 | 4.16e-03 | -6.26 |
| CR1 | ENSG00000203710 | 1 | rs4844600 | chr1:207,505,962 | 0.022 | -5.54 | 2.96e-08 | -6.17 |
| CR1 | ENSG00000203710 | 1 | rs2093761 | chr1:207,613,197 | 0.019 | -2.54 | 1.11e-02 | -6.15 |
| CR1 | ENSG00000203710 | 1 | rs10779335 | chr1:207,625,349 | 0.018 | -6.06 | 1.38e-09 | -6.14 |
| CR1 | ENSG00000203710 | 1 | rs10863420 | chr1:207,626,529 | 0.018 | -2.59 | 9.55e-03 | -6.14 |
| CR1 | ENSG00000203710 | 1 | rs2296160 | chr1:207,621,975 | 0.018 | -2.60 | 9.27e-03 | -6.14 |

**CR1L (ENSG00000203709):** no credible set returned (no significant cis-eQTL signal in Brain_Cortex at the chr1q32.2 ± 500 kb window).

**Note on credible-set lead SNPs:**
- CR1 CS 1 lead by PIP: **rs679515** (chr1:207,577,223; PIP 0.101; amyloid z = −6.33, p = 2.5×10⁻¹⁰; 58.5 kb from the amyloid GWAS lead rs6656401). rs679515 is the second-strongest amyloid signal in the chr1q32.2 1 Mb window after rs6656401 itself.
- CD46 CS 1 lead by PIP: **rs2724384** (chr1:207,756,858; PIP 0.187; amyloid z = −1.26, p = 0.21; 238.2 kb from rs6656401). The credible set sits well outside the amyloid signal envelope (top-5 amyloid SNPs all within chr1:207,512,620 – 207,629,207).
