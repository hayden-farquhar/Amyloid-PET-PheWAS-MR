# =============================================================================
# Script 12: Build the self-contained Causal Atlas HTML browser
# =============================================================================
#
# Generates a single-file HTML (results/atlas_browser.html) that any researcher
# can open in a browser to search, filter, and inspect the Amyloid-PET Causal
# Atlas. This is the §38.5 binding deliverable (parquet + static HTML/PDF +
# Zenodo DOI; binding regardless of primary-paper outcome).
#
# The HTML embeds the atlas data as JSON (no external data file needed) and
# uses DataTables.js (loaded via CDN) for the sortable/filterable table.
# Per-row expansion reveals the full per-method estimates + coloc evidence.
#
# Total page size: ~5-10 MB for a 7,870-row atlas. Loads fast on any modern
# browser. Works offline once cached. No server required.
#
# Run as:
#   Rscript R/12_build_atlas_browser.R [atlas_path] [output_html]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
})

# --- Paths -------------------------------------------------------------------
proj_root <- Sys.getenv("AMYLOID_PHEWASMR_ROOT", unset = here::here())
setwd(proj_root)

args <- commandArgs(trailingOnly = TRUE)
atlas_fp <- if (length(args) >= 1) args[1] else
            file.path(proj_root, "data/processed/amyloid_pet_causal_atlas.parquet")
out_html <- if (length(args) >= 2) args[2] else
            file.path(proj_root, "results/atlas_browser.html")

dir.create(dirname(out_html), showWarnings = FALSE, recursive = TRUE)
if (!file.exists(atlas_fp)) stop(sprintf("Atlas parquet not found: %s", atlas_fp))

atlas <- as.data.table(read_parquet(atlas_fp))
cat(sprintf("[%s] Atlas: %d rows × %d cols\n", Sys.time(), nrow(atlas), ncol(atlas)))

# --- Prepare display data ----------------------------------------------------

# Round numeric columns sensibly for display
fmt_p <- function(x) ifelse(is.na(x), NA, signif(x, 3))
fmt_b <- function(x) ifelse(is.na(x), NA, round(x, 4))

display <- atlas[, .(
  phecode = phecode,
  panel = panel,
  phenotype = phenotype,
  category = ifelse(is.na(category), "(uncategorised)", category),
  n_cases = num_cases,
  n_iv = n_iv,
  b_ivw = fmt_b(b_ivw),
  se_ivw = fmt_b(se_ivw),
  p_ivw = fmt_p(pval_ivw),
  fdr_q = fmt_p(ivw_fdr_q_within_panel),
  b_wm = fmt_b(b_weighted_median),
  p_wm = fmt_p(pval_weighted_median),
  b_egger = fmt_b(b_egger),
  p_egger = fmt_p(pval_egger),
  egger_int_p = fmt_p(egger_intercept_pval),
  q_stat = fmt_b(q_stat),
  q_p = fmt_p(q_pval),
  robust = ifelse(is.na(robust_hit), FALSE, robust_hit),
  coloc_pp_h4 = if ("best_pp_h4" %in% names(atlas)) fmt_b(best_pp_h4) else NA_real_,
  coloc_gene = if ("top_coloc_gene" %in% names(atlas)) top_coloc_gene else NA_character_,
  coloc_tissue = if ("top_coloc_tissue" %in% names(atlas)) top_coloc_tissue else NA_character_,
  coloc_strengthened = if ("coloc_strengthened" %in% names(atlas))
    ifelse(is.na(coloc_strengthened), FALSE, coloc_strengthened) else FALSE
)]

# Order by IVW p-value
display <- display[order(p_ivw, na.last = TRUE)]

# Category palette for badges
categories <- sort(unique(display$category))
panels_unique <- sort(unique(display$panel))

# Convert to JSON (compact, NA → null)
data_json <- toJSON(display, na = "null", null = "null", digits = 6)
categories_json <- toJSON(categories)
panels_json <- toJSON(panels_unique)

# Quick stats for the header
n_rows <- nrow(display)
n_phecodes <- uniqueN(display$phecode)
n_robust <- sum(display$robust, na.rm = TRUE)
n_coloc <- sum(display$coloc_strengthened, na.rm = TRUE)

# --- Atlas metadata ----------------------------------------------------------

last_updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
osf_doi <- "10.17605/OSF.IO/HVEDJ"
project_title <- "Amyloid-PET Phenome-Wide Mendelian Randomization Atlas"

# --- Emit HTML --------------------------------------------------------------

html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>', project_title, '</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://cdn.datatables.net/2.0.8/css/dataTables.dataTables.min.css">
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "SF Pro", Helvetica, Arial, sans-serif;
         margin: 1.5em; color: #222; line-height: 1.45; max-width: 1500px; }
  h1 { font-size: 1.6em; margin-bottom: 0.2em; }
  h2 { font-size: 1.1em; margin-top: 1.8em; color: #444; }
  .subtitle { color: #666; margin-bottom: 1.5em; }
  .stats { display: flex; gap: 2em; margin: 1em 0; flex-wrap: wrap; }
  .stat { background: #f0f3f6; padding: 0.7em 1.1em; border-radius: 5px; }
  .stat .big { font-size: 1.5em; font-weight: 600; }
  .stat .lbl { font-size: 0.8em; color: #666; }
  table.dataTable { font-size: 0.85em; width: 100% !important; }
  table.dataTable thead th { font-size: 0.78em; }
  td.dt-control { cursor: pointer; }
  .badge { display: inline-block; padding: 0.1em 0.5em; border-radius: 3px;
           font-size: 0.7em; font-weight: 600; }
  .robust { background: #d4edda; color: #155724; }
  .coloc { background: #d1ecf1; color: #0c5460; }
  .filters { background: #f9fafb; padding: 1em; border-radius: 6px; margin-bottom: 1em; }
  .filters label { margin-right: 1.5em; font-size: 0.9em; }
  .detail { background: #fafbfc; padding: 1em; margin: 0.5em 0; border-left: 3px solid #aaa; }
  .detail h4 { margin: 0.3em 0 0.5em 0; }
  .detail dl { display: grid; grid-template-columns: max-content auto;
               gap: 0.2em 1em; font-size: 0.85em; }
  .detail dt { font-weight: 600; color: #555; }
  .pval-strong { font-weight: 600; color: #a00; }
  .footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #ddd;
            font-size: 0.85em; color: #666; }
  code { background: #f0f3f6; padding: 0.1em 0.3em; border-radius: 3px; font-size: 0.85em; }
</style>
</head>
<body>

<h1>', project_title, '</h1>
<p class="subtitle">Pre-registered phenome-wide Mendelian randomization of amyloid-PET cortical
deposition against 1,574 FinnGen R12 disease endpoints. Pre-registration:
<a href="https://doi.org/', osf_doi, '">OSF DOI ', osf_doi, '</a>.
Last updated: ', last_updated, '.</p>

<div class="stats">
  <div class="stat"><div class="big">', n_rows, '</div><div class="lbl">total rows (phecode × panel)</div></div>
  <div class="stat"><div class="big">', n_phecodes, '</div><div class="lbl">unique phecodes</div></div>
  <div class="stat"><div class="big">', n_robust, '</div><div class="lbl">robust hits</div></div>
  <div class="stat"><div class="big">', n_coloc, '</div><div class="lbl">coloc-strengthened</div></div>
</div>

<div class="filters">
  <label>Panel: <select id="panel-filter"><option value="">All</option></select></label>
  <label>ICD chapter: <select id="cat-filter"><option value="">All</option></select></label>
  <label><input type="checkbox" id="robust-only"> Robust hits only</label>
  <label><input type="checkbox" id="coloc-only"> Coloc-strengthened only</label>
  <label><input type="checkbox" id="sig-only"> p<sub>IVW</sub> < 1e-4 only</label>
</div>

<table id="atlas" class="display compact" style="width:100%">
<thead><tr>
  <th></th>
  <th>Phecode</th>
  <th>Phenotype</th>
  <th>Category</th>
  <th>Panel</th>
  <th>n IVs</th>
  <th>β IVW</th>
  <th>SE IVW</th>
  <th>p IVW</th>
  <th>FDR q</th>
  <th>PP.H4</th>
  <th>Flags</th>
</tr></thead>
</table>

<h2>How to read this table</h2>
<ul>
  <li><strong>Panel:</strong> instrument set used. Tier 1 wide-excluded is the Pivot A primary
      (3 IVs excluding the wide APOE LD block); Tier 2 panels are the relaxed-significance
      MR-RAPS primary; Tier 1 APOE-inclusive is descriptive only (APOE dominates).</li>
  <li><strong>β IVW + p IVW:</strong> inverse-variance-weighted causal estimate of amyloid-PET
      effect on the log-odds of the phecode. β > 0 = amyloid burden increases risk.</li>
  <li><strong>FDR q:</strong> Benjamini-Hochberg corrected within panel across ~1,574 phecodes.</li>
  <li><strong>PP.H4:</strong> colocalisation posterior probability of a shared causal variant
      between amyloid-PET and the best brain eQTL at the lead locus.</li>
  <li><strong>Robust flag:</strong> requires IVW p < 1e-4, weighted-median direction match,
      non-significant Egger intercept, non-significant Cochran\'s Q. See pre-reg §28.</li>
  <li><strong>Click any row</strong> to expand and see all four MR methods + heterogeneity
      diagnostics + colocalisation top hit.</li>
</ul>

<h2>Methodological caveats</h2>
<ul>
  <li>Pivot A is active: Tier 1 wide-excluded has only 3 IVs (CR1, chr14q21.3, chr19p13.3 distal).
      Sparse instruments increase weak-instrument bias risk. The PGS-MR layer (PRS-CS-weighted)
      is the manuscript\'s primary inferential layer for non-APOE inference; the locus-level
      results here are sensitivity-suite outputs.</li>
  <li>The genome-wide GTEx download ships significant-pairs only (full-pairs are 50-200 GB per tissue).
      For the chr1q32.2 CR1 locus we additionally fetched full-pairs nominal eQTL stats from the
      EBI eQTL Catalogue mirror (~9 MB byte-range pull) and ran per-gene SuSiE fine-mapping
      (CR1 / CR1L / CD46) plus cross-cohort coloc.abf replication (PP.H4 = 0.977 NHW EUR-only
      ≈ 0.982 multi-ethnic re-run ≈ 0.983 original). For other lead loci, PP.H4 > 0.5
      hits in the Atlas table remain a coloc <em>screen</em>; full-pairs follow-up at the EBI
      Catalogue is the reproducible next step.</li>
  <li>FinnGen R12 outcomes use Finnish founder-population LD; harmonisation against EUR
      exposure IVs may miss minor proxy-driven effects. See pre-reg §6.8 for details.</li>
  <li>The atlas is locked at the manuscript submission cycle, but the pre-registration\'s
      <a href="https://doi.org/10.17605/OSF.IO/HVEDJ">amendment log</a> remains the authoritative
      record of any protocol deviations.</li>
</ul>

<div class="footer">
  <p>Citation: Farquhar H. (2026). <em>', project_title, '</em>. Open Science Framework.
  <a href="https://doi.org/', osf_doi, '">https://doi.org/', osf_doi, '</a></p>
  <p>Data downloads: <a href="amyloid_pet_causal_atlas.parquet">Parquet</a> |
  <a href="amyloid_pet_causal_atlas.tsv">TSV</a> (place files alongside this HTML).</p>
  <p>Built ', last_updated, ' from ', n_rows, ' atlas rows.</p>
</div>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/2.0.8/js/dataTables.min.js"></script>
<script>
const ATLAS = ', data_json, ';
const CATEGORIES = ', categories_json, ';
const PANELS = ', panels_json, ';

function isSig(x) { return x !== null && x < 1e-4; }
function pvalCell(x) {
  if (x === null) return "";
  const cls = isSig(x) ? "pval-strong" : "";
  return `<span class="${cls}">${x.toExponential(2)}</span>`;
}
function flagCell(row) {
  let parts = [];
  if (row.robust) parts.push("<span class=\\"badge robust\\">robust</span>");
  if (row.coloc_strengthened) parts.push("<span class=\\"badge coloc\\">coloc</span>");
  return parts.join(" ");
}
function detailRow(row) {
  return `<div class="detail"><h4>${row.phenotype} <small>(${row.phecode}, ${row.panel})</small></h4>
    <dl>
      <dt>IVW</dt><dd>β=${row.b_ivw ?? "NA"} ± ${row.se_ivw ?? "NA"} (p=${pvalCell(row.p_ivw)})</dd>
      <dt>Weighted Median</dt><dd>β=${row.b_wm ?? "NA"} (p=${pvalCell(row.p_wm)})</dd>
      <dt>MR-Egger</dt><dd>β=${row.b_egger ?? "NA"} (p=${pvalCell(row.p_egger)}), intercept p=${pvalCell(row.egger_int_p)}</dd>
      <dt>Cochran Q</dt><dd>${row.q_stat ?? "NA"} (p=${pvalCell(row.q_p)})</dd>
      <dt>Within-panel FDR q</dt><dd>${row.fdr_q ?? "NA"}</dd>
      <dt>Cases / IVs</dt><dd>${row.n_cases ?? "NA"} cases; ${row.n_iv} instruments used</dd>
      ${row.coloc_pp_h4 !== null ? `<dt>Top coloc</dt><dd>PP.H4=${row.coloc_pp_h4} at <code>${row.coloc_gene ?? "?"}</code> (${row.coloc_tissue ?? "?"})</dd>` : ""}
    </dl>
  </div>`;
}

function populateSelect(id, items) {
  const sel = document.getElementById(id);
  for (const v of items) {
    const opt = document.createElement("option");
    opt.value = v; opt.textContent = v;
    sel.appendChild(opt);
  }
}
populateSelect("panel-filter", PANELS);
populateSelect("cat-filter", CATEGORIES);

$(document).ready(function() {
  const table = $("#atlas").DataTable({
    data: ATLAS,
    pageLength: 25,
    order: [[8, "asc"]],
    columns: [
      { className: "dt-control", orderable: false, data: null, defaultContent: "▸" },
      { data: "phecode" },
      { data: "phenotype" },
      { data: "category" },
      { data: "panel" },
      { data: "n_iv" },
      { data: "b_ivw", render: (x) => x === null ? "" : x.toFixed(3) },
      { data: "se_ivw", render: (x) => x === null ? "" : x.toFixed(3) },
      { data: "p_ivw", render: (x) => pvalCell(x) },
      { data: "fdr_q", render: (x) => pvalCell(x) },
      { data: "coloc_pp_h4", render: (x) => x === null ? "" : x.toFixed(2) },
      { data: null, render: (d, t, r) => flagCell(r) }
    ]
  });

  // Custom filters
  $.fn.dataTable.ext.search.push(function(settings, data, idx) {
    const row = table.row(idx).data();
    const pf = document.getElementById("panel-filter").value;
    const cf = document.getElementById("cat-filter").value;
    if (pf && row.panel !== pf) return false;
    if (cf && row.category !== cf) return false;
    if (document.getElementById("robust-only").checked && !row.robust) return false;
    if (document.getElementById("coloc-only").checked && !row.coloc_strengthened) return false;
    if (document.getElementById("sig-only").checked && !(row.p_ivw < 1e-4)) return false;
    return true;
  });
  $("#panel-filter, #cat-filter, #robust-only, #coloc-only, #sig-only").on("change", () => table.draw());

  // Row expansion
  $("#atlas tbody").on("click", "td.dt-control", function() {
    const tr = $(this).closest("tr");
    const row = table.row(tr);
    if (row.child.isShown()) { row.child.hide(); tr.removeClass("shown"); this.textContent = "▸"; }
    else { row.child(detailRow(row.data())).show(); tr.addClass("shown"); this.textContent = "▾"; }
  });
});
</script>

</body>
</html>')

writeLines(html, out_html)

# Get file size
size_mb <- file.size(out_html) / 1024 / 1024
cat(sprintf("[%s] Wrote %s (%.2f MB)\n", Sys.time(), out_html, size_mb))
cat(sprintf("Open in any browser to use: file://%s\n", normalizePath(out_html)))

# Also write a TSV companion for download
tsv_fp <- file.path(dirname(out_html), "amyloid_pet_causal_atlas.tsv")
fwrite(atlas, tsv_fp, sep = "\t")
cat(sprintf("[%s] Wrote TSV companion: %s\n", Sys.time(), tsv_fp))
