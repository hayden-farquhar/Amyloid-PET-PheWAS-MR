#!/usr/bin/env python3
"""Build manuscript/references.bib from CrossRef + PubMed lookups.

Used in lieu of a reference manager — Hayden does not use Zotero / EndNote /
Mendeley. Per portfolio memory feedback_no_reference_manager.md, this is
Claude-side in-context work during submission prep.

Strategy:
  1. Curated list of (citation_key, DOI-or-search-query, citation_kind)
  2. For each: if DOI known, fetch CrossRef directly; if not, search CrossRef
     for the title+author+year and take the top hit.
  3. Format response as BibTeX entry.
  4. Write to manuscript/references.bib.

The DOI guesses below come from well-known reference databases (Ge 2019 PRS-CS
on Nat Commun 10:1776, etc.); they are verified against the CrossRef return
metadata in the run output so any wrong-DOI hit shows up immediately.
"""

from __future__ import annotations
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Reference catalogue
# Format: citation_key -> dict with EITHER 'doi' (preferred) OR 'query'
# (fallback for cases where the DOI isn't easily known).
# ---------------------------------------------------------------------------
REFS: dict[str, dict] = {
    # --- Drug trials ---
    "vandyck-2023-clarity": {"doi": "10.1056/NEJMoa2212948"},
    "sims-2023-trailblazer": {"doi": "10.1001/jama.2023.13239"},
    # --- AD GWAS ---
    "bellenguez-2022-adgwas": {"doi": "10.1038/s41588-022-01024-z"},
    "lambert-2013-adgwas": {"doi": "10.1038/ng.2802"},
    "ali-2023-amyloidpet": {"doi": "10.1186/s40478-023-01563-4"},
    "kim-2025-koreanamyloid": {"doi": "10.1038/s41467-025-57751-4"},
    "farrer-1997-apoe": {"doi": "10.1001/jama.1997.03550160069041"},
    "liu-2013-apoe": {"doi": "10.1038/nrneurol.2012.263"},
    # --- CSF biomarkers ---
    "timsina-2026-csf": {"doi": "10.1038/s41467-026-71682-8"},
    "jansen-2022-csf": {"doi": "10.1007/s00401-022-02454-z"},
    # --- MR methods ---
    "burgess-2013-iv": {"doi": "10.1002/gepi.21758"},  # Mendelian randomization with multiple genetic variants
    "burgess-2017-mregger": {"doi": "10.1007/s10654-017-0255-x"},  # MR-Egger interpretation
    "bowden-2015-egger": {"doi": "10.1093/ije/dyv080"},
    "bowden-2016-weightedmedian": {"doi": "10.1002/gepi.21965"},
    "hartwig-2017-mbe": {"doi": "10.1093/ije/dyx102"},
    "verbanck-2018-mrpresso": {"doi": "10.1038/s41588-018-0099-7"},
    "zhao-2020-mrraps": {"doi": "10.1214/19-aos1866"},
    "daveysmith-hemani-2014-mr": {"doi": "10.1093/hmg/ddu328"},
    "skrivankova-2021-strobemr": {"doi": "10.1001/jama.2021.18236"},
    "lee-2012-liability": {"doi": "10.1002/gepi.21614"},
    # --- Colocalisation ---
    "wakefield-2009-bf": {"doi": "10.1002/gepi.20359"},
    "giambartolomei-2014-coloc": {"doi": "10.1371/journal.pgen.1004383"},
    "hukku-2021-colocldmismatch": {"doi": "10.1016/j.ajhg.2020.11.012"},
    # --- CR1 / APOE pleiotropy biology ---
    "crehan-2012-cr1": {"doi": "10.1016/j.imbio.2011.07.017"},
    "itzhaki-2018-hsv": {"doi": "10.3389/fnagi.2018.00324"},
    "mahley-rall-2000-apoe": {"doi": "10.1146/annurev.genom.1.1.507"},
    # --- PRS-CS ---
    "ge-2019-prscs": {"doi": "10.1038/s41467-019-09718-5"},
    # --- Added 2026-05-22 after targeted CrossRef DOI verification ---
    "holmes-2017-cardiometabolic-mr": {"doi": "10.1038/nrcardio.2017.78"},     # lipid-MR analog ref for Discussion §4
    "burgess-thompson-2021-mrbook":   {"doi": "10.1201/9780429324352"},        # MR textbook 2nd ed; covers R² formula used in Methods Steiger §86
    "schilling-2005-apoe-bone":       {"doi": "10.1359/JBMR.041101"},          # APOE-bone pleiotropy for Discussion §37 limitation 1 (ii)
    "klunk-2004-pib":                 {"doi": "10.1002/ana.20009"},            # PiB-PET tracer for Intro footnote ¹
    "wong-2010-florbetapir":          {"doi": "10.2967/jnumed.109.069088"},    # florbetapir tracer for Intro footnote ¹
    # NOTE: Cauley 2007 (APOE-bone) and Lim 2008 (APOE-cognitive ageing) remain
    # deferred. Both are tangential supporting citations in the APOE-pleiotropy
    # paragraph (Intro §7 + Discussion §37 hypothesis ii). CrossRef searches
    # return unrelated papers under those author surnames at those years; the
    # specific papers Hayden intended need user-side identification at
    # submission packaging. Mahley & Rall 2000 + Itzhaki 2018 + Schilling 2005
    # already cover APOE-pleiotropy in the same paragraphs, so removing the
    # Cauley/Lim inline citations is a defensible alternative.
}

USER_AGENT = "P81-AmyloidPETPheWASMR/1.0 (mailto:hayden.farquhar@icloud.com)"
CROSSREF_DOI = "https://api.crossref.org/works/{doi}"
CROSSREF_SEARCH = "https://api.crossref.org/works?query={q}&rows=3"

def fetch(url: str) -> dict | None:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"  ERROR: {e}", file=sys.stderr)
        return None

def lookup_doi(doi: str) -> dict | None:
    data = fetch(CROSSREF_DOI.format(doi=doi))
    return data["message"] if data else None

def lookup_query(q: str) -> dict | None:
    data = fetch(CROSSREF_SEARCH.format(q=urllib.parse.quote(q)))
    if not data or not data.get("message", {}).get("items"):
        return None
    return data["message"]["items"][0]  # top hit

def fmt_authors(authors: list[dict] | None) -> str:
    if not authors:
        return ""
    parts = []
    for a in authors:
        family = a.get("family", "")
        given = a.get("given", "")
        if family:
            parts.append(f"{family}, {given}" if given else family)
    return " and ".join(parts)

def to_bibtex(key: str, m: dict) -> str:
    # entry type
    t = m.get("type", "journal-article")
    entry_type = {
        "journal-article": "article",
        "book": "book",
        "reference-book": "book",
        "edited-book": "book",
        "monograph": "book",
        "book-chapter": "incollection",
        "proceedings-article": "inproceedings",
        "posted-content": "misc",
    }.get(t, "article")

    fields = {}
    fields["author"] = fmt_authors(m.get("author"))
    title = m.get("title", [""])
    fields["title"] = title[0] if title else ""
    container = m.get("container-title", [""])
    fields["journal"] = container[0] if container else ""
    issued = m.get("issued", {}).get("date-parts", [[None]])[0]
    if issued and issued[0]:
        fields["year"] = str(issued[0])
    if m.get("volume"):
        fields["volume"] = str(m["volume"])
    if m.get("issue"):
        fields["number"] = str(m["issue"])
    if m.get("page"):
        fields["pages"] = m["page"].replace("-", "--")
    if m.get("DOI"):
        fields["doi"] = m["DOI"]
    if m.get("URL"):
        fields["url"] = m["URL"]

    # Drop empty fields
    fields = {k: v for k, v in fields.items() if v}

    lines = [f"@{entry_type}{{{key},"]
    max_k = max(len(k) for k in fields) if fields else 8
    for k, v in fields.items():
        # Escape braces and protect title casing
        v_str = str(v).replace("{", "\\{").replace("}", "\\}")
        if k == "title":
            v_str = "{" + v_str + "}"  # protect against BibTeX lowercasing
        lines.append(f"  {k:<{max_k}} = {{{v_str}}},")
    # Strip trailing comma on last field
    lines[-1] = lines[-1].rstrip(",")
    lines.append("}")
    return "\n".join(lines)

def main() -> int:
    proj_root = Path("/Users/haydenfarquhar/Documents/Research Projects/81 Amyloid PET PheWAS MR")
    out_fp = proj_root / "manuscript" / "references.bib"

    entries = []
    failures = []
    for i, (key, spec) in enumerate(REFS.items(), 1):
        method = "DOI" if "doi" in spec else "search"
        target = spec.get("doi") or spec.get("query")
        print(f"[{i}/{len(REFS)}] {key} via {method}: {target}")
        m = lookup_doi(spec["doi"]) if "doi" in spec else lookup_query(spec["query"])
        if not m:
            print(f"  FAIL")
            failures.append((key, target))
            continue
        # Diagnostic: show what was fetched
        title = (m.get("title") or [""])[0][:80]
        year = m.get("issued", {}).get("date-parts", [[None]])[0][0]
        first_author = (m.get("author") or [{}])[0].get("family", "?")
        print(f"  -> {first_author} {year}: {title}")
        entries.append(to_bibtex(key, m))
        time.sleep(0.3)  # be polite to CrossRef

    header = (
        "% references.bib — P81 Amyloid PET PheWAS MR manuscript\n"
        "% Built via CrossRef API lookups (Claude session 2026-05-22).\n"
        "% Generator: scripts/build_references_bib.py\n"
        f"% Total entries: {len(entries)} (of {len(REFS)} requested)\n"
        "% Re-run the script after adding/correcting citation keys in the script.\n"
        "\n"
    )

    out_fp.write_text(header + "\n\n".join(entries) + "\n")
    print(f"\nWrote {out_fp} with {len(entries)} entries.")
    if failures:
        print("\nFailed lookups (need manual entry):")
        for k, t in failures:
            print(f"  {k}: {t}")
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
