#!/usr/bin/env bash
# Sweep FinnGen R12 sumstats for 1,574 phecodes via TABIX byte-range queries.
#
# Strategy: each FinnGen R12 .gz is bgzip-compressed with a .gz.tbi index
# served alongside it on GCS. We hand tabix the remote URL + a regions file
# (one BED entry per instrument SNP) and get back only the rows we want
# (~17 rows / ~1 KB) per phecode in ~24 seconds — vs ~270 seconds streaming
# the full 818 MB file.
#
# Run as:
#   bash scripts/sweep_finngen_r12.sh OUT_DIR INSTRUMENT_TSV [PARALLEL_WORKERS]
# Where:
#   OUT_DIR        = destination, e.g. "$HOME/.../_shared_mirror/finngen_R12/"
#   INSTRUMENT_TSV = TSV with rsid, chr, hg38_pos columns (cols 1, 2, 3)
#   PARALLEL_WORKERS = optional, default 4
#
# Resumable: existing valid slices are skipped on rerun.
#
# Wall-clock estimate: ~2.5 hours with 4 parallel workers.

set -uo pipefail

OUT_DIR="${1:?Usage: $0 OUT_DIR INSTRUMENT_TSV [PARALLEL_WORKERS]}"
INSTR_TSV="${2:?Usage: $0 OUT_DIR INSTRUMENT_TSV [PARALLEL_WORKERS]}"
WORKERS="${3:-4}"
PHECODE_CSV="osf/phecode_filter_prespec.csv"

command -v tabix >/dev/null || { echo "ERROR: tabix not found; install via 'brew install htslib'"; exit 1; }

mkdir -p "$OUT_DIR/sliced" "$OUT_DIR/logs" "$OUT_DIR/.tabix_cache"
# Tabix downloads remote .tbi indexes to the CURRENT WORKING DIRECTORY by default,
# which spammed 1,574 .tbi files (~2.3 GB) into the project root on first run.
# Isolate them in a dedicated cache dir that gets purged at end of run.
TABIX_CACHE="$OUT_DIR/.tabix_cache"
trap 'rm -f "$TABIX_CACHE"/*.tbi 2>/dev/null' EXIT

# Build regions BED from instrument TSV (cols: rsid chr hg38_pos)
# Output BED is 0-based half-open per tabix convention.
# tabix requires the regions file sorted by chr then start position; sort here
# to ensure correctness regardless of the input TSV's ordering.
REGIONS="$OUT_DIR/.iv_regions.bed"
awk 'NR>1 {pos=$3; print $2"\t"(pos-1)"\t"pos}' "$INSTR_TSV" \
  | sort -k1,1n -k2,2n > "$REGIONS"
N_REGIONS=$(wc -l < "$REGIONS" | tr -d ' ')

# Build the (phecode, url) work list for xargs
WORK_LIST="$OUT_DIR/.work_list.txt"
tail -n +2 "$PHECODE_CSV" | tr -d '\r' | while IFS=',' read -r phecode rest; do
  url=$(echo "$rest" | awk -F',' '{print $NF}' | tr -d '"')
  echo -e "$phecode\t$url"
done > "$WORK_LIST"
N_PHECODES=$(wc -l < "$WORK_LIST" | tr -d ' ')

LOG="$OUT_DIR/logs/sweep.log"
RUN_TAG="run-$(date -u +%Y%m%d-%H%M%S)"
echo "[$(date -u +%FT%TZ)] ====== $RUN_TAG: $N_PHECODES phecodes × $N_REGIONS instruments × $WORKERS workers ======" >> "$LOG"
echo "Tabix sweep: $N_PHECODES phecodes × $N_REGIONS instruments × $WORKERS workers"
echo "Run tag: $RUN_TAG (use 'tail -f $LOG' to follow)"

# Per-phecode worker: tabix-fetch + validate + gzip-write + log
export OUT_DIR REGIONS LOG TABIX_CACHE
fetch_one() {
  local phecode="$1"
  local url="$2"
  local out="$OUT_DIR/sliced/${phecode}.tsv.gz"

  # Resume: skip if existing slice is valid
  if [[ -s "$out" ]] && gunzip -c "$out" 2>/dev/null | head -1 | grep -q '^#chrom'; then
    echo "[$(date -u +%FT%TZ)] SKIP $phecode" >> "$LOG"
    return 0
  fi

  # Tabix query: header + regions. Capture both stdout and exit code.
  # Subshell + cd isolates the downloaded .tbi to TABIX_CACHE rather than CWD.
  local tmp="$OUT_DIR/sliced/.${phecode}.partial"
  if ( cd "$TABIX_CACHE" && tabix -h "$url" -R "$REGIONS" ) 2>/dev/null > "$tmp"; then
    # Validate: at minimum the FinnGen header must be present
    if head -1 "$tmp" 2>/dev/null | grep -q '^#chrom'; then
      local nrows
      nrows=$(wc -l < "$tmp" | tr -d ' ')
      gzip -c < "$tmp" > "$out" && rm -f "$tmp"
      echo "[$(date -u +%FT%TZ)] OK $phecode ($nrows rows)" >> "$LOG"
    else
      echo "[$(date -u +%FT%TZ)] FAIL_HEADER $phecode" >> "$LOG"
      rm -f "$tmp" "$out"
    fi
  else
    echo "[$(date -u +%FT%TZ)] FAIL_TABIX $phecode" >> "$LOG"
    rm -f "$tmp" "$out"
  fi
}
export -f fetch_one

# Run in parallel via xargs (-P workers, -L 1 row per invocation)
cat "$WORK_LIST" \
  | xargs -P "$WORKERS" -L 1 bash -c 'fetch_one "$0" "$1"'

# Summary — count only this run's entries (since the run-start marker)
ok=$(awk "/$RUN_TAG/{seen=1} seen && / OK /{n++} END{print n+0}" "$LOG")
skip=$(awk "/$RUN_TAG/{seen=1} seen && / SKIP /{n++} END{print n+0}" "$LOG")
fail=$(awk "/$RUN_TAG/{seen=1} seen && / FAIL_(TABIX|HEADER) /{n++} END{print n+0}" "$LOG")

echo "[$(date -u +%FT%TZ)] Sweep complete: ok=$ok skip=$skip fail=$fail (of $N_PHECODES)" >> "$LOG"
echo
echo "Sweep complete."
echo "  Total phecodes: $N_PHECODES"
echo "  OK:             $ok"
echo "  Skipped (prior): $skip"
echo "  Failed:         $fail"
echo "  Output dir:     $OUT_DIR/sliced/"
echo
echo "Re-run this script to retry FAILs; OK slices are skipped."
