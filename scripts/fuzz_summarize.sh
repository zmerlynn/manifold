#!/bin/bash
# Summarize results from scripts/fuzz_overnight.sh.
# Usage: scripts/fuzz_summarize.sh [OUTDIR]
# Default OUTDIR=/tmp/fuzz_overnight

OUTDIR="${1:-/tmp/fuzz_overnight}"
SUMMARY="$OUTDIR/seed_summary.csv"

if [[ ! -f "$SUMMARY" ]]; then
  echo "ERROR: $SUMMARY not found" >&2
  exit 1
fi

echo "=== Overall totals ==="
awk -F, 'NR>1 {
  seeds++; pierces+=$4; fixed+=$5; reduced+=$6; reg+=$7;
  unch+=$8; wors+=$9; fb+=$10
} END {
  printf "seeds=%d  pierces=%d  fixed=%d (%.1f%%)  reduced=%d  regularized=%d\n",
    seeds, pierces, fixed, pierces ? 100.0*fixed/pierces : 0, reduced, reg
  printf "  unchanged=%d  worsened=%d  fallback=%d  fix-rate=%.2f%%\n",
    unch, wors, fb, pierces ? 100.0*(fixed+reduced+reg)/pierces : 0
}' "$SUMMARY"

echo ""
echo "=== Worsened cases (= pipeline INCREASED pierce count) ==="
awk -F, 'NR>1 && $9>0 {print $0}' "$SUMMARY" | head -20

echo ""
echo "=== Top 20 seeds by failure count (unchanged + worsened + fallback) ==="
awk -F, 'NR>1 {fails=$8+$9+$10; if (fails>0) print fails","$0}' "$SUMMARY" \
  | sort -t, -k1 -nr | head -20

echo ""
echo "=== Saved adversarial cases by class ==="
ls "$OUTDIR/cases" 2>/dev/null \
  | sed -E 's/.*_(class[0-9]+)_seed[0-9]+_(pre[0-9]+)_(post[0-9]+|fallback).obj/\1 \3/' \
  | sort | uniq -c | sort -rn | head -25

echo ""
echo "=== Worst (highest post-pierce) saved cases ==="
ls "$OUTDIR/cases" 2>/dev/null \
  | grep "_post" \
  | sed -E 's/.*_pre([0-9]+)_post([0-9]+).obj.*/\2 \1 &/' \
  | sort -k1 -nr | head -10 \
  | awk '{print "  post=" $1 " pre=" $2 " " $3}'

echo ""
echo "=== Total saved files ==="
ls "$OUTDIR/cases" 2>/dev/null | wc -l
