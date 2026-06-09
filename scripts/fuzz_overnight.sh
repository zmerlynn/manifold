#!/bin/bash
# Overnight fuzz wrapper for overlap3d_proto --advfuzz.
# Loops a range of seeds, captures per-seed summaries, saves
# adversarial inputs (= unchanged / worsened / fallback cases) for
# later analysis.
#
# Usage:
#   scripts/fuzz_overnight.sh START END OUTDIR
# Defaults: START=1, END=5000, OUTDIR=/tmp/fuzz_overnight

set -u
START="${1:-1}"
END="${2:-5000}"
OUTDIR="${3:-/tmp/fuzz_overnight}"

PROTO=/home/zml/src/manifold-overlap/build/extras/overlap3d_proto
SUMMARY="$OUTDIR/seed_summary.csv"
PROGRESS="$OUTDIR/progress.log"

mkdir -p "$OUTDIR/cases"
echo "seed,total,valid,pierces,fixed,reduced,regularized,unchanged,worsened,fallback" > "$SUMMARY"

START_TIME=$(date +%s)
for seed in $(seq "$START" "$END"); do
  out=$(OVERLAP3D_FUZZ_SEED="$seed" OVERLAP3D_FUZZ_SAVE=1 \
        OVERLAP3D_FUZZ_OUTDIR="$OUTDIR/cases" \
        taskset -c 2,3 "$PROTO" --advfuzz 2>&1)

  # Parse the TOTAL summary block from the spike output.
  total=$(echo "$out" | grep -oE "TOTAL self-pierce cases: [0-9]+/[0-9]+" \
          | head -1 | tr -d '\n')
  pierces=$(echo "$out" | grep -oE "TOTAL self-pierce cases: [0-9]+" \
            | head -1 | awk '{print $NF}')
  valid=$(echo "$out" | grep -oE "TOTAL self-pierce cases: [0-9]+/[0-9]+" \
          | head -1 | awk -F/ '{print $NF}')
  fixed=$(echo "$out" | grep -oE "fixed \(post=0\):\s+[0-9]+" \
          | head -1 | awk '{print $NF}')
  reduced=$(echo "$out" | grep -oE "reduced \(post<pre\):[0-9]+" \
            | head -1 | tr -dc 0-9)
  regularized=$(echo "$out" | grep -oE "regularized-to-empty.*: [0-9]+" \
                | head -1 | awk '{print $NF}')
  unchanged=$(echo "$out" | grep -oE "unchanged:\s+[0-9]+" \
              | head -1 | awk '{print $NF}')
  worsened=$(echo "$out" | grep -oE "worsened:\s+[0-9]+" \
             | head -1 | awk '{print $NF}')
  fallback=$(echo "$out" | grep -oE "fallback \(.*\): [0-9]+" \
             | head -1 | awk '{print $NF}')

  echo "$seed,2495,${valid:-0},${pierces:-0},${fixed:-0},${reduced:-0},${regularized:-0},${unchanged:-0},${worsened:-0},${fallback:-0}" >> "$SUMMARY"

  # Progress every 10 seeds.
  if (( seed % 10 == 0 )); then
    elapsed=$(( $(date +%s) - START_TIME ))
    saved=$(ls "$OUTDIR/cases" 2>/dev/null | wc -l)
    echo "[$(date +%H:%M:%S)] seed=$seed elapsed=${elapsed}s saved-cases=$saved" \
         | tee -a "$PROGRESS"
  fi
done

echo "[$(date +%H:%M:%S)] DONE. Summary: $SUMMARY" >> "$PROGRESS"
