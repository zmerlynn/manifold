#!/usr/bin/env bash
#
# Runs each CrossSectionFuzz target for FUZZ_PER (default 3m) back-to-back,
# tolerating per-target failures (Clipper2 third-party UBSan noise, etc.)
# so one bad target doesn't abort the sweep. Logs per-target output and
# emits a master.log with "==> target start" / "<== target exit=N" markers
# that downstream tooling (the boolean2-robustness-iteration skill) parses
# to find new failures since its last iteration.
#
# Used by:
#   - the boolean2-robustness-iteration skill, as the bootstrap when no
#     fuzz process is running
#   - manual local invocation when you want background fuzz without the
#     CI overhead
#
# Differs from fuzz-cross-section.sh: that one runs ctest and aborts on
# the first failure (suitable for CI smoke). This one drives the binary
# directly, loops over all targets, and ignores per-target failures.

set -u

BINARY="${CROSS_SECTION_FUZZ_BINARY:-./build/codex-fuzz-ubsan-debug-local/test/cross_section_fuzz}"
CORPUS="${FUZZ_CORPUS:-./build/fuzz-corpus}"
STAMP="${FUZZ_STAMP:-$(date +%Y-%m-%d)}"
LOG_DIR="${FUZZ_LOG_DIR:-/tmp/fuzz-log-${STAMP}}"
# libFuzzer wants seconds, not "3m"-style. Default 180 = 3 min.
FUZZ_PER_SECONDS="${FUZZ_PER_SECONDS:-180}"

mkdir -p "$LOG_DIR" "$CORPUS"

export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_container_overflow=0:halt_on_error=0:abort_on_error=0}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=0:print_stacktrace=0}"

# Targets are hardcoded (matches .github/workflows/fuzz_continuous.yml
# matrix; keep in sync if a new FUZZ_TEST is added to
# test/cross_section_fuzz.cpp).
targets=(
  CrossSectionFuzz.BooleanRobustness
  CrossSectionFuzz.OffsetRobustness
  CrossSectionFuzz.ManifoldExtrudeRoundTrip
  CrossSectionFuzz.ManifoldSimpleExtrudeRoundTrip
  CrossSectionFuzz.BooleanExtrudeRoundTrip
  CrossSectionFuzz.PrismBooleanMatchesCrossSection
  CrossSectionFuzz.ManifoldTransformedExtrudeRoundTrip
  CrossSectionFuzz.DecomposeComposeAndHull
  CrossSectionFuzz.BatchBooleanSeparated
  CrossSectionFuzz.DecomposedExtrusionsRecompose
  CrossSectionFuzz.TranslatedExtrusionSliceMatchesCrossSection
  CrossSectionFuzz.RotatedExtrusionSliceMatchesCrossSection
  CrossSectionFuzz.HoledBooleanExtrudeRoundTrip
  CrossSectionFuzz.WarpAffineEquivalence
  CrossSectionFuzz.MirrorExtrudeRoundTrip
  CrossSectionFuzz.SimpleBooleanIdentities
  CrossSectionFuzz.ApexSkipNearLine
  CrossSectionFuzz.TranslationInvariance
  CrossSectionFuzz.IterateToFixedPointConverges
  CrossSectionFuzz.SubtractInvariants
  CrossSectionFuzz.BooleanCommutativity
  CrossSectionFuzz.BooleanAssociativity
  CrossSectionFuzz.BooleanDistributivity
  CrossSectionFuzz.ScaleInvariance
  CrossSectionFuzz.RotationInvariance
  CrossSectionFuzz.OffsetIdentityAtZero
  CrossSectionFuzz.EmptyIdentities
  CrossSectionFuzz.DoubleMirrorIdentity
  CrossSectionFuzz.DecomposeRecomposeWithHoles
  CrossSectionFuzz.OffsetInverseConvex
  CrossSectionFuzz.HullIdempotence
  CrossSectionFuzz.VertexMergeIdempotence
  CrossSectionFuzz.BVHPairEnumerationMatchesBruteForce
  CrossSectionFuzz.CanonicalSubEdgeIdempotence
  CrossSectionFuzz.SimplePositiveOffset
  CrossSectionFuzz.OffsetExtrudeRoundTrip
  CrossSectionFuzz.ManifoldHoledExtrudeRoundTrip
)

# Optional override: pass target names as args to restrict to a subset.
if [[ $# -gt 0 ]]; then
  targets=("$@")
fi

if [[ ! -x "$BINARY" ]]; then
  echo "ERROR: fuzz binary not found at $BINARY" >&2
  echo "Build it first, e.g.:" >&2
  echo "  cmake --build build/codex-fuzz-ubsan-debug-local --target cross_section_fuzz" >&2
  exit 1
fi

echo "Sweep config:"
echo "  binary:    $BINARY"
echo "  corpus:    $CORPUS (per-target subdirs)"
echo "  log dir:   $LOG_DIR"
echo "  per-target: ${FUZZ_PER_SECONDS}s"
echo "  targets:   ${#targets[@]}"
echo

# libFuzzer invocation: positional corpus dir per target (accumulates
# across runs), -max_total_time in seconds, -artifact_prefix per target.
# Matches the CI pipeline's per-target corpus layout so we can share
# corpus across local and CI (in principle - they don't actually swap
# corpora today, but the format is the same).
for t in "${targets[@]}"; do
  echo "==> $t ($(date +%H:%M:%S))"
  mkdir -p "$CORPUS/$t" "$LOG_DIR/$t-artifacts"
  "$BINARY" --fuzz="$t" \
    "$CORPUS/$t" \
    -max_total_time="$FUZZ_PER_SECONDS" \
    -timeout=60 \
    -artifact_prefix="$LOG_DIR/$t-artifacts/" \
    > "$LOG_DIR/$t.log" 2>&1
  echo "<== $t exit=$? ($(date +%H:%M:%S))"
done
echo "DONE $(date +%H:%M:%S)"
