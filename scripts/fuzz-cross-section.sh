#!/usr/bin/env bash

set -euo pipefail

binary="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz/test/cross_section_fuzz}"
corpus="${FUZZ_CORPUS:-build/fuzz-corpus}"
duration="${FUZZ_FOR:-30s}"

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
  CrossSectionFuzz.SimplePositiveOffset
  CrossSectionFuzz.OffsetExtrudeRoundTrip
  CrossSectionFuzz.ManifoldHoledExtrudeRoundTrip
)

if [[ $# -gt 0 ]]; then
  targets=("$@")
fi

mkdir -p "${corpus}"

for target in "${targets[@]}"; do
  echo
  echo "==> ${target} (${duration})"
  ASAN_OPTIONS="${ASAN_OPTIONS:-detect_container_overflow=0}" \
    "${binary}" \
    --corpus_database="${corpus}" \
    --fuzz="${target}" \
    --fuzz_for="${duration}"
done
