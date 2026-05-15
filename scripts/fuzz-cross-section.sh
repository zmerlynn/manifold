#!/usr/bin/env bash

set -euo pipefail

binary="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz/test/cross_section_fuzz}"
corpus="${FUZZ_CORPUS:-build/fuzz-corpus}"
duration="${FUZZ_FOR:-30s}"
binary_id="${FUZZ_BINARY_ID:-$(basename "${binary}")}"
coverage_root="${FUZZ_COVERAGE_ROOT:-${corpus}/${binary_id}}"
log_dir="${FUZZ_LOG_DIR:-}"

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

mkdir -p "${corpus}" "${coverage_root}"
if [[ -n "${log_dir}" ]]; then
  mkdir -p "${log_dir}"
fi

quote_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

runner=()
if command -v stdbuf >/dev/null 2>&1; then
  runner=(stdbuf -oL -eL)
fi

echo "FuzzTest corpus database: ${corpus}"
echo "FuzzTest coverage output root: ${coverage_root}"
echo "  --corpus_database reads replay inputs from the database."
echo "  FUZZTEST_TESTSUITE_OUT_DIR writes new coverage-increasing inputs."
if [[ -n "${log_dir}" ]]; then
  echo "FuzzTest log directory: ${log_dir}"
fi

for target in "${targets[@]}"; do
  target_coverage="${coverage_root}/${target}/coverage"
  mkdir -p "${target_coverage}"
  cmd=("${binary}" "--corpus_database=${corpus}" "--fuzz=${target}"
       "--fuzz_for=${duration}")
  env_cmd=("ASAN_OPTIONS=${ASAN_OPTIONS:-detect_container_overflow=0}"
           "FUZZTEST_TESTSUITE_OUT_DIR=${target_coverage}")
  if [[ -n "${UBSAN_OPTIONS:-}" ]]; then
    env_cmd+=("UBSAN_OPTIONS=${UBSAN_OPTIONS}")
  fi

  echo
  echo "==> ${target} (${duration})"
  echo "coverage out: ${target_coverage}"
  echo "replay with: ${binary} --corpus_database=${corpus} --replay_corpus=${target}"
  echo -n "command: "
  quote_cmd env "${env_cmd[@]}" "${runner[@]}" "${cmd[@]}"

  if [[ -n "${log_dir}" ]]; then
    log_file="${log_dir}/${target}.log"
    echo "log: ${log_file}"
    env "${env_cmd[@]}" "${runner[@]}" "${cmd[@]}" 2>&1 | tee "${log_file}"
    status=${PIPESTATUS[0]}
    if [[ ${status} -ne 0 ]]; then
      exit "${status}"
    fi
  else
    env "${env_cmd[@]}" "${runner[@]}" "${cmd[@]}"
  fi
done
