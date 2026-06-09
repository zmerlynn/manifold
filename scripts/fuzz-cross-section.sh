#!/usr/bin/env bash

set -euo pipefail

binary="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz/test/cross_section_fuzz}"
corpus="${FUZZ_CORPUS:-build/fuzz-corpus}"
duration="${FUZZ_FOR:-30s}"
binary_id="${FUZZ_BINARY_ID:-$(basename "${binary}")}"
coverage_root="${FUZZ_COVERAGE_ROOT:-${corpus}/${binary_id}}"
log_dir="${FUZZ_LOG_DIR:-}"

# Targets: positional args override; otherwise discover dynamically from the
# binary via --list_fuzz_tests. Keeps this entry point in sync with whatever
# FUZZ_TEST registrations are actually compiled (the previous hand-curated
# list drifted to 19 of 43 registered targets before being caught in review).
if [[ $# -gt 0 ]]; then
  targets=("$@")
else
  if [[ ! -x "${binary}" ]]; then
    echo "ERROR: fuzz binary not found at ${binary}; cannot enumerate targets" >&2
    exit 1
  fi
  mapfile -t targets < <(
    "${binary}" --list_fuzz_tests 2>/dev/null \
      | sed -n 's/^\[\*\] Fuzz test: //p'
  )
  if [[ ${#targets[@]} -eq 0 ]]; then
    echo "ERROR: --list_fuzz_tests returned no targets from ${binary}" >&2
    exit 1
  fi
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
