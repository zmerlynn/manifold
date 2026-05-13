#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export CROSS_SECTION_FUZZ_BINARY="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz-ubsan-debug/test/cross_section_fuzz}"
export FUZZ_CORPUS="${FUZZ_CORPUS:-build/fuzz-corpus-ubsan}"
export FUZZ_FOR="${FUZZ_FOR:-30s}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-print_stacktrace=1}"

exec "${script_dir}/fuzz-cross-section.sh" "$@"
