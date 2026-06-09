#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export CROSS_SECTION_FUZZ_BINARY="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz/test/cross_section_fuzz}"
export FUZZ_CORPUS="${FUZZ_CORPUS:-build/fuzz-corpus}"
export FUZZ_FOR="${FUZZ_FOR:-1m}"
export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_container_overflow=0}"

exec "${script_dir}/fuzz-cross-section.sh" "$@"
