#!/usr/bin/env bash

set -euo pipefail

binary="${CROSS_SECTION_FUZZ_BINARY:-./build/fuzz/test/cross_section_fuzz}"
filter="${GTEST_FILTER:-CrossSectionFuzz.*}"

ASAN_OPTIONS="${ASAN_OPTIONS:-detect_container_overflow=0}" \
  "${binary}" \
  --gtest_filter="${filter}" \
  "$@"
