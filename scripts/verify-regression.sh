#!/usr/bin/env bash
#
# scripts/verify-regression.sh - clean-rebuild + verify a single regression
# test, with paranoia about stale state. Run this BEFORE claiming a fix is
# landed or that a test passes/fails.
#
# Why: incremental builds with broken dep-tracking can silently leave object
# files stale - the test binary ends up compiled against older source than
# what's on disk, and gtest filters can hit ghost symbols (e.g. a still-
# DISABLED_-prefixed name from before a rename). That happened on 2026-05-18
# in the DecomposeRecomposeOuterStarWithSmallHole verification: I told the
# user and the consumer the fix was broken when actually my binary was just
# stale, wasted a chunk of consumer time.
#
# Guardrails this script enforces:
#   1. Refuse to run on dirty working tree (commit/stash first, or --force).
#   2. Wipe a dedicated build/verify-regression dir; cmake reconfigures
#      from scratch. No inherited dep-tracking corruption.
#   3. Verify the named test actually exists in the rebuilt binary before
#      running it. Catches the "I think I renamed it but didn't" class.
#   4. Print HEAD sha + binary mtimes + result, so the paper trail is clear.
#
# Usage:
#   scripts/verify-regression.sh --filter <gtest_filter> [--branch <name>] [--force]
#
# Examples:
#   scripts/verify-regression.sh --filter "CrossSection.DecomposeRecomposeOuterStarWithSmallHole"
#   scripts/verify-regression.sh --filter "CrossSection.DISABLED_Foo*" --branch pr/boolean2-tests
#
# Flags:
#   --filter <pattern>  gtest filter (required). --gtest_also_run_disabled_tests
#                       is always on, so DISABLED_-prefixed tests run too.
#   --branch <name>     documentation only - records which branch was tested.
#                       Script does NOT switch branches; do that yourself.
#   --force             skip the dirty-tree check.

set -euo pipefail

FILTER=""
BRANCH_HINT=""
FORCE=""
BUILD_DIR="build/verify-regression"

usage() {
  sed -n '/^# Usage:/,/^# Flags:/{/^# Flags:/!p}' "$0" | sed 's/^# *//' >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --filter) FILTER="${2:-}"; shift 2;;
    --branch) BRANCH_HINT="${2:-}"; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) usage;;
    *) echo "unknown arg: $1" >&2; usage;;
  esac
done

[[ -z "$FILTER" ]] && usage

# Guardrail 1: dirty TRACKED source files only (untracked files in
# test/polygons/ etc. don't affect the build, ignore them).
if [[ -z "$FORCE" ]] && ! git diff-index --quiet HEAD -- 'src/' 'test/' 'include/' 'CMakeLists.txt' 2>/dev/null; then
  echo "ERROR: working tree has uncommitted source changes. Commit/stash or pass --force." >&2
  git diff --stat HEAD -- 'src/' 'test/' 'include/' 'CMakeLists.txt' >&2
  exit 1
fi

# Record state
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
CURRENT_SHA=$(git rev-parse HEAD)
echo "verify-regression: filter='$FILTER'"
echo "  branch:        $CURRENT_BRANCH (hint: ${BRANCH_HINT:-none})"
echo "  HEAD sha:      $CURRENT_SHA"

# Warn if local is behind/ahead of origin
git fetch origin --quiet "$CURRENT_BRANCH" 2>/dev/null || true
if git rev-parse --verify --quiet "origin/$CURRENT_BRANCH" >/dev/null; then
  REMOTE_SHA=$(git rev-parse "origin/$CURRENT_BRANCH")
  if [[ "$REMOTE_SHA" != "$CURRENT_SHA" ]]; then
    echo "  NOTE: local HEAD differs from origin/$CURRENT_BRANCH:"
    echo "    local:  $CURRENT_SHA"
    echo "    origin: $REMOTE_SHA"
  fi
fi

# Guardrail 2: wipe build dir for fresh configure
echo
echo "wiping $BUILD_DIR for fresh configure..."
rm -rf "$BUILD_DIR"

cmake -B "$BUILD_DIR" \
  -DMANIFOLD_CROSS_SECTION=ON \
  -DMANIFOLD_CROSS_SECTION_BACKEND=boolean2 \
  -DMANIFOLD_TEST=ON \
  -DMANIFOLD_DEBUG=ON \
  . 2>&1 | tail -3

echo
echo "building manifold_test..."
cmake --build "$BUILD_DIR" --target manifold_test -j"$(nproc)" 2>&1 \
  | tail -5

BINARY="$BUILD_DIR/test/manifold_test"
if [[ ! -x "$BINARY" ]]; then
  echo "ERROR: manifold_test not built at $BINARY" >&2
  exit 1
fi

LIB_MTIME=$(stat -c '%y' "$BUILD_DIR/src/libmanifold.so" 2>/dev/null || echo "unknown")
BIN_MTIME=$(stat -c '%y' "$BINARY")
echo "  libmanifold.so built: $LIB_MTIME"
echo "  manifold_test built:  $BIN_MTIME"

# Guardrail 3: verify the test name actually exists in the binary
echo
echo "checking '$FILTER' is registered in the binary..."
LISTED_COUNT=$("$BINARY" \
                 --gtest_filter="$FILTER" \
                 --gtest_list_tests \
                 --gtest_also_run_disabled_tests 2>/dev/null \
               | grep -c '^[[:space:]]' || true)
if [[ "$LISTED_COUNT" == "0" ]]; then
  echo "ERROR: no test matching '$FILTER' found in $BINARY" >&2
  echo "  (check spelling, the DISABLED_ prefix, and that the test was actually" >&2
  echo "   compiled - if you just renamed it, this script's clean rebuild already" >&2
  echo "   reflects the rename, so a mismatch here is a real source-vs-filter bug.)" >&2
  exit 1
fi
echo "  $LISTED_COUNT test(s) match filter"

# Run it
echo
echo "running test..."
LOG=$(mktemp /tmp/verify-regression-XXXXXX.log)
trap "rm -f $LOG" EXIT
if "$BINARY" \
    --gtest_filter="$FILTER" \
    --gtest_also_run_disabled_tests \
    > "$LOG" 2>&1; then
  RESULT="PASS"
else
  RESULT="FAIL"
fi
# Always show test output (gtest's own pass/fail summary is informative)
tail -20 "$LOG"

# Report
echo
echo "==================================================="
echo "verify-regression RESULT: $RESULT"
echo "  filter:        $FILTER"
echo "  matched count: $LISTED_COUNT"
echo "  branch:        $CURRENT_BRANCH"
echo "  HEAD sha:      $CURRENT_SHA"
echo "  binary mtime:  $BIN_MTIME"
echo "==================================================="

[[ "$RESULT" == "PASS" ]] || exit 1
