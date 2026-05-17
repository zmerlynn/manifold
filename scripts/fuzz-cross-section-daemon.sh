#!/usr/bin/env bash
#
# Long-running local fuzz daemon for CrossSectionFuzz targets.
#
# Lives in its own clone of the repo (default: ~/src/manifold-fuzz-daemon),
# tracks origin/pr/cross-section-fuzz, rebuilds when source changes, and
# loops through every FUZZ_TEST in the binary FUZZ_PER_SECONDS each.
#
# Coordination is git-only: push to origin/pr/cross-section-fuzz from your
# producer clone, the daemon picks it up on the next loop iteration via
# `git fetch && git checkout origin/pr/cross-section-fuzz` followed by an
# incremental `cmake --build`.
#
# ============================================================
# One-time setup (run from your producer clone, not the daemon clone):
# ============================================================
#
#   cd ~/src
#   git clone https://github.com/zmerlynn/manifold.git manifold-fuzz-daemon
#   cd manifold-fuzz-daemon
#   git checkout pr/cross-section-fuzz
#   # libfuzzer-compat + ASan + UBSan + boolean2 backend:
#   cmake -B build/daemon \
#     -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
#     -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#     -DCMAKE_C_FLAGS="-fsanitize=undefined" \
#     -DCMAKE_CXX_FLAGS="-fsanitize=undefined" \
#     -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined" \
#     -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=undefined" \
#     -DMANIFOLD_FUZZ=ON \
#     -DFUZZTEST_COMPATIBILITY_MODE=libfuzzer \
#     -DMANIFOLD_CROSS_SECTION_BACKEND=boolean2 \
#     -DMANIFOLD_TEST=ON .
#   cmake --build build/daemon --target cross_section_fuzz -j$(nproc)
#   # Start the daemon:
#   nohup ./scripts/fuzz-cross-section-daemon.sh > /tmp/fuzz-daemon-stdout.log 2>&1 &
#
# ============================================================
# Lifecycle:
# ============================================================
#
#   start:   nohup ./scripts/fuzz-cross-section-daemon.sh > /tmp/fuzz-daemon-stdout.log 2>&1 &
#   status:  pgrep -af fuzz-cross-section-daemon
#            cat /tmp/fuzz-daemon.pid
#            tail -50 /tmp/fuzz-daemon/$(date +%Y-%m-%d)/master.log
#   pause:   touch /tmp/fuzz-daemon-pause   # daemon checks before each target
#   resume:  rm /tmp/fuzz-daemon-pause
#   stop:    kill $(cat /tmp/fuzz-daemon.pid)
#            (also: pkill -f fuzz-cross-section-daemon)
#

set -u

# Defaults, overridable via env.
DAEMON_CLONE="${FUZZ_DAEMON_CLONE:-$HOME/src/manifold-fuzz-daemon}"
TRACK_BRANCH="${FUZZ_DAEMON_BRANCH:-pr/cross-section-fuzz}"
BUILD_DIR="${FUZZ_DAEMON_BUILD_DIR:-build/daemon}"
CORPUS_DIR="${FUZZ_DAEMON_CORPUS_DIR:-build/fuzz-corpus}"
FUZZ_PER_SECONDS="${FUZZ_PER_SECONDS:-180}"
PAUSE_FILE="${FUZZ_DAEMON_PAUSE:-/tmp/fuzz-daemon-pause}"
PID_FILE="${FUZZ_DAEMON_PID:-/tmp/fuzz-daemon.pid}"
LOG_BASE="${FUZZ_DAEMON_LOG_BASE:-/tmp/fuzz-daemon}"

# Refuse to start if a daemon is already running.
if [[ -f "$PID_FILE" ]]; then
  old_pid=$(cat "$PID_FILE" 2>/dev/null || echo "")
  if [[ -n "$old_pid" ]] && kill -0 "$old_pid" 2>/dev/null; then
    echo "fuzz daemon already running (pid $old_pid). Stop it first or remove $PID_FILE." >&2
    exit 1
  fi
fi
echo $$ > "$PID_FILE"
trap 'rm -f "$PID_FILE"' EXIT

# Move into the daemon clone.
if [[ ! -d "$DAEMON_CLONE/.git" ]]; then
  echo "daemon clone not found at $DAEMON_CLONE. See setup at top of this script." >&2
  exit 1
fi
cd "$DAEMON_CLONE"

export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_container_overflow=0:halt_on_error=0:abort_on_error=0}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=0:print_stacktrace=0}"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

while_paused() {
  while [[ -e "$PAUSE_FILE" ]]; do
    log "paused: $PAUSE_FILE present, sleeping 60s"
    sleep 60
  done
}

log "fuzz daemon starting"
log "  clone:       $DAEMON_CLONE"
log "  branch:      $TRACK_BRANCH"
log "  build dir:   $BUILD_DIR"
log "  corpus dir:  $CORPUS_DIR (per-target subdirs)"
log "  per target:  ${FUZZ_PER_SECONDS}s"
log "  pause flag:  $PAUSE_FILE"
log "  pid:         $$ (recorded in $PID_FILE)"

cycle=0
while :; do
  cycle=$((cycle + 1))
  while_paused
  log "cycle $cycle: fetching origin"
  git fetch origin --quiet

  current=$(git rev-parse HEAD)
  upstream=$(git rev-parse "origin/$TRACK_BRANCH")
  if [[ "$current" != "$upstream" ]]; then
    log "cycle $cycle: source changed ($current -> $upstream); checking out"
    # Detached HEAD on the upstream tip. Daemon doesn't care about branch state.
    git checkout -q "$upstream"
  fi

  log "cycle $cycle: cmake --build (incremental)"
  if ! cmake --build "$BUILD_DIR" --target cross_section_fuzz -j"$(nproc)" \
       >>/tmp/fuzz-daemon-build.log 2>&1; then
    log "cycle $cycle: BUILD FAILED, sleeping 5m before retry"
    tail -20 /tmp/fuzz-daemon-build.log >&2
    sleep 300
    continue
  fi

  binary="$BUILD_DIR/test/cross_section_fuzz"
  if [[ ! -x "$binary" ]]; then
    log "cycle $cycle: binary missing at $binary after build, sleeping 5m"
    sleep 300
    continue
  fi

  # Refresh target list from the binary itself - no hardcoded array.
  mapfile -t targets < <(
    "$binary" --list_fuzz_tests 2>/dev/null \
      | grep -oE 'CrossSectionFuzz\.[A-Za-z_]+' | sort -u
  )
  log "cycle $cycle: ${#targets[@]} fuzz targets discovered"

  log_dir="$LOG_BASE/$(date +%Y-%m-%d)"
  mkdir -p "$log_dir"

  for t in "${targets[@]}"; do
    while_paused
    mkdir -p "$CORPUS_DIR/$t" "$log_dir/$t-artifacts"
    log "cycle $cycle: ==> $t"
    "$binary" --fuzz="$t" \
      "$CORPUS_DIR/$t" \
      -max_total_time="$FUZZ_PER_SECONDS" \
      -timeout=60 \
      -artifact_prefix="$log_dir/$t-artifacts/" \
      > "$log_dir/$t.log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
      log "cycle $cycle: <== $t exit=$rc (likely BUG FOUND; see $log_dir/$t.log)"
    else
      log "cycle $cycle: <== $t exit=0"
    fi
  done

  log "cycle $cycle: completed all ${#targets[@]} targets"
done
