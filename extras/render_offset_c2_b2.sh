#!/usr/bin/env bash
# Copyright 2026 The Manifold Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Diagnostic only (not shipped). Build render_offset_scenarios against both
# CrossSection backends, run each, and compose a self-contained HTML comparison.
# Run from a stack-current worktree so the boolean2 lane reflects the latest
# fixes. Output: <outdir>/offset_c2_b2.html (default extras/renders/).
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
outdir="${1:-${root}/extras/renders}"
mkdir -p "${outdir}"
jobs="$(nproc 2>/dev/null || echo 2)"

build_backend() { # $1 = clipper2|boolean2
  local backend="$1"
  local bd="${root}/build/cs-render-${backend}"
  if [[ ! -f "${bd}/CMakeCache.txt" ]]; then
    cmake -S "${root}" -B "${bd}" -DCMAKE_BUILD_TYPE=Release \
      -DMANIFOLD_CROSS_SECTION=ON \
      -DMANIFOLD_CROSS_SECTION_BACKEND="${backend}" \
      -DMANIFOLD_TEST=ON -DMANIFOLD_PAR=OFF -DBUILD_SHARED_LIBS=ON >/dev/null
  fi
  cmake --build "${bd}" --target renderOffsetScenarios -j"${jobs}" >/dev/null
  local bin
  bin="$(find "${bd}" -name renderOffsetScenarios -type f | head -1)"
  "${bin}" > "${outdir}/${backend}.json"
}

build_backend clipper2
build_backend boolean2
python3 "${root}/extras/compose_offset_renders.py" \
  "${outdir}/clipper2.json" "${outdir}/boolean2.json" \
  "${outdir}/offset_c2_b2.html"
echo "open ${outdir}/offset_c2_b2.html"
