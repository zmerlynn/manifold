// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Stage C of the 3D sweep-plane overlap-removal prototype.
// Design: docs/SweepPlane3D.md, section "Stage C - the sweep".

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>

#include "boolean2.h"
#include "overlap3.h"
#include "shared.h"

namespace manifold {

namespace {

// Section (y,z) coordinates of a 3D point.
static vec2 SectionYZ(vec3 p) { return {p.y, p.z}; }

// Compute the directed section segment for `face` at x = xMid.
// Direction: dot(p1-p0, yz(cross(+x, face.normal))) > 0
// (spec "ORIENTATION AND SIGN CONVENTIONS").
// Returns false if the face does not straddle xMid.
static bool ComputeSectionSegment(const CanonicalFace& face,
                                  const std::vector<MergedVert>& verts,
                                  double xMid, SectionFaceSegment& out) {
  const int v[3] = {face.verts.x, face.verts.y, face.verts.z};
  const vec3 p[3] = {verts[v[0]].pos, verts[v[1]].pos, verts[v[2]].pos};

  // Collect the two edge-crossings at xMid.
  vec2 pts[2];
  int found = 0;
  for (int i = 0; i < 3 && found < 2; ++i) {
    const int j = (i + 1) % 3;
    const double xi = p[i].x, xj = p[j].x;
    if (xi == xj) continue;                       // edge parallel to section
    if ((xi - xMid) * (xj - xMid) > 0) continue;  // same side
    const vec2 yz = Interpolate(p[i], p[j], xMid);
    if (found == 0 || yz.x != pts[0].x || yz.y != pts[0].y) {
      pts[found++] = yz;
    }
  }
  if (found < 2) return false;

  // Orient: dot(p1-p0, yz(cross(+x, face.normal))) > 0.
  // cross((1,0,0), (nx,ny,nz)) = (0,-nz,ny), yz-projection = (-nz, ny).
  const vec2 refDir = {-face.normal.z, face.normal.y};
  const vec2 d = pts[1] - pts[0];
  if (la::dot(d, refDir) < 0) std::swap(pts[0], pts[1]);

  out.faceId = -1;  // set by caller
  out.p0 = pts[0];
  out.p1 = pts[1];
  out.mult = face.mult;
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Stage C: x-criticals, per-slab sections, engine calls.
// ---------------------------------------------------------------------------

StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt) {
  const int nFaces = (int)arr.faces.size();

  // Collect x-criticals from all verts (stage-A + event + triple-point).
  std::vector<double> crits;
  crits.reserve(arr.verts.size());
  for (const auto& v : arr.verts) crits.push_back(v.pos.x);
  std::sort(crits.begin(), crits.end());
  crits.erase(std::unique(crits.begin(), crits.end()), crits.end());

  if (crits.empty()) return StageResult<std::vector<SlabResult>>::Ok({});

  // Exterior sentinel slabs (spec "BRACKETED by two exterior sentinel slabs").
  // Width = 2*eps > eps so they pass the slab-width gate and are built.
  // No face vertex lies in these x-ranges, so no face straddles them:
  // the engine produces empty capture => winding 0 everywhere (correct
  // exterior).
  const double xMin = crits.front();
  const double xMax = crits.back();
  const double sentW = 2.0 * eps;
  crits.insert(crits.begin(), xMin - sentW);
  crits.push_back(xMax + sentW);

  const int nSlabs = (int)crits.size() - 1;
  std::vector<SlabResult> slabs(nSlabs);

  for (int si = 0; si < nSlabs; ++si) {
    SlabResult& slab = slabs[si];
    slab.xLo = crits[si];
    slab.xHi = crits[si + 1];
    slab.xMid = (slab.xLo + slab.xHi) * 0.5;

    if (slab.xHi - slab.xLo <= eps) {
      slab.built = false;
      continue;
    }
    slab.built = true;

    // Collect section verts and edges for all straddling faces.
    std::vector<EdgeM> edges;
    std::vector<vec2> secVerts;
    std::map<std::pair<double, double>, int> vertIdx;

    auto getVid = [&](vec2 p) -> int {
      auto key = std::make_pair(p.x, p.y);
      auto it = vertIdx.find(key);
      if (it != vertIdx.end()) return it->second;
      const int id = (int)secVerts.size();
      secVerts.push_back(p);
      vertIdx.emplace(key, id);
      return id;
    };

    for (int fi = 0; fi < nFaces; ++fi) {
      const CanonicalFace& face = arr.faces[fi];
      const int vs[3] = {face.verts.x, face.verts.y, face.verts.z};
      const double xFaceMin =
          std::min({arr.verts[vs[0]].pos.x, arr.verts[vs[1]].pos.x,
                    arr.verts[vs[2]].pos.x});
      const double xFaceMax =
          std::max({arr.verts[vs[0]].pos.x, arr.verts[vs[1]].pos.x,
                    arr.verts[vs[2]].pos.x});
      if (xFaceMax <= slab.xMid || xFaceMin >= slab.xMid) continue;

      SectionFaceSegment seg;
      if (!ComputeSectionSegment(face, arr.verts, slab.xMid, seg)) continue;
      seg.faceId = fi;
      slab.segments.push_back(seg);

      const int v0 = getVid(seg.p0);
      const int v1 = getVid(seg.p1);
      if (v0 == v1) continue;
      edges.push_back({v0, v1, (int)face.mult, fi});
    }

    // Save raw section data for gate-2 test hooks.
    slab.sectionEdges = edges;
    slab.sectionVerts = secVerts;

    if (edges.empty()) continue;

    int conflictCount = 0;
    SweepWinding(edges, secVerts, WindRule::Add, &slab.pieces, &conflictCount);

    if (conflictCount > 0) {
      cnt.engineIdConflicts += conflictCount;
      return StageResult<std::vector<SlabResult>>::Fatal(
          FatalReason::EngineIdConflict, "engine id conflict in slab");
    }
  }

  return StageResult<std::vector<SlabResult>>::Ok(std::move(slabs));
}

}  // namespace manifold
