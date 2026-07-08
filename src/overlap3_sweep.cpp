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

// Stage C of the 3D sweep-plane overlap-removal prototype: x-criticals,
// per-slab section build, engine calls with capture, and the point-winding
// query helper. Design: docs/SweepPlane3D.md.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "boolean2.h"
#include "overlap3.h"
#include "shared.h"

namespace manifold {

namespace {

// Section (y,z) coordinates of a 3D point.
static vec2 SectionYZ(vec3 p) { return {p.y, p.z}; }

// Directed section segment for face `fi` at x = xMid. The segment endpoints
// are computed by interpolating the face edges at xMid. The direction is chosen
// so that dot(p1-p0, yz(cross(+x, face.normal))) > 0.
// Returns false if the face does not straddle xMid (both endpoints same side).
static bool ComputeSectionSegment(const CanonicalFace& face,
                                  const std::vector<MergedVert>& verts,
                                  double xMid, SectionFaceSegment& seg) {
  // The three face edges. We need the two that straddle xMid.
  const int vs[3] = {face.verts.x, face.verts.y, face.verts.z};
  const vec3 p[3] = {verts[vs[0]].pos, verts[vs[1]].pos, verts[vs[2]].pos};
  const double x[3] = {p[0].x, p[1].x, p[2].x};

  // Find the two edges of the triangle that straddle xMid.
  // A face straddles xMid if at least one vert is < xMid and at least one is >
  // xMid.
  vec3 qA, qB;
  bool foundA = false, foundB = false;
  for (int i = 0; i < 3; ++i) {
    const int j = (i + 1) % 3;
    const double xi = x[i], xj = x[j];
    if ((xi <= xMid && xj >= xMid) || (xi >= xMid && xj <= xMid)) {
      if (xi == xj) continue;  // degenerate edge at xMid
      // Interpolate to get section point
      const vec2 yz = Interpolate(p[i], p[j], xMid);
      const vec3 q = {xMid, yz.x, yz.y};
      if (!foundA) {
        qA = q;
        foundA = true;
      } else if (!foundB && (q.y != qA.y || q.z != qA.z)) {
        qB = q;
        foundB = true;
      }
    }
  }
  if (!foundA || !foundB) return false;

  // Determine direction: dot(p1-p0, yz(cross(+x, face.normal))) > 0.
  // cross(+x, face.normal) = cross((1,0,0), face.normal)
  // = (0*nz - 0*ny, 0*nx - 1*nz, 1*ny - 0*nx)
  // = (0, -nz, ny) ... that is the cross product formula:
  // cross([1,0,0], [nx,ny,nz]) = [0*nz-0*ny, 0*nx-1*nz, 1*ny-0*nx] = [0,-nz,ny]
  // yz of that = (-nz, ny) in (y,z) section space.
  const vec2 refDir = {-face.normal.z,
                       face.normal.y};  // yz(cross(+x, face.normal))
  const vec2 dir = SectionYZ(qB) - SectionYZ(qA);
  if (la::dot(dir, refDir) < 0) {
    std::swap(qA, qB);
  }

  seg.faceId = -1;  // Will be set by caller
  seg.p0 = SectionYZ(qA);
  seg.p1 = SectionYZ(qB);
  seg.mult = face.mult;
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Stage C: build all slab results
// ---------------------------------------------------------------------------

StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt) {
  const int nFaces = (int)arr.faces.size();

  // Collect x-criticals: all vert x-coordinates (stage-A verts + event verts
  // + triple-point verts = all verts in arr.verts).
  std::vector<double> crits;
  crits.reserve(arr.verts.size());
  for (const auto& v : arr.verts) crits.push_back(v.pos.x);
  std::sort(crits.begin(), crits.end());
  crits.erase(std::unique(crits.begin(), crits.end()), crits.end());

  if (crits.empty()) {
    return StageResult<std::vector<SlabResult>>::Ok({});
  }

  // Add exterior sentinels (before min-x and after max-x).
  const double xMin = crits.front();
  const double xMax = crits.back();
  const double margin = std::max(eps * 10.0, (xMax - xMin) * 0.01 + eps);
  // Insert sentinels as boundary criticals.
  crits.insert(crits.begin(), xMin - margin);
  crits.push_back(xMax + margin);

  // Build one slab per consecutive pair of criticals.
  // Slab i: (crits[i], crits[i+1]).
  const int nSlabs = (int)crits.size() - 1;
  std::vector<SlabResult> slabs(nSlabs);

  for (int si = 0; si < nSlabs; ++si) {
    SlabResult& slab = slabs[si];
    slab.xLo = crits[si];
    slab.xHi = crits[si + 1];
    slab.xMid = (slab.xLo + slab.xHi) * 0.5;

    if (slab.xHi - slab.xLo <= eps) {
      slab.built = false;  // sub-eps slab: skip
      continue;
    }
    slab.built = true;

    // Find all faces that straddle xMid (strictly between xLo and xHi).
    // A face straddles if NOT all verts are on one side.
    std::vector<EdgeM> edges;
    std::vector<vec2> sectionVerts;
    std::map<std::pair<double, double>, int> vertId;

    auto getVertId = [&](vec2 p) -> int {
      const auto key = std::make_pair(p.x, p.y);
      auto it = vertId.find(key);
      if (it != vertId.end()) return it->second;
      const int id = (int)sectionVerts.size();
      sectionVerts.push_back(p);
      vertId.emplace(key, id);
      return id;
    };

    for (int fi = 0; fi < nFaces; ++fi) {
      const CanonicalFace& face = arr.faces[fi];
      const int vs[3] = {face.verts.x, face.verts.y, face.verts.z};
      const double xv[3] = {arr.verts[vs[0]].pos.x, arr.verts[vs[1]].pos.x,
                            arr.verts[vs[2]].pos.x};
      const double xFaceMin = std::min({xv[0], xv[1], xv[2]});
      const double xFaceMax = std::max({xv[0], xv[1], xv[2]});

      // Face straddles xMid if xFaceMin < xMid < xFaceMax
      // (strictly: not all on one side, and not all at xMid).
      if (xFaceMax <= slab.xMid || xFaceMin >= slab.xMid) continue;

      SectionFaceSegment seg;
      if (!ComputeSectionSegment(face, arr.verts, slab.xMid, seg)) continue;
      seg.faceId = fi;
      slab.segments.push_back(seg);

      // Add edge to the 2D engine input.
      const int v0id = getVertId(seg.p0);
      const int v1id = getVertId(seg.p1);
      if (v0id == v1id) continue;
      edges.push_back({v0id, v1id, (int)face.mult, fi});
    }

    if (edges.empty()) continue;  // no straddling faces

    // Run the 2D engine.
    int conflictCount = 0;
    SweepWinding(edges, sectionVerts, WindRule::Add, &slab.pieces,
                 &conflictCount);

    if (conflictCount > 0) {
      cnt.engineIdConflicts += conflictCount;
      return StageResult<std::vector<SlabResult>>::Fatal(
          FatalReason::EngineIdConflict, "engine id conflict in slab");
    }
  }

  return StageResult<std::vector<SlabResult>>::Ok(std::move(slabs));
}

}  // namespace manifold
