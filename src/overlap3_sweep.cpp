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

// Stage C' of the 3D sweep-native emission prototype.
// Design: docs/SweepEmit3D.md, Stage C'.

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

// Compute the directed section segment for `face` at x = xMid, storing the
// 3D edge-pair origins (va, vb) for each endpoint so stage D' can extend to
// any x in the slab.
//
// Orientation: dot(p1 - p0, yz(cross(+x, face.normal))) > 0.
// cross((1,0,0),(nx,ny,nz)) = (0, -nz, ny), yz-projection = (-nz, ny).
//
// Returns false if the face does not straddle xMid (< 2 distinct crossings).
static bool ComputeSectionSegment(const CanonicalFace& face,
                                  const std::vector<MergedVert>& verts,
                                  double xMid, SectionFaceSegment& segOut,
                                  FaceTrack* trackOut) {
  const int vi[3] = {face.verts.x, face.verts.y, face.verts.z};
  const vec3 p[3] = {verts[vi[0]].pos, verts[vi[1]].pos, verts[vi[2]].pos};

  vec2 pts[2];
  vec3 ePair[2][2];  // ePair[k] = {va, vb} for pts[k]
  int found = 0;

  for (int i = 0; i < 3 && found < 2; ++i) {
    const int j = (i + 1) % 3;
    const double xi = p[i].x, xj = p[j].x;
    if (xi == xj) continue;                       // edge parallel to section
    if ((xi - xMid) * (xj - xMid) > 0) continue;  // same side of xMid
    const vec2 yz = Interpolate(p[i], p[j], xMid);
    if (found == 0 || yz.x != pts[0].x || yz.y != pts[0].y) {
      pts[found] = yz;
      ePair[found][0] = p[i];
      ePair[found][1] = p[j];
      found++;
    }
  }
  if (found < 2) return false;

  // Orient: dot(p1-p0, refDir) > 0.
  const vec2 refDir = {-face.normal.z, face.normal.y};
  const vec2 d = pts[1] - pts[0];
  if (la::dot(d, refDir) < 0) {
    std::swap(pts[0], pts[1]);
    std::swap(ePair[0], ePair[1]);
  }

  segOut.faceId = -1;  // set by caller
  segOut.p0 = pts[0];
  segOut.p1 = pts[1];
  segOut.mult = face.mult;

  if (trackOut) {
    trackOut->faceId = -1;  // set by caller
    trackOut->p0 = pts[0];
    trackOut->p1 = pts[1];
    trackOut->va0 = ePair[0][0];
    trackOut->vb0 = ePair[0][1];
    trackOut->va1 = ePair[1][0];
    trackOut->vb1 = ePair[1][1];
  }
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Stage C': x-criticals, per-slab sections, engine calls, face tracks.
// ---------------------------------------------------------------------------

StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt) {
  const int nFaces = (int)arr.faces.size();

  // Collect x-criticals: all vert x's (stage-A + seam endpoints) plus the
  // vertex-free criticals from stage B' (degenerate contacts, seam-seam
  // crossings).
  std::vector<double> crits;
  crits.reserve(arr.verts.size() + arr.criticalXs.size());
  for (const auto& v : arr.verts) crits.push_back(v.pos.x);
  for (double x : arr.criticalXs) crits.push_back(x);
  std::sort(crits.begin(), crits.end());
  crits.erase(std::unique(crits.begin(), crits.end()), crits.end());

  if (crits.size() < 2) return StageResult<std::vector<SlabResult>>::Ok({});

  // No sentinel slabs: ComputeCap handles null exterior (leftSlab=nullptr or
  // rightSlab=nullptr) as an empty region, which is the correct exterior limit
  // beyond the first/last critical (spec E', "exterior limit is empty region").
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

    std::vector<EdgeM> edges;
    std::vector<vec2> secVerts;
    std::map<std::pair<double, double>, int> vertIdx;

    auto getVid = [&](vec2 p_yz) -> int {
      auto key = std::make_pair(p_yz.x, p_yz.y);
      auto it = vertIdx.find(key);
      if (it != vertIdx.end()) return it->second;
      const int id = (int)secVerts.size();
      secVerts.push_back(p_yz);
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
      FaceTrack track;
      if (!ComputeSectionSegment(face, arr.verts, slab.xMid, seg, &track))
        continue;
      seg.faceId = fi;
      track.faceId = fi;
      slab.segments.push_back(seg);
      slab.faceTracks.push_back(track);

      const int v0 = getVid(seg.p0);
      const int v1 = getVid(seg.p1);
      if (v0 == v1) continue;
      edges.push_back({v0, v1, (int)face.mult, fi});
    }

    slab.sectionEdges = edges;
    slab.sectionVerts = secVerts;

    // Populate seam tracks for class-ii endpoint extension (spec D').
    // For each seam whose x-range spans slab.xMid, record its (y,z) at xMid
    // and the 3D endpoints so caps/strips can use the seam track instead of
    // the face edge track for arrangement-constructed crossing vertices.
    for (const auto& seam : arr.seams) {
      const vec3 vA = arr.verts[seam.vertId0].pos;
      const vec3 vB = arr.verts[seam.vertId1].pos;
      const double xLo3D = std::min(vA.x, vB.x);
      const double xHi3D = std::max(vA.x, vB.x);
      if (slab.xMid < xLo3D - eps || slab.xMid > xHi3D + eps) continue;
      const vec2 yzMid = InterpolateSafe(vA, vB, slab.xMid);
      slab.seamTracks.push_back({yzMid, vA, vB});
    }

    if (edges.empty()) continue;

    int conflictCount = 0;
    SweepWinding(edges, secVerts, WindRule::Add, &slab.pieces, &conflictCount);

    if (conflictCount > 0) {
      cnt.engineIdConflicts += conflictCount;
      return StageResult<std::vector<SlabResult>>::Fatal(
          FatalReason::EngineIdConflict, "engine id conflict in slab");
    }
  }

  // SubEpsFeature guard: any face whose x-extent has no coverage from a built
  // slab and whose area > perimeter * eps is a macroscopic feature in a merged
  // sub-eps critical run - fail closed.
  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];
    const vec3 p0 = arr.verts[face.verts.x].pos;
    const vec3 p1 = arr.verts[face.verts.y].pos;
    const vec3 p2 = arr.verts[face.verts.z].pos;
    const double xFaceMin = std::min({p0.x, p1.x, p2.x});
    const double xFaceMax = std::max({p0.x, p1.x, p2.x});
    if (xFaceMin >= xFaceMax) continue;  // axis-parallel face, OK

    bool hasCoverage = false;
    for (const auto& slab : slabs) {
      if (!slab.built) continue;
      // Built slab covers any part of the face's x-range.
      if (slab.xLo < xFaceMax && slab.xHi > xFaceMin) {
        hasCoverage = true;
        break;
      }
    }
    if (hasCoverage) continue;

    const double area = 0.5 * la::length(la::cross(p1 - p0, p2 - p0));
    const double perim =
        la::length(p1 - p0) + la::length(p2 - p1) + la::length(p0 - p2);
    if (area > perim * eps) {
      return StageResult<std::vector<SlabResult>>::Fatal(
          FatalReason::SubEpsFeature,
          "face in unbuilt x-range with area > perimeter*eps");
    }
  }

  return StageResult<std::vector<SlabResult>>::Ok(std::move(slabs));
}

}  // namespace manifold
