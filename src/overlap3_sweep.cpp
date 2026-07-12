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

// Slabs stage of the 3D sweep-native emission prototype: x-criticals,
// per-slab sections through the 2D engine, extension tracks.
// Design: docs/SweepEmit3D.md, SLABS.

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
// 3D edge-pair origins (va, vb) for each endpoint so the strips stage can
// extend to any x in the slab.
//
// Orientation: dot(p1 - p0, cross(+x, face.normal).yz()) > 0, computed as
// the 2D perp la::cross(1.0, normal.yz()).
//
// Returns false if the face does not straddle xMid (< 2 distinct crossings).
bool ComputeSectionSegment(SectionFaceSegment& segOut, FaceTrack& trackOut,
                           const CanonicalFace& face,
                           const std::vector<vec3>& verts, double xMid) {
  const int vi[3] = {face.verts.x, face.verts.y, face.verts.z};
  const vec3 p[3] = {verts[vi[0]], verts[vi[1]], verts[vi[2]]};

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
  const vec2 refDir = la::cross(1.0, face.normal.yz());
  const vec2 d = pts[1] - pts[0];
  if (la::dot(d, refDir) < 0) {
    std::swap(pts[0], pts[1]);
    std::swap(ePair[0], ePair[1]);
  }

  segOut.faceId = -1;  // set by caller
  segOut.p0 = pts[0];
  segOut.p1 = pts[1];
  segOut.mult = face.mult;

  trackOut.faceId = -1;  // set by caller
  trackOut.p0 = pts[0];
  trackOut.p1 = pts[1];
  trackOut.va0 = ePair[0][0];
  trackOut.vb0 = ePair[0][1];
  trackOut.va1 = ePair[1][0];
  trackOut.vb1 = ePair[1][1];
  return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// Slabs stage: x-criticals, per-slab sections, engine calls, face tracks.
// ---------------------------------------------------------------------------

StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt) {
  const int nFaces = static_cast<int>(arr.faces.size());

  // Collect x-criticals: all vert x's (canonical + seam endpoints) plus the
  // vertex-free criticals from the seams stage (degenerate contacts, seam-seam
  // crossings).
  std::vector<double> crits;
  crits.reserve(arr.verts.size() + arr.criticalXs.size());
  for (const auto& v : arr.verts) crits.push_back(v.x);
  for (double x : arr.criticalXs) crits.push_back(x);
  std::sort(crits.begin(), crits.end());
  crits.erase(std::unique(crits.begin(), crits.end()), crits.end());

  if (crits.size() < 2) return StageResult<std::vector<SlabResult>>::Ok({});

  // No sentinel slabs: ComputeCap handles null exterior (leftSlab=nullptr or
  // rightSlab=nullptr) as an empty region, which is the correct exterior limit
  // beyond the first/last critical (spec CAPS, "exterior limit is empty
  // region").
  const int nSlabs = static_cast<int>(crits.size()) - 1;
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
      const int id = static_cast<int>(secVerts.size());
      secVerts.push_back(p_yz);
      vertIdx.emplace(key, id);
      return id;
    };

    for (int fi = 0; fi < nFaces; ++fi) {
      const CanonicalFace& face = arr.faces[fi];
      const int vs[3] = {face.verts.x, face.verts.y, face.verts.z};
      const double xFaceMin = std::min(
          {arr.verts[vs[0]].x, arr.verts[vs[1]].x, arr.verts[vs[2]].x});
      const double xFaceMax = std::max(
          {arr.verts[vs[0]].x, arr.verts[vs[1]].x, arr.verts[vs[2]].x});
      if (xFaceMax <= slab.xMid || xFaceMin >= slab.xMid) continue;

      SectionFaceSegment seg;
      FaceTrack track;
      if (!ComputeSectionSegment(seg, track, face, arr.verts, slab.xMid))
        continue;
      seg.faceId = fi;
      track.faceId = fi;
      slab.segments.push_back(seg);
      slab.faceTracks.push_back(track);

      const int v0 = getVid(seg.p0);
      const int v1 = getVid(seg.p1);
      if (v0 == v1) continue;
      // Grouped faces seed the GROUP id (spec COPLANAR mechanism 2), so
      // coincident in-plane segments merge under one source id and
      // anti-oriented content cancels without a conflict.
      const int g = arr.face2Group.empty() ? -1 : arr.face2Group[fi];
      const int srcId = g >= 0 ? nFaces + g : fi;
      edges.push_back({v0, v1, static_cast<int>(face.mult), srcId});
    }

    slab.sectionEdges = edges;
    slab.sectionVerts = secVerts;

    // Populate seam tracks for class-ii endpoint extension (spec STRIPS).
    // For each seam whose x-range spans slab.xMid, record its (y,z) at xMid
    // and the 3D endpoints so caps/strips can use the seam track instead of
    // the face edge track for arrangement-constructed crossing vertices.
    // Membership is exact: seam endpoint x's are criticals and xMid is
    // strictly interior to its slab, so a seam either genuinely spans xMid or
    // has no crossing in this slab's section (this also excludes constant-x
    // seams, keeping InterpolateSafe's divisor nonzero).
    for (const auto& seam : arr.seams) {
      const vec3 vA = arr.verts[seam.vertId0];
      const vec3 vB = arr.verts[seam.vertId1];
      const double xLo3D = std::min(vA.x, vB.x);
      const double xHi3D = std::max(vA.x, vB.x);
      if (slab.xMid < xLo3D || slab.xMid > xHi3D) continue;
      const vec2 yzMid = InterpolateSafe(vA, vB, slab.xMid);
      slab.seamTracks.push_back({yzMid, vA, vB});
    }

    if (edges.empty()) continue;

    int conflictCount = 0;
    SweepWinding(edges, secVerts, WindRule::Add, &slab.pieces, &conflictCount);

    // Attribution-only source-id conflicts are BENIGN (spec [WALL-B]).  When
    // two ungrouped near-coplanar or near-degenerate faces coincide in a
    // section, SweepWinding flags srcId as conflicted (-1) but leaves the
    // geometry untouched: the retained boundary is emitted from net
    // multiplicity alone, and srcId has no downstream consumer in the 3D
    // pipeline.  So the conflict is a pure attribution diagnostic - count it
    // and continue rather than fail closed.
    cnt.engineIdConflicts += conflictCount;
  }

  // SubEpsFeature guard: any face whose x-extent has no coverage from a built
  // slab and whose area > perimeter * eps is a macroscopic feature in a merged
  // sub-eps critical run - fail closed.
  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];
    const vec3 p0 = arr.verts[face.verts.x];
    const vec3 p1 = arr.verts[face.verts.y];
    const vec3 p2 = arr.verts[face.verts.z];
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
