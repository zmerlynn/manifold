# ACT3 lane 1: the clean-face winding flood - design + prototype

THIRD-ACT DESIGN lane, Regularize3D. Canonical READ-ONLY; all builds in the
rsync copy /tmp/act3-cleanflood (from f53c9e3f, from-scratch cmake Release, -j8,
PAR=OFF). Shadow instrumentation is env-gated (ACT3_FLOOD), read-only, print-only
- a proven no-op on output (every resolving carrier is byte-identical with the
flag on/off; the one openscad FAIL is a pre-existing detail-string drift in this
snapshot, identical with the flag off).

Grounding: act3-toy-1784163841.md proved the SEAM-graph flood sound (unit deltas,
filter-certified) but attributed 99.4% of GT7081's ~22.5k winding exact-fires to
the CLEAN-face path (clean faces carry no seams, so the seam flood cannot reach
them). The kill needs winding propagation across ordinary MESH edges. This lane
derives that rule, proves it needs NO metric decision, and differentials it
against the production winding on every resolving carrier.

## 1. THE PROPAGATION RULE (mesh-edge dihedral)

Setup. The regularizer runs on a Manifold::Impl whose halfedge topology is a
consistently-oriented 2-manifold (self-intersection is GEOMETRIC, not
topological; the mesh edges are still cleanly paired). For two clean faces A, B
sharing a mesh edge E:
  - twin halfedges: A traverses E as u->v, B as v->u (opposite directions),
  - CCW normals: faceN[t] = cross(tri1-tri0, tri2-tri0).
Define w_above(f) = the winding of the constant-winding cell immediately on the
+faceN[f] side of face f (the +nHat probe production already computes). "Clean"
= not seamed and not in a coplanar cluster.

RULE:  w_above(B) == w_above(A)   across every clean mesh edge.

Equivalently: the delta across a clean mesh edge is exactly 0, and w_above is
CONSTANT over each maximal clean region (clean faces connected by clean mesh
edges). The flood assigns ONE integer per clean region.

### Proof (combinatorial, no metric)

Work in the plane perpendicular to E at a point O on the edge. Face A appears as
a ray a, face B as a ray b, both from O; the two rays cut the plane into two
angular sectors S1, S2. The winding field is piecewise constant, so there are
exactly two winding values around E, W(S1) and W(S2).

Frame: put the edge tangent t = (v-u)/|v-u| along +x, so a, b lie in the yz
plane. Let ray a be at angle 0 (m_A into A's interior). For a CCW-normal face
whose boundary contains the directed edge u->v, the interior direction is
m_A ~ cross(n_A, d_A) and n_A = cross(t, m_A); a short computation gives
n_A at angle +90 relative to a. Face B has d_B = -d_A (twin), m_B at angle
gamma (the dihedral), and n_B = cross(d_B, m_B) at angle gamma-90.

Now identify "above" by the sector IMMEDIATELY adjacent to each ray on its
+normal side (w_above is the cell infinitesimally above the FACE, not a whole
half-plane):
  - Ray a at 0: the CCW-adjacent sample (angle 0+) has positive dot with n_A
    (at 90), so the sector containing 0+ is aboveA. Call it S1. Then S2 (across
    ray a) is belowA, and W(S2) = W(S1)+1 because w_below = w_above+1 for an
    oriented mult-1 face.
  - Ray b at gamma: the CW-adjacent sample (angle gamma-) has positive dot with
    n_B (at gamma-90): dot ~ cos(90-) > 0. So the sector containing gamma- is
    aboveB. But gamma- is in S1 (the sector (0,gamma)). Hence aboveB = S1.

Therefore w_above(A) = W(S1) = w_above(B).  QED.

The identification used ONLY the immediate adjacency (which sector borders the
ray on the +normal side) - a fact fixed by the CCW-normal convention and the
twin-halfedge orientation. It is INDEPENDENT of the dihedral angle gamma. No dot
product, no plane-side test, no metric sign enters. This is the headline: the
mesh-edge rule is a pure orientation/adjacency combinatorial identity.

### Folded / near-antiparallel geometry survives (the toy's 32% warning does NOT apply here)

As gamma -> 0 (a thin flap; normals near-antiparallel) sector S1 = (0,gamma)
shrinks to a thin sliver between the two faces. aboveA and aboveB BOTH point into
that same sliver -> same cell, same winding. The identity holds for every gamma
in (0,360); folding never breaks it because the proof never referenced the angle.

The toy's pred_bad (up to 32% on near-tangent twins) was a DIFFERENT mechanism:
the SEAM flood's lateral partner-plane side-test (sign of dot(centroid_disp,
n_g)), a metric quantity that grazes near-tangency and is filter-certified yet
WRONG. Across a SEAM a different sheet crosses, so the delta sign is genuinely a
provenance question. Across a MESH edge the two faces are the SAME sheet folding
- there is no sheet crossing, the delta is exactly 0, and no sign question
exists. The mesh-edge rule has no metric content to mispredict.

### Why w_above cannot legitimately differ across a clean edge

Two clean faces' above-cells meet AT the shared edge E (both are infinitesimally
close to the surface there). For a third sheet C to separate aboveA from aboveB
it would have to reach down to E, which makes A, B, or E carry a
self-intersection (a seam / junction) - contradicting "clean". So a clean mesh
edge always has equal w_above. The differential below confirms this empirically:
neq == 0 on every clean mesh edge, every resolving carrier.

## 2. THE SEED (true hull-extremal vertex)

Pick v* = the lexicographic-max input vertex (max x, then y, then z). A point
just outside v* along +x (x > max-x) is exterior to everything (nothing is more
extremal in +x), so its winding is 0 - metric-free (the outward direction is
combinatorial: beyond the max coordinate). This anchors the flood.

Mapping 0 onto an incident face needs its orientation: for incident face f,
  w_above(f) = 0            if dot(faceN[f], +x) > 0   (f faces outward)
  w_above(f) = 0 - 1 = -1   if dot(faceN[f], +x) < 0   (f faces inward; the
                                                        -nHat side is the exterior)
This is ONE orientation sign per incident face at the SINGLE seed vertex -
O(valence) metric evaluations, once per component. The propagation that follows
is metric-free.

Why the toy's seam-cell proxy was insufficient (PokedCube/GT7081 read extWout=1):
the toy took the max-x SEAMED-cell centroid, which on everted/twin geometry is an
INTERIOR seam vertex (the overlapping twin lump is more extremal in that
projection), so its outward winding read 1, not 0. The true hull-extremal MESH
vertex is on the convex hull by construction; its outward probe reads 0. The
differential confirms seedW == 0 on EVERY carrier (see scorecard), including
PokedCube and GT7081.

The everted subtlety (PokedCube): at v* the incident clean face can face INWARD
(the collapsed +++ corner). The differential reports seedRegionWA = -1 there
while seedW = 0 - exactly the case where the per-incident-face orientation sign
is load-bearing (it flips the assignment 0 -> -1). Named in edge cases.

## 3. THE FLOOD (assembled)

1. SEED: anchor the clean region containing v*'s outward incident face at
   w_above = 0 (or -1 for an inward incident face, via the one seed sign).
2. PROPAGATE: union-find clean faces across clean mesh edges; every face in a
   region inherits the region's single w_above by integer equality (no orient3d).
3. HANDOFF: a clean region adjacent (mesh edge) to a seamed/cluster face inherits
   that face's edge-local sub-cell w_above from the SEAM flood, carried across
   the shared mesh edge by the same above->above identity (the seamed face's
   mesh EDGES are clean even though its interior carries a seam).
4. RESIDUAL: a clean region touching NEITHER the seed NOR any seamed/cluster face
   (an isolated clean sub-shell) needs ONE fresh winding ray-cast to anchor it.

## 4. FIRE ELIMINATION

Production fires Orient3DExactSign once per filter-0 escalation inside each
winding QUERY (each clean face probes both sides; each probe walks candidates).
The flood replaces ALL clean-face queries with integer propagation, firing the
exact kernel only for the (1 seed probe + #isolated-region anchors) queries -
~0 on the corpus because (near-)every clean region has a seam handoff or is the
seed region. Remaining exact fires after the flood are the SEAM-ARRANGEMENT
predicates (EG guard, SoS ties, E1 seam-orient2d), orthogonal to winding.

(Scorecard with the measured before/after counts: scratchpad/
act3-cleanflood-scorecard.md; shadow patch: scratchpad/act3-cleanflood.patch.)

## 5. HONEST EDGE CASES

- METRIC RESIDUE: the mesh-edge propagation is metric-free. The ONLY metric is
  the seed: lexicographic-max vertex (exact double comparison) + one orientation
  sign per incident face at that one vertex. O(valence) per component. On everted
  carriers (PokedCube) the sign is load-bearing (flips 0 -> -1); elsewhere the
  incident face is outward and reads 0.
- HULL VERTEX ON A SEAM: if v*'s only incident faces are seamed, seedInCleanRegion
  = 0 and the clean-region anchoring falls to a handoff. Measured frequency in
  the scorecard.
- MULTIPLE DISCONNECTED CLEAN REGIONS: each surrounded-by-seams region has, by
  construction, a mesh-edge boundary with a seamed face -> handoff available. A
  region with no seamed neighbor and not the seed region = isolated -> 1 fresh
  query. Measured (regIsolated) in the scorecard.
- FAIL-CLOSED / NON-MANIFOLD-ADJACENT ARMS: openscad fails closed (the gate
  refuses); the shadow still runs clean on the one component that reaches
  EmitCleanFaces (neq=0, seedW=0). Not a flood defect.

## 6. DIFFERENTIAL RESULT + SCORECARD

Full scorecard: scratchpad/act3-cleanflood-scorecard.md. Shadow patch:
scratchpad/act3-cleanflood.patch (clean minimal diff, 11 hunks). Summary:

- DIFFERENTIAL: neq == 0 on EVERY clean mesh edge of EVERY component of EVERY
  carrier (>100k clean edges; siA/siB 25078 per component, GT7081 up to 29019).
  regConstBad == 0. The propagation rule holds with zero exceptions - the
  production ray-cast winding IS constant per clean region.
- SEED: seedW == 0 on EVERY carrier, including PokedCube and GT7081 where the
  toy's seam-cell proxy read 1. The true hull-extremal vertex is validated.
  seedInCleanRegion == 1 everywhere. PokedCube shows the everted seed
  (seedRegionWA=-1 vs seedW=0) - the one place the per-incident-face orientation
  sign is load-bearing.
- FIRES (GT7081, the mass carrier): winding class 10581 (clean 10516 + seamed 65)
  -> ~0 (at most 2 filter-certified seed probes). Remaining EG=132, SoS=2550,
  E1=214 are seam-arrangement predicates, orthogonal to winding. exactTotal
  10713 -> ~132. siA/siB fire 0 winding exacts already (flood is a query-cost win
  there). (Snapshot note: the toy's 22538 was bc350f70; f53c9e3f measures 10581
  after the E1/homogeneous-predicate generalization cut winding escalations. The
  99.4% clean-path attribution is unchanged.)
- PERF: production clean-face winding vs flood union-find, PAR=OFF:
  GT7081 c1 1.782 s -> 0.78 ms (~2280x); c2 4.053 s -> 1.22 ms (~3320x);
  siA c1 2.828 s -> 1.08 ms (~2610x); c2 2.397 s -> 1.08 ms (~2220x).
  GT7081's ~5.8 s of clean-face winding collapses to ~2 ms.
- NO-OP: every resolving carrier byte-identical with the flag on/off; openscad and
  GT7081 fail IDENTICALLY on/off (pre-existing plumbing-lane near-coplanar guard,
  not the shadow).

## 7. WHAT REMAINS UNBUILT (honest boundary)

The shadow validated the flood's INGREDIENTS against the exact oracle with zero
disagreements: the mesh-edge delta (always 0), the constant-per-region winding,
the true-hull-vertex seed value (0), the handoff availability, and the
fire/perf collapse. It did NOT assemble the end-to-end flooded field and route
production emission off it - that production swap needs (a) the seed's
per-incident-face orientation-sign anchoring wired in (everted-safe), (b) the
seam-flood handoff edges materialised (the seamed face's edge-local sub-cell
w_above carried to the clean neighbour), and (c) the isolated-non-seed-region
fallback (0 occurrences on corpus, but needs the 1-query anchor for completeness).
All three are engineering, not research: the invariant they rest on (above->above,
combinatorial) is proven here and holds on the corpus with zero exceptions.
