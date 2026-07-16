# THEORY lane: maximum construction depth / degree of the complete Regularize3D pipeline

Branch explore/sweep-plane-3d-v5, HEAD 7f8aec4e. READ-ONLY canonical except this
notebook. Sibling Codex lane works the same question independently. ASCII;
file:line cites the landed code; every unproven step is marked [UNPROVEN] or
[ASSUMED].

THE QUESTION (owner): for the COMPLETE exact regularization pipeline (valid
oriented 2-manifold self-overlapping component, double-coordinate vertices ->
boundary of {w_S>=1}), what is the maximum construction depth / polynomial degree
any decision can EVER require? Not the current broken path's empirical needs -
what the construction GRAMMAR can produce.

Answer in one line: DEPTH 1 (no arrangement point is ever built from a constructed
point); MAX PREDICATE DEGREE 9 in the plane-coefficient basis / 20 in the
input-coordinate basis (the constructed-point orient2d, instantiation (2)); the
decision layer that fixes the arrangement COMBINATORICS is degree 3 in inputs.


## 0. SCOPE AND GROUND TRUTH

Contract (docs/Regularize3D.md:21-49): the operator maps a valid oriented face
soup to the boundary of {p : w_S(p) >= 1}, per connectivity component. Gate
(overlap3.cpp:350,468): IsManifold && Is2Manifold, so EVERY edge is shared by
EXACTLY two faces - the 2-manifold guarantee this proof leans on. The kernel is
the two landed instantiations of ONE homogeneous form (overlap3.cpp:1068-1112
tripwire; homog-design-1784161597):
 (1) input-point orient3d, degree 3 (W==1, bit-identical to the historical
     adaptive-integer orient3d);
 (2) constructed-point orient2d of three plane-triple points, degree 9 in the
     plane coefficients (sos::HomogOrient2DExact, overlap3.cpp:938-952), on the
     SumSignN<8> accumulator (kAccumLimbs=320, overlap3.cpp:557).

I take the INPUT DOUBLES as ground truth (owner's framing: "exact inputs are the
DOUBLES"). Planes are then EXACT degree-<=3 functions of inputs, not primitives.


## 1. THE CONSTRUCTION GRAMMAR (exhaustive, with closure proof)

### 1a. Objects, generation by generation

LEVEL 0 - INPUT VERTICES. Exact doubles; degree 1. (overlap3.cpp:1727,
A.tri[t][k] = in.vertPos_[...].)

LEVEL 1 - FACE PLANES. A.faceN[t] = cross(t1-t0, t2-t0) (overlap3.cpp:1729-1730);
offset d = n.a with a = t0 (Intersect3Planes computes d_i = n_i . a_i,
overlap3.cpp:2299). Degrees IN INPUTS (this is the load-bearing decomposition):
  - normal n: cross of two edge-vectors = 2x2 determinants of degree-1 diffs =>
    DEGREE 2, per component.
  - offset d = n.a: degree 2 * degree 1 => DEGREE 3.
So a plane is (n[deg 2], d[deg 3]) - "degree <=3 in inputs", the 3 being the
offset. Every plane is fixed by exactly THREE input vertices (its triangle).
[In the plane-COEFFICIENT basis all of n,d are degree-1 generators; the campaign
counts degree there. The two bases are related by n=deg2, d=deg3 in inputs.]

LINES. Two kinds, both from INPUT generators only:
  (L-plane-pair) SEAM lines f/\g and TRANSVERSAL EDGE lines f/\(edge-neighbor).
    The edge line: the 2-manifold gate makes every edge shared by exactly two
    faces f, e2; if f,e2 are NOT coplanar the edge's supporting line = plane(f)
    /\ plane(e2), a plane pair. Verified in code: RecordSeams asserts "Every
    endpoint lands on the plane-i /\ plane-j intersection line, so it is exact in
    both faces' bases" (overlap3.cpp:1873).
  (L-point-pair) COPLANAR-CLUSTER EDGE lines. EDGE CASE the task flags: if the
    two faces sharing the edge ARE coplanar (a flat fold / a cap tessellation),
    plane(f)==plane(e2) and the pair does NOT determine a line. The edge line is
    then the LINE THROUGH ITS OWN TWO INPUT VERTICES - a point-pair line, a
    depth-1 object built from two input points (degree 1 each). It IS in the
    grammar, but see 1c: the coplanar FOLD discharges it before seam work.

POINTS (arrangement vertices). Exhaustive enumeration:
  (P0) INPUT VERTEX on a seam/edge/face. Trivial (degree 1). Predicates:
       input-vertex-on-edge (collinear/between, InputVertexStrictlyOnEdge via
       ExactOrient2DDrop, overlap3.cpp:2325) and input-vertex-on-face (orient3d,
       overlap3.cpp:1151) - both degree 3 in inputs.
  (P1) SEAM-SEAM-IN-FACE crossing = (f/\g) /\ (f/\h) = f/\g/\h = a plane TRIPLE
       point (Cramer, Intersect3Planes overlap3.cpp:2292 / CramerHPoint
       overlap3.cpp:903). Depth 1: three INPUT PLANES -> one point.
  (P2) SEAM-ENDPOINT on a TRANSVERSAL edge. Seam f/\g ends where it meets an edge
       of f; that edge = f/\e2 (plane pair). So the endpoint = f/\g/\e2 = again a
       plane TRIPLE {f,g,e2} (symwalk-1784185422 "seam endpoints = the
       {f,g,edge-neighbor} plane triple = CramerHPoint"; overlap3.cpp:1873). The
       ROUNDED POSITION is computed by the lower-degree SegPlanePoint (edge-
       segment /\ plane, overlap3.cpp:1648) - same point, cheaper formula - but
       its exact IDENTITY is the triple. Depth 1.
  (P3) SEAM-ENDPOINT on a COPLANAR-CLUSTER edge (the point-pair-line case). Seam
       f/\g meets a cluster edge whose neighbor is coplanar => no third plane =>
       the point = (f/\g) /\ line(u,v) = a plane-pair line crossing a point-pair
       line, both in plane f. Exact identity = SegPlanePoint(u,v, plane g): the
       edge-segment (two input points) pierces plane g. A DIFFERENT depth-1 object
       (2 input points + 1 plane), NOT a pure plane triple. Degree: see 2c
       (strictly LOWER than P1/P2). This is the "coplanar/transversal
       entanglement" (docs A1, overlap3.cpp:3558-3559) and the openscad E1
       terminal.

### 1b. Closure theorem (depth 1)

CLAIM: no construction step ever consumes a CONSTRUCTED point to build a new line
or point. Hence every arrangement point is at most ONE intersection-construction
away from input-derived planes/lines: DEPTH 1.

PROOF (by exhaustion over the generators of every object above):
  - Every LINE (L-plane-pair, L-point-pair) is generated by input planes or input
    vertices ONLY. No line is "through two constructed points."
  - Every POINT (P0-P3) is generated by (three input planes) or (two input planes
    + a point-pair line) or (a plane + an input edge-segment) or (an input
    vertex). No point is "intersection of two constructed points/lines."
  - EMISSION does not break this: SplitTouchingSheets (overlap3.cpp:77) operates
    on already-emitted rounded geometry by connectivity, minting no exact
    construction; the winding probes (centroid +/- eps*nHat, overlap3.cpp:3076,
    3350,3451) are PROBE positions, not arrangement vertices, and the crossing
    predicate against them is orient3d on input triangles (WindCrossTri
    overlap3.cpp:1425-1441) - degree 3; the once-only construction rule
    (docs:156-158) mints each point ONCE and references it, never re-intersecting
    a constructed point. []

### 1c. Why the campaign's "nested constructions" are NOT grammar objects

The symwalk/exjunct/stitch lanes observed openscad residue that LOOKS like
depth-2 (skew-line projected crossings, phantom seam-fragments). These are
ARTIFACTS of the double-precision eps-weld, not grammar objects:
  - exjunct-1784189568 measured that the residual crossings are between sub-edges
    "whose carrier LINES share NO plane (3D-SKEW)", so "the crossing is the
    2D-projection near-crossing of two skew lines with no common plane - NOT a
    plane-triple point." In EXACT 3D the skew lines DO NOT cross; the apparent
    crossing exists only after the junction registry eps-snaps a foreign seam
    endpoint within eps of a face's edge (BuildJunctionRegistry, f4-junction).
  - exjunct proved these phantoms have NO exact plane-triple identity within the
    weld: "EVERY exact plane-triple re-derivation ... moves the grazing cross-face
    vertices by 11-1700 eps" (one to four orders beyond eps). A genuine
    arrangement vertex would move 0. So they are not depth-2 points; they are
    non-points manufactured by rounding.
  - stitch-1784195475 pinned the same residue at the boundary-emission level: a
    "crowded near-tangent triple-point cluster ... on DIFFERENT plane triples"
    (all genuine P1/P2 triples, spread ~1e-6..1e-5 >> eps) that the double
    per-face arrangement cannot subdivide two-sided-consistently. The vertices are
    depth-1 triples; the FAILURE is representability of their ORDER in doubles, not
    a deeper construction.
CONCLUSION: the grammar is depth 1; the openscad wall is the double-precision
2D-arrangement rounding of depth-1 objects, exactly the E1 terminal
(docs:580-646).


## 2. DEGREES (both bases; the max over all needed predicates)

Convention: "degree" = total polynomial degree of the (homogeneous) form whose
SIGN the predicate evaluates. Basis (a) = plane coefficients (n,d as generators);
basis (b) = input coordinates (n=deg2, d=deg3 substituted).

### 2a. The decision layer (what the arrangement combinatorially IS)

Every existence/relation decision restructures to orient3d of INPUT points
(docs:184-185 "every combinatorial decision restructures to a level-0 input
predicate"). Callers and their degree in inputs:
  - face-pair straddle / edge-plane crossing / edge-edge z-order: EdgePiercesTri
    (overlap3.cpp:1149), all Orient3DFilterSign => DEGREE 3.
  - vertex-on-face, coplanarity: FacesFilterCoplanar (overlap3.cpp:1228) =>
    DEGREE 3.
  - winding-crossing delta (the coupled integer w_S): WindCrossTri
    (overlap3.cpp:1425) orient3d, delta is +-1 integer (docs:130-132 "FP-safe by
    construction") => DEGREE 3.
  - radial order around a >2-sheet edge (UNBUILT, kernel-tripwire-gated): the
    determinant rule is "level-0 FREE from input-plane normals + the emitted fwd
    bit" (docs:597-599). [ASSUMED degree 3 - it orders half-planes about a line by
    input-normal orientation; no corpus carrier forces it, so this is the spec's
    claim, not re-derived here.]
DECISION-LAYER MAX = DEGREE 3 in inputs. Exact w.r.t. true inputs.

### 2b. The near-tangent emission layer (instantiation (2), the MAX)

The exact per-face 2D arrangement over CONSTRUCTED crossing points (the E1 spec,
docs:600-602). Predicate = orient2d of three plane-triple points P0,P1,P2 =
sign(det[[A0,B0,W0],[A1,B1,W1],[A2,B2,W2]]) * sign(W0 W1 W2)
(HomogOrient2DExact, overlap3.cpp:938-952).

Triple point homogeneous coords (CramerHPoint, overlap3.cpp:903-925):
  W = nF.(ng x nh):       basis(a) DEGREE 3; basis(b) 2+2+2 = DEGREE 6.
  X,Y,Z = dF*(ng x nh) + ...: basis(a) DEGREE 3; basis(b) 3 + (2+2) = DEGREE 7.
orient2d determinant (each Leibniz term picks one entry per row AND per column =
one A[deg7], one B[deg7], one W[deg6]):
  basis(a): 3+3+3 = DEGREE 9.    basis(b): 7+7+6 = DEGREE 20 (tight, HOMOGENEOUS).
W-product sign correction W0 W1 W2:
  basis(a): DEGREE 9.            basis(b): 6*3 = DEGREE 18.
Single-polynomial decision det*(W0 W1 W2) (if one refuses to factor):
  basis(a): DEGREE 18.          basis(b): 20+18 = DEGREE 38.
LOOSE bound (treat every plane coeff as degree 3 in inputs, ignoring n=deg2): det
= 9*3 = 27, single-poly 54. The TIGHT numbers (20 / 38) are the honest ones.

The KERNEL factors: it evaluates det (deg 9 / 20) and each W (deg 3 / 6)
separately, so the widest single polynomial it ever sums is DEGREE 9 (plane) /
DEGREE 20 (input). kAccumLimbs=320 was sized for exactly this (9*53=477 mag bits +
worst-case exponent spread, overlap3.cpp:545-557; homog-design measured 12 limbs
active on mesh data).

### 2c. The coplanar/transversal point (P3), priced

P3 = SegPlanePoint(u,v, plane g) = edge(u,v) /\ plane g. Homogeneous form:
  W = ng.(v-u):           basis(b) 2+1 = DEGREE 3.
  coord = u*(ng.(v-u)) + (dg - ng.u)*(v-u): basis(b) DEGREE 4.
orient2d of three P3-type points: term = 4+4+3 = DEGREE 11 in inputs (< 20).
Mixed with P1/P2 triples it is bounded by the triple case (deg 20). So P3 does NOT
raise the ceiling; it is the "one or two SMALL new instantiations for
point-pair-line predicates" the task anticipated. It is a new GENERATOR MIX (input
points enter the HPoint, not just plane coeffs) but NOT a new degree - it fits
UNDER instantiation (2)'s accumulator and needs only its own (lower) filter
constant. Price: ~30-60 LOC (a point-pair HPoint/EHPoint constructor reusing
HomogOrient2DExact/HomogOrient2DFilter verbatim) + one filter constant. [UNPROVEN:
the LOC and that a differential oracle certifies the degree-11 filter as cleanly
as homog-design's degree-9 - not built or validated in this lane.]

### 2d. Basis verdict (where the kernel SHOULD evaluate)

The DECISION LAYER (2a) must certify in the INPUT basis (degree 3) - it decides
what the arrangement IS, and it must agree with the true input geometry
(coplanarity, existence). It does: FacesFilterCoplanar / EdgePiercesTri are
orient3d on input doubles. This is exact w.r.t. inputs.

The EMISSION LAYER (2b) currently evaluates in the PLANE-COEFFICIENT basis:
CramerHPoint takes ROUNDED double plane coefficients (A.faceN rounds at
overlap3.cpp:1730), so instantiation (2) is exact RELATIVE TO ROUNDED PLANES, not
to true inputs. This is a PLANE-ROUNDING GAP. For a FULLY exact path the ground
truth is the input doubles, so the orient2d should either carry planes EXACTLY
(exact-rational n[deg2]/d[deg3] polys, evaluating the determinant at degree 20 in
inputs) or accept the rounding as a deliberate quantization. The gap only bites in
the sub-eps ordering tail (the openscad E1 residue); the COMBINATORICS are pinned
by the degree-3 input layer, which is why the corpus off openscad resolves
byte-clean with the degree-9-in-rounded-planes form. [UNPROVEN: that carrying
planes exactly closes openscad - e1-plumb showed even the exact-in-rounded-planes
extent does not, and the terminal is the exact-rational 2D arrangement replacing
RemoveOverlaps2D; whether input-exact planes vs rounded-exact planes changes the
E1 outcome is not measured here.]


## 3. THE THEOREM

Let the input be a valid oriented 2-manifold connectivity component with
double-coordinate vertices (the gate, overlap3.cpp:468). Let the pipeline be the
COMPLETE exact regularization: exact-coplanar fold, transversal seam arrangement,
coupled integer winding, and the exact {w_S>=1} boundary emission (the E1 exact
2D arrangement in place of the double-precision RemoveOverlaps2D).

THEOREM (construction depth).
  Every arrangement vertex is one of {input vertex; plane triple f/\g/\h; plane-
  pair /\ point-pair-line point}. No construction consumes a constructed point.
  MAX CONSTRUCTION DEPTH = 1.

THEOREM (polynomial degree).
  (i) The DECISION layer (seam existence, straddle, edge-edge order, vertex-on-
      face, coplanarity, winding-crossing, radial order) is orient3d of input
      points: DEGREE 3 in the input-coordinate basis. Exact w.r.t. inputs.
  (ii) The EMISSION layer's widest predicate is the constructed-point orient2d
       (three plane-triple points): DEGREE 9 in the plane-coefficient basis,
       DEGREE 20 in the input-coordinate basis (tight/homogeneous), with a
       DEGREE-3 (plane) / DEGREE-6 (input) W-sign correction. The fully-cleared
       single-polynomial form is degree 18 (plane) / 38 (input); the kernel
       factors and never forms it.
  MAX PREDICATE DEGREE = 9 (plane basis) / 20 (input basis).

PROOF SKETCH.
  Depth: section 1b (exhaustion over generators; no point-of-points anywhere,
  including emission). Degree (i): section 2a + docs:184-185. Degree (ii): section
  2b (Cramer coords deg 6/7 in inputs, determinant deg 20, homog-design's deg-9
  plane-basis count). The point-pair case P3 is deg <= 11 in inputs (2c), strictly
  under the ceiling. The campaign's apparent depth-2 objects are rounding phantoms
  with no exact identity within eps (1c; exjunct's 11-1700 eps measurement), hence
  not grammar objects. []

GRAMMAR-CLOSURE VERDICT: depth 1 SUFFICES. The landed two-instantiation kernel
(orient3d deg 3 + constructed-point orient2d deg 9) plus ONE small point-pair-line
instantiation (P3, deg <= 11 in inputs, reusing the same homogeneous form under
the same accumulator) is a COMPLETE exact-predicate basis for the whole pipeline.
No third DEGREE is ever forced; the openscad terminal is the exact-rational 2D
ARRANGEMENT ENGINE consuming these predicates (~600-1300 LOC, docs:638-641), not a
deeper predicate.


## 4. IMPLICATIONS

### 4a. Kernel inventory the closure build actually needs

  # | predicate                          | plane-deg | input-deg | callers
  --|------------------------------------|-----------|-----------|--------------
  1 | orient3d, input points (inst. 1)   |     3     |     3     | MANY (straddle,
    |                                    |           |           | z-order, v-on-
    |                                    |           |           | face, coplanar,
    |                                    |           |           | winding, SoS)
  2 | orient2d, 3 plane-triple pts (i.2) |     9     |    20     | ExactSeamsCross
    | + W-sign correction                |     3     |     6     | / HPointStrict-
    |                                    |           |           | lyInTri (E1)
  3 | orient2d, point-pair/\plane pts(P3)|   <=6     |   <=11    | coplanar/trans.
    |   [NEW, small; owner-gated]        |           |           | entanglement
  --|------------------------------------|-----------|-----------|--------------
Count: TWO landed + ONE small new (same form, lower degree, same accumulator). The
SoS convention (Edelsbrunner-Mucke, overlap3.cpp:758) rides #1 only (input-vertex
indices; the constructed sites use level-0 incidence, not perturbation -
homog-design "SoS stays INPUT-POINT-SCOPED"). Constructions (positions, NOT
decisions, tripwire-free): SegPlanePoint (1648), Intersect3Planes (2292),
SegLineIntersect2D (2308) - all depth-1 rational maps of inputs.

RECOMMENDED EVALUATION BASIS: decision layer #1 in the INPUT basis (degree 3,
exact w.r.t. the ground-truth doubles). Emission layer #2/#3 in the plane-
coefficient basis if planes are accepted as rounded (the landed choice, exact-
relative-to-rounded-planes), OR the input basis (degree 20) with exact-rational
planes for a strictly input-exact path (closes the plane-rounding gap; wider
accumulator). [UNPROVEN which is needed to cross the E1 wall - see 2d.]

### 4b. Coplanar-fold edge-case treatment (verified in code)

The fold DISCHARGES the coplanar-adjacent-edge line ambiguity, but only for the
PURE-coplanar case:
  - Ordering: EmitComponentBoundary runs FoldCoplanarClusters FIRST, then the
    transversal seam emit (overlap3.cpp:3560, 3586-3587). Comment:
    "Exact-coplanar clusters first (transversal seams on their faces are the
    entanglement decline)" (3558-3559).
  - Detection: exact coplanarity via orient3d on inputs (FacesFilterCoplanar
    1228, degree 3). Same-cluster pairs are SKIPPED by the seam machinery
    (overlap3.cpp:1868), so a coplanar-adjacent edge NEVER enters the seam
    machinery as an ill-defined plane-pair line - it is resolved as a 2D point-
    pair-line edge inside the fold's within-plane RemoveOverlaps2D arrangement.
    => the L-point-pair edge line (1a) is discharged; no plane-pair ambiguity
    survives for coplanar-coplanar edges.
  - NOT discharged: a TRANSVERSAL seam crossing INTO a coplanar cluster. Its
    endpoint on a cluster edge is the P3 point (plane-pair /\ point-pair-line),
    which the fold's DOUBLE-precision RemoveOverlaps2D cannot place exactly. This
    is the coplanar/transversal entanglement (docs A1) and the openscad E1
    terminal (exjunct: no exact plane-triple within eps). The exact fix is the
    #3 predicate driving an exact 2D arrangement in the cluster plane - the SAME
    engine E1 needs, not a new degree.


## 5. UNPROVEN / ASSUMED LEDGER (honesty rail)

  [A] Radial-order predicate degree = 3: the spec's claim (docs:597-599), not
      re-derived; no corpus carrier forces the >2-sheet branch so it is unbuilt.
  [B] P3 (point-pair-line) filter cleanliness and ~30-60 LOC: priced by analogy
      to homog-design's degree-9 validation, NOT built or oracle-tested here.
  [C] Input-exact planes (degree 20) vs rounded-exact planes (degree 9) closing
      openscad: NOT measured; e1-plumb shows the exact-in-rounded-planes extent
      does not close it, and the terminal is the exact-rational 2D arrangement
      engine regardless of plane basis.
  [D] Degree bounds 20/38 are for the CramerHPoint homogeneous representation
      (unreduced); a spurious common factor could lower the true irreducible
      degree - not checked, but the kernel evaluates the unreduced form, so 20 is
      the honest degree of the polynomial the kernel sums.
  [E] "No point-of-points in emission" (1b) checked by reading SplitTouchingSheets
      / winding-probe / once-only construction; not a mechanized proof over every
      emit path.


## VERDICT

Max construction DEPTH = 1 (no arrangement point is built from a constructed
point; the campaign's skew-crossing "nested constructions" are eps-weld phantoms
with no exact identity within eps, per exjunct's 11-1700 eps measurement). Max
predicate DEGREE = 9 (plane-coefficient basis) / 20 (input-coordinate basis,
tight) for the constructed-point orient2d, + a degree-3/6 W-sign; the decision
layer that fixes the arrangement combinatorics is degree 3 in inputs. Kernel
inventory for a complete exact path: orient3d deg-3 (inst. 1, many callers) +
constructed-point orient2d deg-9 (inst. 2, E1) + ONE small point-pair-line orient2d
(deg <= 11 in inputs, same homogeneous form, same accumulator) for the
coplanar/transversal entanglement - THREE predicates, ONE form, no third degree
ever forced. The coplanar fold discharges the coplanar-coplanar edge-line
ambiguity (runs first, same-cluster seam-skip at overlap3.cpp:1868) but NOT the
transversal-seam-into-cluster endpoint (P3), which is the openscad E1 wall and
needs the exact 2D arrangement engine, not a deeper predicate.
