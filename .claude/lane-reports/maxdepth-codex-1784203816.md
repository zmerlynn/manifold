# THEORY lane (2 of 2, Codex operator): max construction depth / degree, Regularize3D

Codex 0.139.0, exec --dangerously-bypass-approvals-and-sandbox over a read-only
rsync of src/ + docs/ + lane-reports (f4-junction/symwalk/exjunct/graze/stitch/
e1-plumb/homog-design) at /tmp/codex-maxdepth/tree.  Raw log:
/tmp/codex-maxdepth/codex-log.txt.  Distilled raw:
.claude/lane-reports/scratchpad/codex-maxdepth-findings.txt.  I triaged Codex's
derivation MYSELF against overlap3.cpp at file:line and hand-checked its degree
arithmetic on 4 entries.  ASCII; READ-ONLY canonical (nothing built, nothing pushed).

## CODEX'S THEOREM (verbatim core)

- Construction grammar is CLOSED AT DEPTH 1: every real arrangement object is a
  polynomial in the RAW INPUTS (input vertices + face planes), never in other
  constructed objects.  Seam endpoints, plane triples, coplanar-edge crossings -
  all depth 1.
- MAX EXACT TOPOLOGICAL PREDICATE (raw input basis) = degree 20: the general
  homogeneous orient2d of three plane-triple points.  Each plane-triple point has
  homogeneous coords X,Y,Z degree 7, W degree 6; the 2D determinant term is
  7+7+6 = 20.
- The code/doc "degree 9" is correct IN THE PLANE-COEFFICIENT GENERATOR BASIS
  (Cramer coords degree 3 in plane coeffs; 3x3 orient2d = 9).  Same primitive,
  raw basis = 20.
- SHIPPED production predicate (ExactSeamsCross via HPointStrictlyInTri) maxes at
  degree 8: ONE constructed triple point vs TWO input triangle vertices.
- Skew-line "crossings of chords whose endpoints are constructed" (symwalk /
  exjunct / graze) are ARTIFACTS of eps-injection / rounded-foot placement, NOT
  grammar objects.  If exactized, an eps distance-collapse would be degree 26
  (numerator 7+6=13, squared) - but the pipeline keeps that decision in doubles.

## MY TRIAGE (verified against overlap3.cpp; hand-checked arithmetic)

VERDICT: Codex is CORRECT on every load-bearing claim.  No errors found in its
degree arithmetic or its adjudication.  Refinements, not corrections, below.

Hand-checked degrees (all agree with Codex):
- plane = (normal deg 2 = cross of two deg-1 edge vectors; offset d = n.p0 deg 3).
  code: faceN = cross(b-a,c-a) at overlap3.cpp:1729-1730; Intersect3Planes uses
  d_i = dot(n_i, a_i) at 2299.
- plane triple (Cramer): W = nF.(nGxnH) = 2.(4) = deg 6; X = dF.(nGxnH)_x =
  3.4 = deg 7.  -> (X,Y,Z,W)=(7,7,7,6).  code CramerHPoint at 903-926.
- orient2d of 3 triples: det columns {A:7,B:7,W:6} -> 20.  code HomogOrient2DExact
  938-952, accumulated on SumSignN<8> (the "degree-9" NMag=8 form, 512 mag-bits).
- in-triangle (triple vs 2 input verts): rows [(7,7,6),(1,1,0),(1,1,0)] ->
  max term 8.  This is the SHIPPED form: HPointStrictlyInTri (2360-2374) calls
  HomogOrient2DFilter/Exact with ETrivialHPoint(a), ETrivialHPoint(b) (input,
  W==1) and the constructed eT.  So production NEVER feeds three constructed
  points; degree-20 is the FORM'S CAPACITY, degree-8 is the shipped decision.
- seam endpoint two representations: as plane triple {f,g,neighbor} = (7,6);
  as input-edge INT plane = (4,3) (SegPlanePoint 1648).  Codex's 4/3 is the
  minimal rep; reconciles the graze report's "degree-6 constructed-point-on-plane"
  (= seam-endpoint(4/3) . plane(2/3) = 6, raw basis).

Adjudication CONFIRMED in code + reports:
- depth-1 closure: the triple point is Intersect3Planes on the REP planes
  (2599-2601); rep = min-index face per planeId (2457-2465).  No point is built
  from constructed points for any DECISION.
- coplanar-adjacent-faces degenerate case: planeId collapses exactly-coplanar
  faces to one id; a seam whose partner shares planeId is SKIPPED (2530-2534,
  "collinear/degenerate").  The shared edge line then comes from its input
  endpoints (degree 1) - the degenerate case LOWERS degree, as the task predicted.
- skew crossings = artifacts: exjunct/symwalk measured the phantom crossings'
  carrier lines share NO plane (3D-skew), 11-1700 eps from any exact triple;
  graze removes them by degree-0 carrier membership, "no nested-construction
  predicate, no new form."  The pipeline FAILS CLOSED on them (it never
  constructs the depth-2 point).
- sub-eps collapse is a DOUBLE eps threshold (uniform hash grid cell=eps, 251;
  la::length(pos-endpoint)<=eps decline, 2620-2624), NOT an exact predicate.
  This is the ONE genuinely-inexact decision.  Codex's degree-26 is HYPOTHETICAL
  (only if this were exactized).

Cited sites spot-checked REAL (no hallucinated line numbers): 1648 SegPlanePoint,
1729 faceN, 2292 Intersect3Planes, 2599 triple-emit, 1144 EdgePiercesTri,
2620 sub-eps decline, 2814/2834 JunctionSplitsOnSegment foot, 3009 pos2in
backstop, 2360-2374 shipped mixed orient2d.

## THE AGREED TABLE (raw input basis; plane-coeff basis in parens)

object            depth  homog coords         predicate                   degree
input vertex      0      (1,1,1;0)            orient3d(4 input)           3
face plane        1      n:2, d:3             pt-on-plane(input)          3
seam/edge line    1      dir:4, moment:5      -                          -
  (coplanar adj)  0      dir:1, moment:2      (degenerate; from endpoints)
seam endpoint     1      (4,4,4;3) or (7,6)   endpoint-on-plane           6
plane triple      1      (7,7,7;6)            triple-on-plane             9 (4)
  (Cramer)                                    triple in-triangle          8   <-SHIPPED max
                                              order-on-line (2 triples)   13
                                              orient2d(3 triples)         20 (9)  <-FORM cap
winding           -      (sum of signs)       pierce orient3d             3
sub-eps collapse  -      NOT exact            eps double threshold        - (26 if exactized)

MAX EXACT PREDICATE DEGREE:
  raw input basis        = 20  (form capacity; orient2d of 3 plane triples)
  plane-coeff basis      = 9   (the code's "degree-9 constructed-point orient2d")
  SHIPPED decisions only = 8   (constructed triple vs input triangle edges)

## MINIMAL EXACT INVENTORY (2 forms, per code = per doc)

1. HomogOrient3DSign<true> = orient3d(4 input pts), degree 3 raw.  Covers pierce/
   enumeration, coplanarity, input in-triangle, winding crossing.  + single-global
   SoS tie-break (same width-3 accumulator, K==0 group IS the exact orient3d).
2. HomogOrient2DExact = homogeneous in-face orient2d, W-sign corrected, degree 9
   plane-coeff (20 raw at full generality; 8 as shipped).  Covers in-face
   subdivision / seam straddle / cell walk.  Input verts = W==1 instantiation.
A 3rd form (constructed-point-on-plane, deg 6/9) was DELIBERATELY AVOIDED (graze:
degree-0 combinatorial carrier membership instead).  Design doc: adding any form
beyond these two is an explicit OWNER DECISION (Regularize3D.md, kAccumLimbs=320).
