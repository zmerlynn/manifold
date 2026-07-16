# VERIFY lane: symbolic-extent plumbing (E1), Regularize3D

Adversarial push-gate for the E1 plumbing landing.  Branch
explore/sweep-plane-3d-v5.  Range f53c9e3f (base, pushed) .. 3e11ed8a (HEAD, the
code commit under review + research notebooks).  Grounding:
e1-plumb-1784174600, nomerge-1784160552, homog-design-1784161597,
docs/Regularize3D.md.  Canonical READ-ONLY except this notebook (my only commit).
DID NOT push.

## SETUP
- From-scratch Release builds, MANIFOLD_PAR=ON, -j8, both SHAs:
  /tmp/plumb-verify (HEAD, from the clean canonical tree) and
  /tmp/plumb-verify-base (base f53c9e3f, git-archived source).  Both compile
  clean: 0 warnings, overlap3.cpp included; identical binary size.
- Code diff in range (non-notebook): src/overlap3.cpp (+187/-100 net) and
  docs/Regularize3D.md only.  The change is confined to EnumerateTriplePoints +
  its call site + the kernel tripwire comment.

## 1. SUITE + HEAVY CARRIERS (both trees, PAR)
- Overlap3 suite (PAR, 6GB ulimit): HEAD 32/32 PASSED, BASE 32/32 PASSED.
- GT7081 default RESOLVES at HEAD under the tight `ulimit -v 4000000; timeout
  900` (75.7s) and under load in the suite (89.9s).  BASE resolves too.
- openscad fail-closed both SHAs (NonManifoldEmission "unresolvable sheet
  contact"); large-comp openEdges=21 on BOTH SHAs (F4B_DUMP census).  Triple
  count 157 (BASE) -> 156 (HEAD): the one-phantom refinement, opens unchanged.

## 2. CROSS-BUILD FNV (independent harness, the deferred compare)
Harness links the built libmanifold.so of each tree, resolves each carrier,
FNV-1a over vertPos raw bytes + halfedge Start ids (the BitIdenticalMesh
identity).  base vs HEAD, BYTE-IDENTICAL on every carrier:
  siA / siB / GT7863 / GT7081 all RESOLVED, identical FNV; openscad identical
  fail-closed.  Mechanism directly confirmed: GT7081 default F4B_TRIPLES
  distinct=0 incidences=0 on both dirty components (empty seamTriples ->
  EnumerateTriplePoints emits nothing -> provable downstream no-op).

## 3. EXTENT ORACLE DIFFERENTIAL (ported)
Standalone harness: the PRODUCTION sos:: kernel copied verbatim
(src/overlap3.cpp 518-1036), HPointStrictlyInTri/DominantAxis reproduced;
generates weighted triples (25% random, 37.5% near-parallel wedge, 37.5%
near-boundary X-within-1e-6-of-edge) and dumps each case + the kernel verdict as
bit-exact %a hex.  An INDEPENDENT Python exact-rational oracle (fractions)
recomputes the strict-interior verdict.
  120,000 cases; 8,282 filter escalations (the exact HomogOrient2DExact path
  IS exercised); interior counts agree (28,077 both sides); DISAGREEMENTS = 0.

## 4. CODE REVIEW of 3e11ed8a
- (i) Only the two landed instantiations.  Every exact-arithmetic site routes to
  HomogOrient3DSign (inst.1, via Orient3DExactSign / Orient3DSoS /
  ExactOrient2DDrop) or HomogOrient2DFilter/Exact (inst.2).  ExactSeamsCross /
  HPointStrictlyInTri are filter-first callers of inst.2 only.
  ECramerHPoint/CramerHPoint are pre-existing POINT constructions (not
  predicates).  No third exact FORM.
- (ii) Default is refinement-only: `straddle && ExactSeamsCross(clip=3)` is a
  subset of the rounded gate (can only drop/confirm), and a subset of the
  pre-plumb X-in-f default (X in f,g,h implies X in f).  Never adds.
- (iii) E1_PURE cleanly isolated: `straddle = kE1Pure ? true : ExactSegProperCross(..)`
  is the ONLY env dependence of the crossing decision; the default path (neither
  env set) is unaffected by E1_PURE.  clip = kE1Off ? 1 : 3; guard gated on !kE1Off.
- (iv) Aliasing guard: `la::length(pos - seam-endpoint) <= eps` against the four
  LEVEL-0 seam endpoints (k1/k2 p0/p1), all sharing X's seam carrier - the
  nomerge SITE-1 within-construction eps collapse it explicitly blesses as a
  provenance-consistent proxy (NOT the SITE-2 cross-provenance metric merge it
  prohibits).  eps = the component weld radius.  Empirically load-bearing only
  under E1_PURE (catches SelfIntersectB's 2e-15 alias; does NOT over-catch
  GT7081's 1.4e-7 above-eps crossings - confirmed by the E1_PURE mutation).
  NOTE: it is literally a distance-to-level-0-vertex <= weld-eps test; the
  "level-0 incidence" framing is the endpoint it keys on, not a combinatorial
  predicate.  Consistent with the design; not a defect.
- (v) Fail-closed arms intact (WallCensus return-false on open edges;
  NonManifoldEmission "unresolvable sheet contact") - untouched by the diff.
  Tripwire comment + docs updated to name ExactSeamsCross/HPointStrictlyInTri as
  inst.2's production caller; inventories consistent.
- (vi) All instrumentation env-gated (E1_MEASURE / E1_OFF / E1_PURE / F4B_*);
  default prints nothing.  No FNV/memcmp PRINTED (line-260 FNV is a pre-existing
  coordinate hash, not in this diff).  0 non-ASCII bytes in src, docs-changed
  lines, and the commit message.

## 5. MUTATIONS (levers bite)
- E1_PURE, GT7081: regresses to fail-closed, F4B_CENSUS openEdges=4 (the
  documented decisive negative); default GT7081 resolves.
- E1_PURE, openscad: distinct triples 156 -> 170, test FAILS (worse).
- E1_OFF, openscad: distinct 157, openEdges=21 (pre-plumb straddle-and-X-in-f).
- E1_MEASURE, openscad large comp: default=468 prePlumb=469 pureExact=536
  phantomsDropped=1 extentClipRemovals=5234 - exactly the lane's STAGE-2 census.

## 6. DETERMINISM
- siA: 1t == 4t (tbb::global_control), identical FNV, matches default.
- GT7081: 1t == 4t == default, FNV=7f15398ea159523e.
- openscad census stable across two runs (openEdges=21 both).

## VERDICT: PASS (push-safe)
All five claims reproduce independently; suite green both trees PAR; cross-build
byte-identical on the resolving corpus; oracle 0-disagree at 120k; mutations
anchor the levers; PAR-deterministic.  One wording nuance on the aliasing guard
(a weld-eps distance to level-0 endpoints, the blessed SITE-1 proxy) - design
correct, not a blocker.  DID NOT push.
