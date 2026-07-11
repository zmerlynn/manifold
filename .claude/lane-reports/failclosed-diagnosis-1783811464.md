# Fail-closed diagnosis crucible - corpus single-mesh cases

Branch explore/sweep-plane-3d-v3 @ fc9f8df7. Solo, clean context. Sole editor.
Suite green at start: 601 = 600 pass + 1 documented Gate4c skip.

## Mission

Four real-mesh single-inputs (self-overlap removal) all terminate at
FatalReason::NonManifoldEmission, cause never attributed. DIAGNOSE each:

- test/models/Offset1.obj (2287 verts, ~1.6s) - Corpus_Offset1_Recorded
- openscad-nonmanifold-crash.obj (717 verts, ~10s) - Corpus_OpenscadNonmanifold_Recorded
  (DELIBERATELY non-manifold input - record input status)
- self_intersectA.obj (8582 verts, ~35s) - Corpus_SelfIntersectA_Recorded
- self_intersectB.obj (8582 verts, ~35s) - Corpus_SelfIntersectB_Recorded

Classify each: WALL A / WALL B-feeding / NEW CLASS <name>, with a re-runnable
evidence chain.

## Known wall signatures (from the journals)

WALL A (emission-closure junction clusters; docs SweepEmit3D "Emission-closure
crucible" + "3D-IDENTITY EXTENSION"): near-degenerate arr.verts clusters at cap
planes; cap critical sits sub-eps (0.3-0.5 eps) from the true 3D vertex V; a
steep track (slopes ~6-125) amplifies that offset to 60-200 eps; symptoms =
unbalanced fans at micro-edges, sliver cap triangles, weld-frozen-vs-track twin
divergence, retained loops splitting hairs. FATAL ARM: "unresolvable sheet
contact" (SplitTouchingSheets). Carriers: Havocglass8 (1F/2B), GenericTwin7863
(1F/0B).

WALL B (near-coplanar grouping miss): fires EngineIdConflict at BuildSlabs - a
DIFFERENT fatal code than these four report. If wall-B structure appears here it
must feed a downstream failure = a finding itself.

## The six NonManifoldEmission arms (overlap3.cpp), by detail string

1. "strip chain count disagrees with piece count" (EmitStrips :675)
2. "empty strip chain" (EmitStrips :681)
3. "cap boundary walk failed to close" (ComputeCap/OutEdgesToPolygons :841)
4. "cap triangulation returned invalid indices" (ComputeCap/TriangulateCap :844)
5. "unresolvable sheet contact" (SplitTouchingSheets :1114) = WALL A signature
6. "emitted triangulation not 2-manifold" (Is2Manifold gate :1127)

SplitTouchingSheets sub-arms (all -> arm 5): unbalanced (fwd!=bwd, :938),
sliver (len==0, :960), tie (gap<1e-9, :984), non-alternating (cur.fwd==nxt.fwd,
:985), unpaired-after (p<0, :999).

FIRST STEP per case: capture the exact detail string (which arm), THEN
instrument that arm.

## Method (established patterns)

Probe build: g++ -std=c++17 -O2 -DMANIFOLD_PAR=1 /tmp/foo.cpp -I../include -I../src
-o /tmp/foo ./src/libmanifold.so -Wl,-rpath,$PWD/src -ltbb  (from build/).
Single-mesh compose = ReadOBJ(fin) -> Manifold::Impl(mesh) -> eps =
EpsilonFromScale(impl.bBox_.Scale(),1000) -> RemoveOverlaps3D(impl, eps).
Instrumentation TEMP (env-gated fprintf, marked "// TEMP DEBUG", removed before
commit). SPLITDEBUG (failing edge + incident tris at %.17g + displaced-partner
search) and CAPDUMP (full cap input/output at target x) patterns per
coplanar-crucible Steps 12-13; junction-cluster characterization (cluster spread
in eps, dist crit->nearest arr.vert, track slopes) per identity-extension.

Cap: ~4 instrumented runs/case, else record UNATTRIBUTED honestly. Diagnosis
only - NO fixes, NO behavior changes. Record fix-shape hypotheses only.

## Log

### Step 1: capture the exact detail strings (DONE)

ABOUT-TO: /tmp/fcd_detail.cpp - run all four singles, print fatal reason +
detail + counters, plus openscad input status. FIRST STEP mandated per case.

OBSERVED (one probe run, all four):
```
Offset1.obj      inVerts=2287 inStatus=0 eps=4.50e-08
    FATAL reason=3 "unresolvable sheet contact" capArr=1599 idConflicts=0 subEpsDropped=18040
openscad         inVerts=717  inStatus=0 eps=3.52e-10
    FATAL reason=3 "unresolvable sheet contact" capArr=3493 idConflicts=0 subEpsDropped=9590
self_intersectA  inVerts=8582 inStatus=0 eps=2.75e-12
    FATAL reason=3 "unresolvable sheet contact" capArr=8933 idConflicts=0 subEpsDropped=62071
self_intersectB  inVerts=8582 inStatus=0 eps=2.75e-12
    FATAL reason=3 "unresolvable sheet contact" capArr=8927 idConflicts=0 subEpsDropped=61162
```

CONCLUSION:
- ALL FOUR die at the SAME arm: "unresolvable sheet contact" =
  SplitTouchingSheets (overlap3.cpp:1114) = the WALL A top-level signature
  (Havocglass8/GenericTwin7863 die here too).
- idConflicts=0 on all four -> NO EngineIdConflict anywhere, so WALL B is NOT
  feeding through (no near-coplanar grouping-miss signal reached the engine).
- openscad-nonmanifold-crash imports as inStatus=0 = NoError (Error enum:
  NoError==0). Despite the filename, Manifold's importer accepts it as a valid
  closed 2-manifold. DIAGNOSIS FRAME UNCHANGED: to RemoveOverlaps3D the input
  is a clean manifold; the failure is genuinely in OUR emission, same as the
  other three. (The "nonmanifold-crash" name is an OpenSCAD-crash provenance
  label, not an input-validity claim against us.)
- These are the SAME four the perf-campaign made terminate; that journal
  ASSUMED "steep-track junction-spread residual" (WALL A) without diagnosing
  the sub-arm. Top-level arm now confirmed == WALL A's arm; the sub-arm and
  the geometry still need to be proven WALL A vs a new class surfacing here.

NEXT: instrument SplitTouchingSheets sub-arms (unbalanced / sliver / tie /
non-alternating / unpaired) + dump the first-failing edge geometry, to decide
WALL A (unbalanced micro-edge at a cap plane, near-degenerate cluster) vs a
new class (e.g. macro material overlap = non-alternating on genuine
self-intersection geometry).

### Step 2: instrument SplitTouchingSheets sub-arms + first-failure dump (DONE)

ABOUT-TO: env-gated histogram in SplitTouchingSheets (overlap3.cpp) counting
every fan/boundary edge into UNBAL / SLIVER / TIE / NONALT / ringOk, plus a
first-failure geometry dump (edge verts, len/eps, dx, incident tris at
%.17g). OVERLAP3_SPLITDBG on, OVERLAP3_EPS per-case.

OBSERVED (one run, all four; counts are over the WHOLE emitted triangulation):
```
Offset1   nTri=246047  edges=369629  balMultiFan=184 ringOk=27 | UNBAL=1495 (boundary1he=1491) SLIVER=0 TIE=126 NONALT=31
openscad  nTri=1055600 edges=1583986 balMultiFan=155 ringOk=32 | UNBAL=1508 (boundary1he=1498) SLIVER=0 TIE=102 NONALT=21
self_A    nTri=4156680 edges=6237602 balMultiFan=40  ringOk=19 | UNBAL=5246 (boundary1he=5245) SLIVER=0 TIE=11  NONALT=10
self_B    nTri=3881890 edges=5824765 balMultiFan=101 ringOk=64 | UNBAL=4066 (boundary1he=4064) SLIVER=0 TIE=18  NONALT=19
```
First-fail edges (what the real code returns on, in map order):
- Offset1  TIE   edge(v501,v502)   len/eps=291.7 dx/eps=55.1 (4 tris, one third-vert v503 appears twice FWD+BWD = a fold)
- openscad NONALT edge(v3,v37)     len/eps=42.4  dx=0 (cap edge; tri33 third on-plane, tris 4752/4754/4756 third 2340 eps away in x = strips)
- self_A   TIE   edge(v1629,v1630) len/eps=15.0  dx=0 (cap edge; 4 strip tris from both neighbor slabs, 2 tangent)
- self_B   NONALT edge(v6621,v6622) len/eps=8.4  dx=0 (cap edge; 4 strip tris, non-alternating)

CONCLUSION:
- The DOMINANT structural defect in ALL FOUR is UNBAL, and ~99.7-99.98% of
  UNBAL are BOUNDARY edges (hes=1, single halfedge) = OPEN HOLES. Thousands of
  them per case (1491 / 1498 / 5245 / 4064). The emitted surface is riddled
  with holes; ANY one is fatal. This is NOT the isolated single-junction
  picture of WALL A's pair fixtures (Havocglass8 1F/2B = one fan; GenericTwin
  1F/0B = one strip-less cap edge) - it is the same emission-closure disease
  but PERVASIVE (100s-1000s of sites) on dense real geometry.
- The pipeline TRIPS first on a TIE (Offset1, self_A) or NONALT (openscad,
  self_B), not on the boundary edges, purely because map order (sorted vertid
  pair) reaches a ring-failing edge first. But even if that first edge were
  repaired the surface still has thousands of holes -> the boundary edges are
  the real disease.
- The first-fail edges are MICRO cap edges (dx=0, len 8-292 eps) where 4
  strips from adjacent slabs meet tangentially (TIE) or with non-alternating
  material (NONALT) and NO cap triangle bridges them - the exact cap/strip
  junction-closure signature of WALL A, at micro-eps scale. Offset1's is a
  fold (v503 as both a FWD and a BWD third-vert).
- SLIVER=0 everywhere (no third-vert-on-edge-line degeneracy).

OPEN QUESTION for classification: are the thousands of boundary holes MICRO
(sub-~200-eps cap-closure gaps = WALL A at scale) or are there MACRO holes
(whole missing faces/caps = a distinct emission-COMPLETENESS class)? That
scale measurement decides WALL A-at-scale vs NEW CLASS. -> Step 3.

### Step 3: boundary-edge (hole) length-scale distribution (DONE)

ABOUT-TO: extend the histogram with boundary-hole len/eps buckets + cap(dx=0)/
strip split, and dump a few macro holes (len/eps>1e3) with incident-triangle
classification (CAPtri = all verts same x; STRIPtri = third vert at neighbor x).

OBSERVED (boundary-hole len/eps buckets [<1,<10,<1e2,<1e3,<1e4,<1e6,>=1e6]):
```
Offset1  [0,0,11,50,60,196,1174] cap=1491 strip=0   (1174 MACRO >=1e6, 11 micro)
openscad [0,2,12,2,8,132,1342]  cap=1498 strip=0   (1342 MACRO, ~16 micro)
self_A   [0,0,3,0,0,0,5242]     cap=5245 strip=0   (5242 MACRO, 3 micro)
self_B   [0,0,6,0,0,0,4058]     cap=4064 strip=0   (4058 MACRO, 6 micro)
```
EVERY boundary hole is a cap-plane edge (dx=0); ZERO strip holes. Macro-hole
incident-triangle classification:
- Offset1: ALL STRIPtri - strips present at the cap plane, cap FILL missing.
  8 dumped holes share ONE cap plane x=4747.0582796500457, a chain
  v4316..v4327 = one cap-loop boundary left open.
- self_A: ALL STRIPtri, one cap plane x=6.05e-4, chain v796695..v796703.
- self_B: ALL STRIPtri, one cap plane x=-0.60051, chain v1223611..v1223619.
- openscad: MIXED - some CAPtri (cap present, strips missing: v779-v780 etc.)
  AND some STRIPtri (strips present, cap missing: v781-v71744 etc.).

CONCLUSION - this REFUTES the perf-campaign's assumption (it labelled all four
"steep-track junction-spread residual" = WALL A). The dominant defect is MACRO
cap-plane boundary NON-CLOSURE (len/eps 1e6-1e9, thousands per case), all in
cap planes, forming open cap-loop chains. WALL A is MICRO (60-200 eps isolated
junction spread). Only a tiny tail here is micro (3-16 holes + the TIE/NONALT
first-fails) = WALL A-family minority. Scale + pervasiveness + all-cap-plane +
strips-present/cap-missing dominant = a DISTINCT class from WALL A.

The dominant pattern (strips present, cap fill missing) is the INVERSE of the
M4-close recorded dead-zone (a) ("emits a macro cap whose corners no strip
carries" = cap present, strips missing = the openscad CAPtri minority). So the
dominant macro holes are NOT dead-zone (a); they are a new, un-recorded shape.

OPEN FORK (design-relevant): are these macro-hole cap planes CANONICAL
criticals where ComputeCap RAN but produced fill that misses the strips (a
closure-argument breakdown), or SKIPPED non-canonical criticals (ci != li+1,
interior to a sub-eps run = dead-zone at scale)? -> Step 4 (CAPAUDIT/CAPDUMP).

### Step 4: cap-emission audit at macro-hole planes (DONE)

ABOUT-TO: instrument EmitCaps (canonical/skipped per critical + global counts)
and ComputeCap (CAPDUMP of pieces/edges/tris at a target x + global zeroFill
audit). Target x per case = a macro-hole plane. Plus CONTROLS: run the
resolving Offset2/3/4.

OBSERVED:
- The macro-hole cap planes are CANONICAL criticals (ci==li+1) that BRIDGE two
  built slabs separated by a sub-eps RUN of unbuilt slabs. Offset1 target:
  ci=596 canonical, li=595(nP=69), ri=604(nP=78); criticals 597..604 all
  non-canonical. The cap RAN over 69+78=147 pieces (294 raw verts) but emitted
  capPlusEdges=0 capMinusEdges=0 capTris=0. self_A/self_B identical shape
  (0 fill over 694/796 pieces). openscad target emitted a TINY partial fill
  (6 edges, 4 tris over 271 pieces).
- Non-canonical fraction is HUGE: Offset1 8287/9886=84%, openscad 5302/8795=60%,
  self_A 54157/63090=86%, self_B 48468/57395=84% of criticals are interior to
  sub-eps runs.
- GLOBAL zeroFill audit: Offset1 1599/1599 (100%) canonical caps emit empty
  fill; self_A 8933/8933 (100%); self_B 8927/8927 (100%); openscad 1861/3493
  (53%). Nearly all with >=10 pieces/side.
- CONTROL (decisive): the RESOLVING Offset2/3/4 ALSO have 100% zeroFill caps
  (Offset2 584/584, Offset3 16/16, Offset4 22/22) and resolve cleanly. So an
  EMPTY cap is the architecture's NORMAL (retained region genuinely unchanged
  across most slab boundaries; closure is by strip-to-strip weld, not cap).

CONCLUSION (corrects an in-flight over-claim): the defect is NOT "cap collapse"
- empty caps are normal and resolving cases have them too. The defect is
  STRIP-TO-STRIP WELD FAILURE across the sub-eps runs. The run bridged at
  Offset1's target spans crits[596]..crits[604] = 4747.0582796500457 ..
  4747.0582797133748 = 6.33e-8 in x = ~1.4 eps (eps=4.5e-8). A multi-sliver
  run whose TOTAL width EXCEEDS eps: the LEFT slab's strips end at x=crits[596]
  and the RIGHT slab's strips start at x=crits[604], >eps apart in x, so
  BuildImpl's eps-weld cannot fuse them -> the shared cap-plane loop splits
  into open boundary edges. This is the RECORDED M4-close dead-zone (c)
  ("Multi-sliver runs (total gap > eps) weld strip-to-strip across more than
  eps"), documented as "unreachable in the suite" - here shown PERVASIVELY
  REACHABLE on dense near-x-perpendicular real geometry (offset walls /
  self-intersections pack many criticals into eps-thin x-bands spanning macro
  yz). NOT WALL A (micro steep-track). Must VERIFY the weld-gap directly (not
  infer - the cap-collapse over-claim was just refuted by control). -> Step 5.

### Step 5: verify the weld-gap directly (displaced-partner search) (DONE)

ABOUT-TO: for each dumped macro hole endpoint pa, search all output verts for
the nearest one in yz that sits >eps away in x (the un-welded strip endpoint on
the far side of the run). Report partnerDx/eps and partnerYz/eps.

OBSERVED - EVERY macro hole has a displaced partner:
```
Offset1  partnerDx/eps=1.406  partnerYz/eps=0  (matches run width crits[604]-crits[596]=1.4 eps)
openscad partnerDx/eps=1.25..2.48 partnerYz/eps=0
self_A   partnerDx/eps=1.685  partnerYz/eps=0
self_B   partnerDx/eps=1.177  partnerYz/eps=0
```
Controls Offset2/3/4: ZERO macro holes (they resolve).

CONCLUSION - DEAD-ZONE (c) CONFIRMED, directly, not inferred: every macro-hole
endpoint has a partner vertex at the SAME (y,z) (partnerYz=0) but 1.2-2.5 eps
away in x (> eps). That partner is the far-side strip endpoint of a sub-eps RUN
whose TOTAL x-width exceeds eps; BuildImpl's eps-weld (3D distance <= eps)
cannot fuse two points >eps apart in x, so the shared cap-plane loop splits
into open boundary edges. This is the RECORDED M4-close dead-zone (c) ("Multi-
sliver runs (total gap > eps) weld strip-to-strip across more than eps"),
documented as "inherited by construction and unreachable in the suite." It is
PERVASIVELY REACHABLE on these dense real meshes and is the dominant cause of
all four fatals. The (correctly empty) canonical cap cannot bridge it: caps
fill L!=R symmetric differences, but here region(L)==region(R) (same loop) at
two x-planes >eps apart, so no cap is emitted and closure would have to come
from a weld that eps forbids.

## FINAL CLASSIFICATION

All four cases: NOT WALL A, NOT WALL B. They are the RECORDED M4-close
DEAD-ZONE (c) - multi-sliver-run strip-weld gap - made pervasively reachable on
dense near-x-perpendicular real geometry (offset walls, self-intersection
sheets pack many vertex-x criticals into eps-thin x-bands spanning macro yz;
84-86% of criticals collapse into sub-eps runs, many summing to >eps).

- Offset1                    : DEAD-ZONE-C (clean: 100% empty caps, all-STRIPtri holes, partnerDx~1.4 eps)
- openscad-nonmanifold-crash : DEAD-ZONE-C (richer: 53% empty caps, mixed CAPtri+STRIPtri, partnerDx~1.3-2.5 eps; input imports valid, inStatus=0)
- self_intersectA            : DEAD-ZONE-C (clean: 100% empty caps, all-STRIPtri, partnerDx~1.7 eps)
- self_intersectB            : DEAD-ZONE-C (clean: 100% empty caps, all-STRIPtri, partnerDx~1.2 eps)

EVIDENCE CHAIN (re-runnable): the SPLITDBG histogram + macro-hole/partner dump
in SplitTouchingSheets and the CAPAUDIT/CAPDUMP in EmitCaps/ComputeCap (all
TEMP DEBUG, env OVERLAP3_SPLITDBG / OVERLAP3_CAPAUDIT / OVERLAP3_CAPDUMP_X /
OVERLAP3_EPS), driven by /tmp/fcd_detail.cpp (single-mesh compose, four cases +
Offset2/3/4 controls). Removed before commit; the notebook records the numbers.

WHY THE PERF-CAMPAIGN LABEL WAS WRONG: it assumed "steep-track junction-spread
residual" (WALL A) without instrumenting the arm. WALL A is MICRO (60-200 eps
isolated junction, unbalanced fan at a near-degenerate arr.verts cluster). This
is MACRO (len/eps 1e6-1e9) and PERVASIVE (thousands of holes), driven by run
TOTAL-width > eps, not by track slope. A WALL A-family MICRO tail co-exists
(the TIE/NONALT first-fail edges + 3-16 micro boundary holes per case), but it
is a small minority; the defining/dominant defect is dead-zone (c).

FIX-SHAPE HYPOTHESES (recorded, NOT implemented - diagnosis only):
1. Bound the run: when a sub-eps run's TOTAL width exceeds eps, do not merge it
   into one gap - either build an interior section (defeats the sub-eps skip)
   or emit an explicit SubEpsFeature/named guard when total-width > eps and the
   flanking sections differ. (The SubEpsFeature guard today only catches a
   single face living ENTIRELY inside a run; it misses content distributed
   across straddling faces - which is why these fail as NonManifoldEmission,
   not SubEpsFeature.)
2. Weld across the run explicitly: bind the far-side strip endpoints to the
   near-side by run-provenance (same (y,z), the partner search finds them at
   partnerYz=0) rather than by eps-distance, so the >eps x-gap does not block
   the weld. This is a targeted assembly-weld change, orthogonal to the WALL A
   3D-identity snap.
3. Note the co-located WALL A-family micro tail will still need the 3D-identity
   work; dead-zone (c) is the dominant blocker to address first.

CONSTRUCTED MINIMAL FIXTURE: deferred (diagnosis-only; the four corpus singles
are recorded reproducers). Construction shape for the fix-arc: a closed mesh
with a near-x-perpendicular wall whose vertices give 3+ criticals each spaced
< eps in x but summing to > eps, with macro yz content changing across the run
so the flanking slabs' strips are the same loop but > eps apart in x. Not
attempted here to respect the do-not-spiral cap.

### Step 6: cleanup + record contracts + doc + commit (DONE)

- Restored src/overlap3.cpp to clean HEAD (all TEMP DEBUG was instrumentation
  only; verified 0 remnants of TEMP DEBUG/SPLITDBG/CAPAUDIT/CAPDUMP/g_capsTotal
  in src/, and the rebuilt .so has 0 debug symbols). git diff HEAD --
  src/overlap3.cpp is EMPTY -> zero behavior change.
- Updated the four Corpus_*_Recorded comments with the DEAD-ZONE (c)
  classification + doc pointer (assertions untouched); corrected the shared
  block comment that had misattributed to the steep-track WALL A residual.
- Appended "Corpus fail-closed attribution" to docs/SweepEmit3D.md (by kind,
  magnitudes only per the no-counts house rule).
- clang-format clean on the changed test file; FULL suite green: 601 = 600
  pass + 1 skip (Gate4c) = baseline.
- Single commit (comments + doc + notebook only), NOT pushed.

DEAD ENDS / corrections recorded (crucible discipline):
- Mid-Step-4 I nearly concluded "wholesale cap-fill collapse is the defect"
  (100% of canonical caps emit empty fill on the failing cases). The RESOLVING
  control (Offset2/3/4) REFUTED it: they too have ~100% empty caps and resolve
  cleanly. Empty caps are the architecture's normal; the real defect is the
  strip-weld gap. Verified directly in Step 5 (displaced-partner search) rather
  than left as inference.
- Considered but did NOT build a synthetic minimal fixture (respecting the
  do-not-spiral cap; the four corpus singles are recorded reproducers and this
  maps to the already-documented dead-zone (c)). Construction shape recorded
  in Step 5 for the fix-arc.
