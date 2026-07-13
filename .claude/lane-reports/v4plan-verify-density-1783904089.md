# V4 IMPL crucible round-1 - ADVERSARIAL DENSITY LANE (R-TRIPLE + R-DENSITY + R-ONCE-VS-STRICTFP)

Artifact under attack: docs/V4ImplPlan.md (stage A, the RISKS) grounded on the
stage-0 skeleton (.claude/lane-reports/v4impl-stage0-1783901754.md). The skeleton
measured the identity table at the VERTEX-JUNCTION scope only (~18.5k junctions on
GT7081, 0.02s). The plan's R-TRIPLE names the unmeasured heavier object: pass 1 must
resolve the CROSSING set into junction identities, not just canonicalize verts. This
lane measures the crossing set before the build bets on it.

Workspace: /tmp/v4plan-density (rsync copy of canonical @ 2d404d57, fresh cmake -B
vbuild Release -j4). Canonical repo READ-ONLY except this notebook. Copy reverted at
wrap. Memory discipline: every heavy run under (ulimit -v 6000000; timeout 900),
output to scratchpad files, bounded reads, one heavy case at a time.

## Probe mechanism (copy-only, env-gated, reverted)

The plan's crossing/critical set is x-ONLY as stored (arr.criticalXs holds only the
x-value of each seam-seam crossing / degenerate contact - overlap3.h:184). To resolve
crossings to 3D JUNCTIONS I instrumented FindSeams to capture, for EVERY criticalX it
emits, the full 3D point + an order-free input-side combinatorial key + a category:
  - cat 0 = DEGENERATE CONTACT (seamLen <= eps face-pair near-tangency; key = face
    pair) - overlap3.cpp:426-445
  - cat 1 = M1 SEAM-SEAM TRIPLE (two seams sharing a face cross; key = the sorted 3
    distinct faces = the input-side triple-point identity) - overlap3.cpp:451-470
  - cat 2 = COPLANAR-GROUP SKELETON crossing (key = sorted 4 vert-ids of the two
    crossing edges) - overlap3.cpp:472-492
Added RemoveOverlaps3D_CaptureCrossings (canonicalize + FindSeams ONLY, no BuildSlabs,
so it dodges the GT7081 ArrangementBudget blowup). Probe test V4DensityProbe
(env OV3_V4DENSITY=<fixture>, OV3_V4DENSITY_PERM=1) clusters the captured 3D points by
an eps-ball grid union-find (cell=eps; dense >4000-pt cells star-unioned to bound cost
- one such blob on GT7081), lex-min canonical per cluster (coordinate-determined),
measures cost + peak RSS (VmHWM), distinct input-keys, the Q2 slab reduction, and C1
under input permutation (vertex remap + triangle reorder). CAPTURE VALIDATION: captured
count == totCritXs EXACTLY on every fixture (1,374,649 on GT7081) - the instrumentation
sees every critical push and nothing else.

## RESULTS

| fixture | nVerts | nSeams | totCritXs | deg / M1 / skel      | CONSTRUCTED junctions | cluster cost      |
|---------|--------|--------|-----------|----------------------|-----------------------|-------------------|
| GT7081  | 18643  | 3043   | 1,374,649 | 1,374,634 / 0 / 15   | 16,889                | 9.2s / 134MB      |
| selfA   | 8920   | 338    | 124,142   | 124,142 / 0 / 0      | 8,546                 | 0.04s / 25MB      |
| selfB   | 8920   | 338    | 122,324   | 122,324 / 0 / 0      | 8,555                 | 0.03s / 25MB      |
| Havoc   | 212    | 89     | 1,431     | 1,422 / 1 / 8        | 101                   | <0.01s / 10MB     |
| GT7863  | 596    | 391    | 4,476     | 4,384 / 3 / 89       | 330                   | <0.01s / 10MB     |

Input-side distinct keys: GT7081 687,332 (M1 triples 0); selfA 62,071; Havoc 720.
C1 permutation (seeds 1,7): GT7081 nJ=16889 both, countDiff=0, worstDevEps=0.000689,
exactMiss=1/16889. selfA/selfB/Havoc/GT7863: countDiff=0, worstDevEps=0, exactMiss=0
(BITWISE).
Q2 slabs (before-raw / before-epsMerge / after-bundle-raw / after-bundle-epsMerge):
GT7081 231,616 / 19,530 / 34,850 / 19,112 ; selfA 63,089 / 8,942 / 17,365 / 8,938 ;
selfB 57,394 / 8,936 / 17,320 / 8,928 ; Havoc 532 / 154 / 223 / 154.

## Q1 - R-TRIPLE : VERDICT SURVIVE (and the feared object is ABSENT)

THE decisive finding: on all three heavy carriers the ~1.37M / ~124k "criticals" are
99.99%+ DEGENERATE CONTACTS (near-tangent face pairs, seamLen <= eps), with ZERO M1
seam-seam TRIPLE points (GT7081 M1=0, selfA/B M1=0) and a negligible skeleton tail
(15 / 0 / 0). The triple-point completion R-TRIPLE feared would explode the table does
NOT exist to be resolved on these carriers.

1. The crossing/contact set resolves to BUNDLE/VERTEX scale, NOT crossing scale.
   GT7081: 1,374,649 -> 16,889 junctions (~81x), tracking nVerts=18,643 not
   critXs=1.37M. selfA: 124,142 -> 8,546 (nVerts 8,920). So even FOLDING THE FULL
   CROSSING/CONTACT SET into junction identities keeps the table bounded at ~vertex
   scale - the skeleton's R-DENSITY bound holds WITH crossings included, not just for
   verts. Cost is modest: 9.2s / 134MB on GT7081, sub-second elsewhere.

2. C1 order-invariance HOLDS at density for CONSTRUCTED crossings. The partition is
   BITWISE order-invariant (GT7081 countDiff=0, same 16,889 clusters under two
   permutations; selfA/B/Havoc/GT7863 exactly identical). The lex-min POSITION is
   order-invariant to within 0.000689 eps (ONE junction of 16,889 wobbles sub-milli-eps
   on GT7081; bitwise on the others). This is the R-ORDER-UNMEASURED-AT-DENSITY gap
   closed FOR CROSSINGS: consuming constructed contact coordinates does not break C1 -
   the near-tangent contacts are either clearly within eps or clearly not, so the
   0.0007-eps FP wobble flips no partition membership. The tiny position residual is
   exactly what C1's condition-1 "exact-predicate re-derivation from input" zeroes
   (derive each bundle's canonical from the input face planes, not the constructed
   point).

3. Input-side COMBINATORIAL derivation does NOT achieve bundle-scale. The order-free
   face-pair/triple key gives 687,332 distinct identities on GT7081 where geometry has
   16,889 junctions (~40x over-count; ~7x on selfA/Havoc). So the design's "derive
   identity input-side" CANNOT mean a pure input-combinatorial key - that over-counts
   the bundle. It must mean: GROUP by the constructed points (empirically order-safe,
   measured above) and derive each bundle's POSITION input-side. The plan's C1 wording
   ("Junction::pos is a PURE FUNCTION of the member INPUT coordinates") is satisfiable
   for the position, but the GROUPING consumes constructed coordinates and its
   order-safety is now MEASURED, not assumed - the plan should say so rather than imply
   the whole identity is input-combinatorial.

SURVIVE, with a factual correction the plan should absorb: the density mechanism is
DEGENERATE CONTACTS, not seam-seam crossings. The plan (1.3 "crossing bundles"),
R-TRIPLE, R-DENSITY, and the wallb-verify notebook all attribute the ~1.37M to "seams
crossing at a combinatorial pile" - mechanically WRONG (M1=0 on GT7081). The plan's
mechanisms still transfer (degenerate-contact x's are consumed x-only, same as crossing
x's, so crossing-bundle thinning's safety argument applies verbatim), but R-TRIPLE's
specific fear (triple-point completion adds ~1.37M junctions) is MOOT with evidence:
there are no triple points to complete on the heavy carriers, and the plan's decision
to scope stage A to vertex junctions + defer triple-point completion is VINDICATED.

## Q2 - R-DENSITY downstream : VERDICT SURVIVE (localizes to thinning; ~order-of-magnitude shrink CONFIRMED)

The plan's claim that the density localizes to v3's existing thinning/budget territory
holds, and the stronger conjecture (the CRITICAL SET itself shrinks under
bundle-clustering) is CONFIRMED:
- GT7081 slabs 231,616 -> 19,530 under eps-merge (~12x); bundle-clustering to junction
  canonical x gives the same floor (19,112). selfA 63,089 -> 8,942 (~7x); selfB 57,394
  -> 8,936. So YES: the bundle-interior criticals collapse to the junction's canonical
  x (the validated thinning result, now with the identity-table justification), and the
  critical set shrinks an ORDER OF MAGNITUDE. Stage A's clustering SUBSUMES the thinning
  decision - a plan amendment worth this evidence.
- PROJECTED (not measured - I did not re-run BuildSlabs with thinned crits): wallb-verify
  measured ~30M retained pieces at 231,616 slabs; a ~12x slab reduction concentrated in
  exactly the dense near-tangent bands where contacts pile should cut retained pieces
  proportionally, plausibly under the 4M ArrangementBudget floor - i.e. thinning could
  flip GT7081 from budget-refusal to sectionable. This is a projection; the load-bearing
  BuildSlabs-on-thinned-crits run is the follow-up.
- ArrangementBudget still warranted as the resource backstop (the plan is right): the
  ENDPOINT vertex clusters stay dense (skeleton maxCluster=77) after crossing thinning;
  the budget bounds THAT, not the (now-thinned) contact criticals.

## Q3 - R-ONCE-VS-STRICTFP : VERDICT NEED-CHANGE (specify + pin the per-slab junction MATCH)

Read from v3 code: BuildSlabs (overlap3_sweep.cpp:95) + ComputeSectionSegment (:41,
Interpolate at slab.xMid) + SlabResolver (overlap3.cpp:519, Extend via InterpolateSafe,
:569) + ComputeCap (resolver.Extend(piece.from/to, xCap), :785).

- Each built slab computes its section verts by Interpolate(edge, xMid) - a per-slab
  CONSTRUCTED coordinate. SlabResolver::Extend then re-derives each section vert's
  cap-plane image by InterpolateSafe along a per-eps-cluster Track to the critical x.
  THIS per-slab self-location is the wall-A locus: two adjacent slabs each Extend their
  OWN track to the shared critical -> two images of one junction -> the twin.
- What pass 2 must STOP recomputing, precisely: the Extend self-location of any
  JUNCTION-imaging section vert. Both adjacent slabs must take the junction's canonical
  from the identity table (yz at stage B; full 3D point incl. x at stage C), NOT each
  InterpolateSafe its own track.
- Does the plan say so? For the POSITION, YES (stage B pin 4a: two cap images
  bit-identical yz; stage C: one 3D point; section 1.4: self-location DELETES for
  junction-imaging verts). And it is HONEST that per-slab section stays for NON-junction
  geometry (1.4: "non-junction endpoints may keep track subdivision"; stage B: "fallback
  to the resolver for non-junction verts"). That is defensible: a non-junction section
  vert lies on ONE edge, and both adjacent slabs InterpolateSafe the SAME edge to the
  SAME critical x -> identical point, no twin. C2 constrains the CLUSTER canonicalization
  (done once in pass 1); a non-junction vert is not a cluster, so leaving its per-slab
  section is not a C2 violation.
- THE GAP R-ONCE-VS-STRICTFP correctly names: pass 2 introduces a NEW per-slab,
  per-consumer decision - "is this section vert junction J's image?" (the MATCH). The C2
  pin (canonical computed once) covers J's POSITION, NOT the per-slab match. v3's
  resolver matches section verts to tracks POSITIONALLY (la::length(pt - ft.p0) <= eps,
  :585-593). If pass 2's junction-match is likewise positional, it is a per-slab decision
  that CAN diverge across the two adjacent slabs (one matches J, the other falls through
  to track-extension) - and the twin RELOCATES into the match, exactly the variant-(ii)
  "which section vertex is the image of V" kill (MaintainedEmission3D probe 1). If the
  match is COMBINATORIAL (the track's defining edge is incident to arr.vert J), it is
  once-only-safe. The plan does NOT specify which and does NOT pin cross-slab match
  consistency.
- NEED-CHANGE: stage A/B must specify the junction-imaging MATCH as combinatorial /
  once-only (or prove the positional match is cross-slab-consistent) and PIN it. The
  current pins guarantee the canonical is computed once but leave the per-slab match -
  the residual self-location - unspecified. This is mild (a specification + one pin), not
  a BREAK: the identity table's C2 is sound; the gap is the CONSUMPTION-side match pass 2
  adds, which the plan currently under-specifies.

## WRAP / fence

- Copy reverted: git checkout src/overlap3.cpp src/overlap3.h test/overlap3_test.cpp in
  /tmp/v4plan-density (probe instrumentation + hook gone; diff empty).
- Canonical repo: NOTHING landed except this notebook.
- Logs: scratchpad gt7081_core.log, gt7081_perm.log, selfA_perm.log, selfB_perm.log.
- Dense-cap caveat: GT7081's 16,889 uses a star-union on one 324,721-point eps-cell
  (a massively degenerate near-coincident face stack = one bundle); the true eps-ball
  count is within a small factor and still vertex-scale - the order-of-magnitude claim
  (17k vs 1.37M) is robust to the approximation. selfA/B/Havoc/GT7863 used no cap
  (exact) and corroborate.
