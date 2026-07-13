# reg3d NEAR-COPLANAR research lane: closing B's last refusal class (design memo)

Branch explore/sweep-plane-3d-v5. READ-ONLY on the canonical repo (stage-4 lane
editing concurrently); all builds/probes in /tmp/reg3d-ncr (rsync snapshot),
fresh Release cmake, -j4. This file is a DRAFT design-memo section for the
orchestrator to fold into docs/Regularize3D.md - NOT committed to docs here.
ASCII; kill-table compliant; zero-oracle-wrong absolute in every adjudication.

Probe sources (scratchpad, all against libmanifold.so in the snapshot):
  probe/char_nearcoplanar.cpp - per-pair coplanarity classification + chaining
  probe/thin_gate_simplify.cpp - re-gate tolerance of thin/coincident sheets + Simplify
  probe/simplify_decisive.cpp  - Simplify on a real self-intersector

## THE PROBLEM (settled; restated, not re-derived)

Faces within eps of coplanar but not exactly coplanar are DECIDABLE at level-0
(orient3d certified nonzero, margins measured huge) but their arrangement cells
are physically THIN. B classifies them correctly (exact winding), but EMISSION
double-rounds the thin bounding triangles to slivers: the assembly eps-weld then
either merges two exact-distinct verts (re-manufactures the R1 twin / a self-fold
the re-gate's shares-vertex skip cannot see, docs/Regularize3D.md R1+R2) or
collapses the sliver to a hole. Output-REPRESENTABILITY, not decision. This is the
strict-FP sliver root cause, the wall-A twin anatomy, and the killed weld-bump,
re-encountered inside candidate B rather than the v3 sweep.

## WHERE THE HOLE ACTUALLY LIVES (measurement, corrects the brief's premise)

The brief pointed at "GT7863's near-coplanar residue." Measured (char probe):
  - GT7863 (456-tri composed soup): 438 EXACTLY-coplanar overlap pairs
    (max vertex-plane distance == 0.0 bitwise), ZERO near-coplanar pairs.
    => GT7863 is stage-4's EXACT-coplanar FOLD target, NOT a near-coplanar
    carrier. Its fail-closed is exact-zero-tie SoS + entanglement, not the thin
    band. (Honest correction, not a failure to probe it.)
  - hull (body+mask, 13440-tri soup): 14 DECIDABLE-thin near-coplanar OVERLAP
    pairs (orient3d certifies non-coplanar; gap 5.5e-5 .. 3.6e-3 eps, i.e.
    2e-14 .. 1.3e-12 absolute at eps=3.5e-10) + 2 filter-uncertain. THIS is the
    genuine near-coplanar hole geometry: two large flat facets within eps of
    coplanar.
  - siA / siB: transversal self-crossings (S3 dihedrals ~1e-1), ZERO
    near-coplanar - which is why they RESOLVE through B today (~1.5s).

REACHABILITY (the honest scope). B only runs on a DIRTY (self-intersecting)
single component. Of the corpus: siA/siB reach B and are transversal (no hole);
hull / GT7081 / Offsets DECOMPOSE into clean components and early-exit (never
reach B); GT7863-as-one-soup reaches B and is EXACT-coplanar (fold). So NO landed
single-dirty-component carrier currently triggers the near-coplanar-in-B hole. It
is a SPECIFIED-but-untriggered class (like the razor band, SweepEmit3D [D3-R5]),
that goes live when B runs on an imported/composed self-overlapping shell whose
internal overlap is near-coplanar. The hull proves the GEOMETRY is real; stage-5
must CONSTRUCT the red-first single-shell carrier.

## CANDIDATE (a) INPUT-SIDE PLANARIZATION - RECOMMENDED

REFRAME (the key simplification): (a) is not a new pre-pass. It is WIDENING the
stage-4 exact-coplanar fold's cluster-detection threshold from
"all six cross orient3d filter-signs == 0" (gap < ~1 ULP) to "max vertex-plane
distance < eps" (the near band), gated by a GLOBAL-PLANARITY guard. The fold
machinery (DetectCoplanarClusters + FoldCoplanarClusters: project members to a
shared plane -> RemoveOverlaps2D arrangement -> per-cell REAL 3D winding classify
-> emit at canonical 3D -> eps-weld) already exists and is GREEN on synthetic
exact-coplanar carriers (Coplanar_SameOriented/ThreeFaceGroup/MixedOrientation/
EpsChain/InvertedStacking all pass in the snapshot). Near-coplanar reuses it
verbatim; the projection onto the shared plane becomes a <= eps input
perturbation instead of a <1-ULP no-op.

WHY IT DIFFERS FROM THE KILLED SNAPS (variant E; identity-extension Voronoi snap).
Those snapped an EMITTED / constructed vertex position at EMISSION time, DOWNSTREAM
of the arrangement, so they fought decisions the arrangement had already made:
merging two twin images loses the lump, leaving them distinct reopens the fan, and
placement is not classification (ExactArrangement3D variant-E kill). (a) perturbs
the INPUT plane BEFORE B's enumeration / coupled winding / cell-classify / emit, so
B RE-DERIVES the entire arrangement from the snapped input - the thin cell CEASES
TO EXIST (both sheets fold to one sheet). Coordinated-by-construction, not a
downstream snap. It re-attempts NO kill in ExactArrangement3D's table: not the
weld-bump (no constant-radius weld), not variant E (no emission-time placement),
not the naked snap (collapse happens WITH re-emission, from re-derived input).

EPS-VALIDITY + winding self-check (why the answer is right).
  - The snap moves cluster verts by <= globalPlanarDev <= eps (guarded below),
    inside the standing epsilon-valid contract.
  - Winding stays over the ORIGINAL un-snapped soup (WindingAt, mult-correct). The
    fold classifies each folded cell at cen +/- eps*n. Because the true wedge is
    < eps thick, the eps-nudge STRADDLES both original sheets, so w_above/w_below
    read the OUTER strata and m == w_below - w_above holds by construction (the
    stage-4 self-check survives). The thin doubly/edge-covered wedge collapses to
    the correct single-cover {w_S>=1} boundary within eps. (Reasoned; the eps
    nudge >> globalPlanarDev makes the filter certify the crossings.)

RISKS ATTACKED (probe evidence).
  1. TRANSITIVE CHAINING / CURVATURE (the wall-B hull hazard). char-probe
     union-find over near-coplanar+overlapping pairs at K in {1,2,5,10}*eps on the
     hull: ONE multi-face near-cluster, 15 faces, spatial extent ~2.6e11 eps
     (~90 units, half the model) - BUT its GLOBAL planarity deviation (max vertex
     to area-weighted fitted plane) is only 3e-4 eps. So it is a genuinely FLAT
     facet pair, not a curved chain; the chaining hazard did NOT materialize here,
     and DIAMETER is a red herring. The decisive, cheap GUARD is GLOBAL PLANARITY:
     fit a plane to the cluster, fail closed if any vertex deviates > eps. Hull
     passes with 3000x headroom; a curved near-tangent tessellation (dev > eps)
     fails closed honestly (strictly narrower than today's blanket refusal).
  2. SNAP-INDUCED NEW INTERSECTIONS. Bounded by globalPlanarDev <= eps: any face a
     snapped vertex newly crosses was already within eps of it, so the crossing is
     sub-eps and absorbed by BuildImpl's eps-weld (the SweepEmit3D [D3-R5]
     delegate-to-vert-merge precedent). No new macro intersection can appear.
  3. VERTEX CONSISTENCY (a riser wall sharing verts with a cluster face). The fold
     emits cluster sub-faces at projected positions and risers at original ones;
     the mismatch <= globalPlanarDev <= eps is welded shut by the same eps grid.
  4. The 2 filter-uncertain hull pairs currently raise boundaryTouch (SoS
     fail-closed); widening the fold to the near band SUBSUMES them (they become
     folded, not refused).

VERDICT (a): ADVANCE as the stage-5 front-runner. Reuses a landing fold; the one
new element is a global-planarity guard that is both the curvature safety net and
cheap to compute.

## CANDIDATE (b) EXACT-SPACE THIN-CELL MERGE - KILL (dominated by a)

Collapsing cells thinner than eps in the exact 3D cell complex requires BUILDING
that global 3D arrangement cell complex first. B deliberately does NOT build one:
reg3d-s2r2 established THE BUILD dissolves into PER-FACE 2D arrangements
(RemoveOverlaps2D per crossed face) + a real 3D coupled-winding classify, never a
global 3D cell complex (the S-mult kill + S-cells pivot). (b) reintroduces exactly
the machinery B was designed to avoid, for the SAME eps-validity payoff (a) gets
by reusing the per-face fold. Same outcome, much larger surface. KILL - it
relocates (a)'s effect into a global-3D operator with no added correctness.

## CANDIDATE (c) EPSILON-VALID EMISSION (coincident-but-separate sheets) - KILL

Probe (thin_gate_simplify): a thin CLOSED box (top+bottom sheets t apart) is a
valid 2-manifold AND passes IsSelfIntersecting (selfInt=0) at EVERY thickness
t in {2, 1, 0.5, 0.1, 0.01, 0}*eps, including t=0 (coincident sheets). Two
readings, both against (c):
  - The re-gate TOLERATES the thin/coincident posture. That does not help (c); it
    HURTS it - a WRONG thin emission would pass the re-gate silently
    (oracle-wrong, violating zero-oracle-wrong-absolute). The gate cannot be the
    backstop for thin geometry (this is R1/R2 restated empirically).
  - (c) also emits the WRONG boundary. Two near-coplanar overlapping sheets are
    INTERIOR to {w_S>=1} (w>=1 on both sides), so B's threshold read DISCARDS them
    (ExactArrangement3D R3-i: interior 1|2, 2|1 sheets are not boundaries).
    Emitting them as coincident-but-separate manufactures a thin void / retains a
    double cover = not the {w_S>=1} boundary. And it never touches the R1 weld
    corruption that is the actual failure. Emit-and-accept is not regularization.
KILL.

## CANDIDATE (d) POST-PASS SIMPLIFY - KILL as a fix (optional cosmetic only)

Decisive probe (simplify_decisive): siA.Simplify(eps) returns the mesh UNCHANGED
- tris 17160 -> 17160, vol 0.147953 -> 0.147953, IsSelfIntersecting 1 -> 1.
Simplify is topology-preserving DECIMATION; it does not resolve self-overlap or
winding, so it CANNOT stand in for B. As a POST-B pass it is a blunt GLOBAL
operator: thin_gate_simplify shows Simplify(eps) collapses a thin box to vol=0 for
t <= 0.1 eps while preserving t >= 0.5 eps - i.e. it collapses ANY sub-tolerance
thin feature indiscriminately and cannot tell a sliver ARTIFACT from a thin
FEATURE. The identity-extension arc already measured the consequence: dropping
load-bearing slivers OPENS MACRO HOLES (len ~70). boolean2's Simplify-at-tolerance
precedent is 2D; on 3D slivers it is manifold-safe on ISOLATED slivers but
hole-prone on load-bearing ones and perturbs clean geometry globally. KILL as the
fix; keep only as an OPTIONAL caller-side cosmetic outside the zero-oracle-wrong
contract.

## RECOMMENDATION (sized for stage-5)

1. DetectCoplanarClusters: relax `coplanar(i,j)` from all-six-orient3d-filter-0 to
   maxVertexPlaneDistance(i,j) < eps AND overlap2D (already present). This routes
   the near band to the fold instead of to seams.
2. GLOBAL-PLANARITY GUARD (new, the only net-new element): per detected cluster,
   fit a plane (centroid + area-weighted normal); if any member vertex deviates
   > eps, FAIL CLOSED with a distinct named reason (curved near-coplanar chain,
   not the transversal SoS residue). This is the curvature/extent guard; it turns
   the wall-B chaining hazard into an honest, strictly-narrower refusal.
3. FoldCoplanarClusters is reused verbatim: project to the fitted plane, arrange,
   classify by REAL 3D winding over the ORIGINAL soup, emit at canonical 3D,
   eps-weld. Keep the m == winding-jump self-check (holds via the eps-nudge).
4. RED-FIRST carriers stage-5 must build: (i) a near-coplanar dirty single-shell
   (two near-coplanar overlapping plates tilted sub-eps, closed into a
   self-intersecting shell, made gate-dirty by a transversal riser) -> resolve
   oracle-true (MC/GWN volume band + one component + tol-invariance); (ii) a
   CURVED near-coplanar chain (near-tangent tessellated spheres) -> global-
   planarity guard fail-closed; (iii) mult-1 EXACT-coplanar pin stays BITWISE
   (the near widening must not perturb stage-4's exact path).
5. Land as a THRESHOLD DELTA on top of stage-4's fold, after stage-4 merges, to
   avoid churn - coordinate so the exact path is unchanged and the near band is
   purely additive.

## HONEST RESIDUE

- No current single-dirty-component carrier triggers the hole; stage-5 owns
  constructing the red-first shell. The hull proves the geometry but decomposes
  clean, so it is a stress fixture (one-soup), not a natural B input.
- The guard threshold (exactly eps vs a fraction) needs a synthetic sweep to tune:
  the hull's 3e-4-eps dev gives huge headroom, but the razor band (2-10 eps,
  SweepEmit3D [D3-R5], recorded UNTESTED) is the boundary case that decides
  fold-vs-fail-closed. Too tight re-refuses valid flat facets; too loose admits a
  chain that snaps eps-invalidly.
- ENTANGLEMENT unchanged: a near-coplanar cluster face ALSO transversally pierced
  by a non-coplanar face stays the stage-4 entanglement decline (fail closed).
- SUBTRACTION / negative winding under a near-coplanar fold is untested, exactly as
  it is for the exact fold today.
- The winding self-check holding via the eps-nudge is REASONED, not yet run
  end-to-end (blocked on the red-first carrier); it is the first thing stage-5's
  fixture (i) should empirically confirm.
