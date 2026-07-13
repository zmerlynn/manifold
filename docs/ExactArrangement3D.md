# Exact-Arithmetic Arrangement for the Hard Carriers (design)

Status: DESIGN, load-bearing-probe-validated core, one named open risk. This is
the campaign's final open question made concrete. It is grounded on three
converged results (do not re-derive them here):

  1. The hard carriers' correct topology is a PURE FUNCTION of the input soup:
     the oracle union is {p : w_S(p) >= 1}, w_S the signed winding of the
     oriented input faces. No ambiguity witness exists
     (.claude/lane-reports/a0-verify-witness-...).
  2. Every eps-quantized channel is DEGENERATE on those carriers (V4ImplPlan.md
     stage-A0 RECORDED RESULT + verification: incidence constant, fold bijective,
     positions non-monotone).
  3. Therefore THE BARRIER IS PRECISION OF DECISIONS, NOT INFORMATION.

"Exact arithmetic" here means EXACT-VALUED DECISIONS, not infinite-precision
storage: predicates on INPUT data (doubles are exact rationals; orient3d and
winding are evaluable exactly by adaptive precision or rationals), with
constructions either avoided (decisions re-derived from input predicates) or
represented exactly where unavoidable. Output COORDINATES may stay double /
eps-noisy; the TOPOLOGY must be exact. eps remains the honest weld radius for
constructed output geometry.

This doc is the exact-decision counterpart to docs/MaintainedEmission3D.md (which
converged on variant iii = the RSI-#3 arrangement completion at EPS precision,
and named exactness as the upgrade path). It re-adjudicates the one FP-premised
RSI-#3 kill under exactness, and adjudicates where exact-valued decisions enter
the pipeline.


## The load-bearing probe (run first; it decides whether the arc is live)

Question: is the material the escape kill LOST present in the exact input, or is
the input genuinely ambiguous there? If w_S = 0 at the lost lump, the whole arc
is dead. Exact rational winding (Python/Fractions, segment-crossing on input
verts, throwaway driver; exact query points recorded in the notebook) on
Havocglass8's composed soup S = faces(left) ++ faces(right):

  - Faithfulness gates PASS. Both operands are closed oriented 2-manifolds
    (every directed edge paired) with positive signed volume, so orientation is
    preserved and OBJ-direct winding is faithful to the pipeline's soup (winding
    is a topological invariant of the oriented surface; ComposeImpl/GetMeshGL64
    merges nothing here - no intra-operand near-duplicate verts).

  - At the escape kill's lost-lump coordinate: w_S = 2 (inside both operands).
    MATERIAL. The oracle keeps it.
  - Correctly-retained control region: w_S = 1. MATERIAL, kept.
  - Outside control: w_S = 0.
  - The right operand pokes a large wedge OUT of the left: most of its vertices
    carry right-only material (w_left = 0, w_right = 1, w_S = 1). The lost ~336
    lump is a slice of this wedge.

VERDICT: the theory holds emphatically. The lost region and the correctly-
retained regions carry the SAME KIND of w_S value (both >= 1, both kept). eps
cannot tell them apart from quantized emitted positions; exact w_S treats them
identically. The 336 delta is a pure TOPOLOGY-SELECTION error at eps, which exact
w_S >= 1 decides trivially. The decision A0 proved unreachable at eps is made
correctly from input data alone. THE ARC IS LIVE.


## The degeneracy is emission-manufactured, not intrinsic (the pivot)

Measured on the same soup: the input has NO sub-eps degeneracy. The closest
left/right vertex pairs are two parallel sheets offset by ~7e-4 (tens of
thousands of eps); every right vertex sits tens of eps or more from the nearest
left face plane. The sub-eps "twin" (~1.7 eps cap images) that the wall-A arc
names is a SWEEP-EMISSION ARTIFACT: steep-track per-plane self-location amplifies
the modest input offset into sub-eps emitted cap images, which then must collapse
(lose the lump) or not (reopen the fan).

Exact arithmetic HAS NO EPS. It resolves the input offset and the vertex/plane
gaps with certainty and never manufactures the twin, so it never faces the
sub-eps collapse decision. Sampling exact w_S up a column through the near-
coincident sheets stratifies the region into clean integer shells - crossing
upward, w_S steps 2 -> 1 -> 0, each step exactly one sheet's signed multiplicity.
The thin right-only sliver BETWEEN the sheets (w_S = 1) is the lump; the
{w_S >= 1} boundary is the single outer sheet, well-defined. This is manifold-by-
construction, locally, on the real geometry: no collapse decision, both sheets
kept distinct, the sliver survives. The eps-collapse of the manufactured twin is
what deletes it.

Consequence: the hard carrier's near-degeneracy is a property of the eps sweep-
emission MACHINERY, not of the input. This is the strongest structural argument
for an exact arrangement - it sidesteps the amplification entirely.


## The reopen (the one sanctioned re-adjudication)

The June RSI-#3 memo killed option 3 - "complete the arrangement + classify each
sub-face independently by a 3D interior-point winding probe = manifold-by-
construction" - as FALSE: independent FP raycasts from different origins have no
cross-origin symbolic-perturbation coupling, so adjacent sub-faces across a
shared edge can disagree, leaving an unpaired halfedge. The kill is CORRECT under
FP.

Under EXACT evaluation with ONE global symbolic-perturbation convention the
premise is removed. All probes evaluate the SAME exact function w_S:
  - Each cell's classification is a well-defined integer; two adjacent cells
    cannot disagree about w_S.
  - The kept boundary is the topological boundary of the solid {w_S >= 1}, a
    closed surface whose halfedges pair by construction. Around every arrangement
    edge the winding jumps sum to zero exactly (conservation of the signed
    multiplicities of the input faces meeting there), so the fill transitions
    pair up. This is the cross-edge coupling the DCEL provided for free, now
    provided by exactness.

REOPEN VERDICT: the FP-disagreement kill does NOT survive exact winding with a
global convention. It does not, however, resurrect a shortcut. Exactness cures
the CLASSIFICATION coupling; "manifold by construction" still requires the
COMPLETE EXACT ARRANGEMENT (all intersection curves constructed, all faces
split). Corrected option-3 therefore collapses INTO the full exact arrangement
(variant B below). The kill's mechanism dies; its deeper wisdom - you cannot skip
the arrangement - survives.

SCOPE. Only the FP-raycast-disagreement kill is reopened. The memo's other kills
stand and are NOT FP-disagreement kills: reusing boolean2's per-plane 2D winding
(the per-plane 2D winding is not the 3D solid winding - wrong quantity, not an FP
tie), and the serial seam-assembly matcher (unspecified triple-point rule). Every
measured kill in the campaign (the seven wall-A kills, the escape kill's 336
invariance, A0, the witness theorem, the naked snap) stands.


## In-tree machinery inventory (what an exact design can and cannot reuse)

Grepped src/: there is NO exact-predicate kernel. No orient3d/orient2d, no
adaptive/Shewchuk expansions, no interval arithmetic, no rationals, no int128.
The robustness surface is:
  - shared.h Shadows(p,q,dir): a 1D FP predicate that breaks EXACT ties by a
    `dir` sign - manifold's entire "one global perturbation convention", used
    pervasively in boolean3.cpp's sweep.
  - eps = kPrecision * bBox.Scale() welding; FP shoelace areas; FP AABB collider.

So an exact-decision architecture must INTRODUCE adaptive orient3d + exact winding
(new surface area, not a reuse). It CAN and MUST reuse the Shadows PATTERN - one
global symbolic-perturbation convention - which is exactly what the reopen
requires so all probes agree at shared boundaries. This is a design constraint,
not a free lunch.


## Variant space and verdicts

Adjudicated against: correctness composition (the relocation test - does an eps
decision downstream undo the exact one?), cost (should scale with DEGENERACY, not
input size), integration (what of v3 survives), and the reopen.

  (A) EXACT DECISION CORE in the existing eps pipeline. The eps pipeline runs as
      today; the named decisions (fan pairing, survive-vs-collapse, cancelled-edge
      agreement) are re-posed as exact input predicates consumed at the existing
      sites.
      VERDICT: KILL - relocates. The load-bearing survive-vs-collapse call is made
      at EMISSION (EmitCaps weld on constructed twin positions), DOWNSTREAM of the
      arr.vert sites where an exact answer would be injected. The eps weld undoes
      the exact answer (merge = lose lump / keep = reopen fan). To reach the
      deciding site you must exactify emission, which is C or B. Subsumed.

  (B) FULL EXACT ARRANGEMENT. Exact intersection points (rational constructions),
      complete arrangement, exact per-cell winding classification (the reopened
      option, corrected), extract the {w_S >= 1} boundary, round outputs to double
      with a validity gate.
      VERDICT: ADVANCE - correct, expensive. Topology all exact; the only eps
      decision is the output-coordinate snap, gated by Is2Manifold / self-
      intersection (snap-rounding literature applies; manifold's eps-weld already
      IS a snap-round). Passes the relocation test for topology. STRIKE: cost
      scales with INPUT size and arrangement complexity, not with degeneracy.
      Integration: DIES most of v3; needs a from-scratch exact kernel. Role: the
      correctness backstop / north star, not the recommended first build.

  (C) LOCAL EXACT RESOLVER (hybrid). The eps pipeline runs; where it fails closed
      (the flagged junction clusters / guards - the flagging machinery already
      exists as FatalReason returns at BuildSlabs / EmitCaps / EmitStrips), an
      exact local resolver takes ONLY the faces incident to the flagged region,
      computes their exact local arrangement + winding, emits the correct local
      patch, and welds to the eps exterior at a boundary in clean territory.
      VERDICT: ADVANCE - FRONT-RUNNER. This is the escape kill's SHAPE WITH the
      information supplied. The escape kill failed because it was positions-only:
      it collapsed clusters on emitted positions and had to guess collapse-or-not,
      an invariant-336 topology-selection error the emitted positions could not
      arbitrate. Exact w_S arbitrates it (probe: the lump is material; the sheets
      stratify cleanly). The kill was information-starved, not mechanism-broken.
      Cost scales with DEGENERACY (the bulk runs eps; only flagged submeshes get
      the exact treatment = a BOUNDED B). Integration: MAX v3 reuse. RISK: the
      patch/exterior weld (named below); mitigated in principle because near-
      degeneracies are localized (probe: the near-coincidences sit in one sheet
      region; elsewhere the operands are macroscopically separated) so a guard-
      clean, input-edge boundary exists. C = "bounded B welded to the eps
      exterior."

  (D) EXACT 2D SECTIONS. Per-slab sections computed as exact functions of input
      for flagged slabs.
      VERDICT: KILL - insufficient. Exactifies the in-plane arrangement but leaves
      the CROSS-PLANE emission (junction self-location across cap planes) in eps.
      The twin is a cross-plane phenomenon - two cap planes at distinct criticals
      emit one junction twice; exact-per-slab does not unify the cap images. The
      twins SURVIVE exactness-per-slab. It does not reach the deciding decision.

RECOMMENDATION: build C (local exact resolver = bounded exact arrangement + clean-
boundary weld), with B as the correctness backstop for any region C cannot bound
with a clean boundary. C and B are the same object at two scopes; C is B made to
pay only for degeneracy.


## Surviving architecture (C): staged sketch

Stage 1 - EXACT PREDICATE KERNEL (the reusable primitive; cost scales with
  degeneracy, not size).
  - Adaptive orient3d (FP filter first; exact expansion only when the FP
    determinant is not sign-certain). Cheap in the common non-degenerate case.
  - Exact w_S(p): segment-crossing count on INPUT vertices (orient3d signs), one
    fixed global symbolic-perturbation convention (the Shadows pattern lifted to
    3D orient) so a ray grazing a shared edge/vertex resolves identically for
    every probe. This is the load-bearing consistency guarantee of the reopen.
  - Pin: a differential test vs the throwaway rational driver on the probe's
    query points (they must agree bit-for-bit on sign).

Stage 2 - FLAG + BOUND. The eps pipeline runs unchanged. On a guard fire (or a
  new near-degeneracy flag), collect the incident INPUT faces and GROW the region
  until its boundary is guard-clean: every boundary edge is an INPUT edge that the
  eps pipeline resolved without a guard (transverse, eps-robust). The region is a
  submesh whose interface with the exterior is shared input geometry (bitwise
  equal on both sides).

Stage 3 - LOCAL EXACT ARRANGEMENT. Compute the exact arrangement of the submesh:
  rational intersection points among its faces, cell decomposition, classify each
  half-face by exact w_S (Stage 1), extract the oriented {w_S >= 1} boundary. The
  topology is exact; constructed intersection coordinates are rounded to double
  only at emission, behind a LOCAL validity gate (2-manifold, no self-intersection
  within the patch).

Stage 4 - WELD. Stitch the local patch to the eps exterior along the shared
  input-edge boundary (boundary vertices are input vertices, shared exactly).
  Measure-zero contacts (touching sheets) are handled by the existing
  SplitTouchingSheets, as today.

Stage 5 - VALIDITY GATE + FALLBACK. If the welded result has unpaired halfedges
  or fails validity (it should not, if the boundary is clean), fall back to B on
  that connected component, or fail closed. The honest fail-closed contract is
  preserved as the floor; the exact resolver only ever UPGRADES a fail-closed
  region to correct.

Second-probe status: the local mechanism's CORE is validated at toy scale on the
real cluster (exact winding stratifies the sheet region cleanly; the lump
survives; the boundary is manifold-by-construction). The UNVALIDATED load-bearing
assumption is Stage 2's clean-boundary growth at corpus scale (see Risks).


## Risks (the attack surface for the adversarial round)

R1. CLEAN-BOUNDARY GROWTH (the load-bearing assumption). C's cost and correctness
    both hinge on being able to grow every flagged region to a guard-clean, input-
    edge boundary. Attack: a carrier where near-degeneracies PERCOLATE - no clean
    boundary exists within a bounded neighborhood, so the flagged region grows to
    the whole mesh and C degenerates to B (input-scaling cost) or the boundary is
    forced through near-degenerate territory (weld relocation). Needed: measure
    flagged-region growth on the corpus; is it bounded, and by what?

R2. WELD RELOCATION. Even with a shared input-edge boundary, a boundary input edge
    may be crossed by a third face near-degenerately, so the eps exterior
    subdivides it at an eps-placed vertex the exact patch places exactly -> an
    unpaired halfedge at the weld. This is the relocation test applied to the
    weld: the near-degeneracy decision must not reappear in the stitching
    predicate. Mitigation is R1's clean-boundary property; the attack is a
    boundary edge with a hidden near-degenerate crossing.

R3. OUTPUT SNAP-ROUNDING. Rounding constructed intersection coordinates to double
    can create NEW near-coincidences at output (two exact-distinct patch vertices
    round together), silently breaking 2-manifoldness below the validity gate's
    resolution. Iterated snap-rounding may be required. Attack: a patch whose
    exact vertices are closer than eps after the arrangement, so no snap is both
    valid and faithful.

R4. ONE CONVENTION ACROSS ALL PREDICATES (the subtle reopen residue). Exactness
    cures the FP-disagreement kill ONLY IF a SINGLE symbolic-perturbation
    convention governs EVERY predicate family (orient3d, the winding segment-
    crossing, the arrangement's intersection ordering). A convention consistent
    within orient3d but different for the winding ray reintroduces cross-boundary
    disagreement - the kill in disguise. Attack: construct a mixed-convention gap.

R5. SUBTRACTION / NEGATIVE-MULT REGIME. Every measurement is union-composed (all
    multiplicities +1). Subtraction inputs carry negative multiplicities; w_S can
    go negative; the {w_S >= 1} rule still holds but the strata and the winding
    jumps differ. Unmeasured. Attack: a subtraction carrier with an inward shell.

R6. COINCIDENT / COPLANAR FACES (measure-zero overlaps). Exact w_S is undefined ON
    a coincident face (measure zero); the fill needs a tie-break there (the
    coincident-wall cleaning the witness lane's Case B relied on). Attack:
    coplanar coincident sheets where w_S jumps by 2 across a doubled wall, or an
    operand wall buried inside another - the tie-break must match the library's.

R7. KERNEL COST + CERTIFICATION. Adaptive orient3d must be certified correct
    (Shewchuk-class); the exact arrangement's rational constructions are not cheap
    even when bounded. Attack: a flagged region that is bounded but internally
    dense (many mutually near-degenerate faces), so the local arrangement is
    combinatorially large despite a small face count.


## What is settled (do not re-litigate)

The probe results and every measured kill across the campaign; A0; the witness
theorem; the escape kill's 336 invariance. The v3/v4 branches stay as they are;
zero-oracle-wrong; correct-or-honest. The only object reopened is the FP-premised
RSI-#3 option-3 kill, which does not survive exact winding with one global
convention - and which, corrected, is the exact arrangement this doc designs C to
pay for only where the eps pipeline fails closed.
