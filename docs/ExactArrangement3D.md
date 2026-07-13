# Exact-Arithmetic Arrangement for the Hard Carriers (design)

Status: DESIGN, load-bearing probe validated by two independent algorithms,
round-1 adversarial review folded (rev 2). The open work is now NAMED, not
hand-waved: a per-region localizer, a selective weld, and clean-boundary growth
measurement (all stage-0, measurable without the exact kernel). This is the
campaign's final open question made concrete. It is grounded on three converged
results (do not re-derive them here):

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

Round-1 review lanes folded here (exact numbers live in the notebooks; this doc
carries magnitudes and cites):
  - probe : .claude/lane-reports/v5-verify-probe-1783913337.md (probe reproduction,
            per-carrier pivot, integer shells, query cost)
  - c     : .claude/lane-reports/v5-verify-c-1783911848.md     (localizer, selective
            weld, downstream care, per-carrier scope, zero silent wrong-resolve)
  - bvc   : .claude/lane-reports/v5-verify-bvc-1783911600.md   (variant E kill,
            B/C merge, reopen phrasing, staged-sketch stage-0)


## The load-bearing probe (run first; it decides whether the arc is live)

Question: is the material the escape kill LOST present in the exact input, or is
the input genuinely ambiguous there? If w_S = 0 at the lost lump, the whole arc
is dead. Exact rational winding on Havocglass8's composed soup S = faces(left) ++
faces(right):

  - Faithfulness gates PASS. Both operands are closed oriented 2-manifolds
    (every directed edge paired) with positive signed volume, so orientation is
    preserved and OBJ-direct winding is faithful to the pipeline's soup (winding
    is a topological invariant of the oriented surface; ComposeImpl/GetMeshGL64
    merges nothing here - no intra-operand near-duplicate verts).

  - At the escape kill's lost-lump coordinate: w_left = 1, w_right = 1, w_S = 2
    (inside both operands). MATERIAL. The oracle keeps it.
  - Correctly-retained control region: w_S = 1. MATERIAL, kept.
  - Outside control: w_S = 0.
  - The right operand pokes a large wedge OUT of the left: most of its vertices
    carry right-only material (w_left = 0, w_right = 1, w_S = 1). The lost lump
    is a slice of this wedge.

INDEPENDENTLY CONFIRMED (round-1, probe lane). Two genuinely different algorithms
reproduced every value bit-for-bit against the original segment-cast driver:
exact solid-angle summation (Van Oosterom-Strackee, high-precision), and exact
rational semi-infinite ray casting over nine independent generic directions
(unanimity across directions is an internal ray-invariance proof). Margins are
astronomically clean (solid-angle residuals from the nearest integer are ~1e-100;
the rays are unanimous). INTERIORITY is exact: these are not boundary-adjacent
points - the lump sits millions of eps deep inside both surfaces, the control
centroids hundreds of millions of eps deep. "Strictly inside" is an exact, robust
statement here, not a sampling artifact.

VERDICT: the theory holds emphatically. The lost region and the correctly-
retained regions carry the SAME KIND of w_S value (both >= 1, both kept). eps
cannot tell them apart from quantized emitted positions; exact w_S treats them
identically. The lost lump is a pure TOPOLOGY-SELECTION error at eps, which exact
w_S >= 1 decides trivially. The decision A0 proved unreachable at eps is made
correctly from input data alone. THE ARC IS LIVE.

INTEGER SHELLS - caveat (probe lane, ATTACK 3). Sampling exact w_S up a column
through the near-coincident sheets stratifies the region into clean integer
shells, and the thin right-only sliver survives - the architectural conclusion is
correct. But the JUSTIFICATION is exact evaluation, NOT the sampling. The shell
sequence is location-dependent (0 -> 1 -> 0, 2 -> 1 -> 0, and 0 -> 1 -> 2 all
occur along different columns) with thicknesses from tens of eps (the edge-edge
surface minimum) to over ten thousand eps. Critically, at the mission-lump column
a shell hundreds of eps thick was SKIPPED by both a 400-sample and a 700-sample
scan (both reported a clean 0 -> 2 jump); only bisecting the two sheet crossings
separately resolved it. A thin shell can hide between samples. So the manifold-by-
construction claim rests on the winding being a clean integer OFF-surface
(verified exactly, margins ~1e-80..1e-100), not on any stratification sample.

QUERY COST - measured (probe lane, ATTACK 4). Exact w_S is O(ntri) PER QUERY. The
Python/high-precision baseline runs in tens of ms per query at Havoc scale (176
tris), which projects to seconds per query at GT7081 scale (~31k tris). Those
absolute numbers are Python-specific; a C++ adaptive-orient3d kernel (FP filter,
exact only when the sign is uncertain) is the production path and would be orders
of magnitude faster in the common case. The O(ntri) full-soup scaling is REAL and
makes the B-scales-with-input / C-bounds-to-submesh adjudication concrete. Ray
degeneracy (a grazing ray needing symbolic perturbation) is BOUNDARY-ONLY: zero
occurrences in thousands of random interior queries, but 100% on points placed
exactly on a face. Because the design classifies CELLS by interior points, the
symbolic-perturbation machinery is almost never invoked at runtime - the SoS load
is light (R4's correctness concern still stands; it is just not a hotspot).


## The pivot, re-scoped per carrier (was: emission-manufactured, not intrinsic)

The rev-1 pivot claimed the hard carriers' near-degeneracy is a property of the
eps sweep-emission MACHINERY, not of the input, and offered that as the general
structural argument for an exact arrangement. Round-1 direct measurement of the
minimum INPUT feature separation per carrier REFUTES the universal claim. The
pivot is HAVOC-ONLY:

  carrier    input degeneracy                                  verdict
  -------    ----------------------------------------------    ------------
  Havoc      no sub-eps input feature; surface min ~50 eps     MANUFACTURED
             (edge-edge), so the ~1.7-eps emitted twin is
             genuinely made by emission
  GT7863     a vertex exactly ON a face (zero separation,      INTRINSIC
             exact coplanar + inside)
  Offset1    exactly-coincident vertices (dx = dy = dz = 0)    INTRINSIC
  openscad   exactly-coincident vertex pairs                   INTRINSIC
  GT7081     over a million near-tangent seam contacts +       INTRINSIC
             hundreds of near-coplanar overlapping faces       (surface
             (surface near-tangency)                            near-tangency)

For HAVOC the pivot holds, but the rev-1 cleanliness numbers were OVERSTATED:
  - "parallel sheets ~7e-4 apart (tens of thousands of eps)" was the VERTEX-vertex
    gap. The true SURFACE-surface minimum is ~50 eps (edge-edge) - hundreds of
    times closer. Emission amplifies ~50 eps into the ~1.7-eps twin (a modest
    factor), which is plausible; the "tens of thousands of eps" headline was a
    vertex artifact.
  - "right verts >= 23 eps from the nearest left face plane" was a weak proxy:
    plane distance, not face distance (a vert can sit tens of eps from a far
    triangle's infinite plane while being tens of thousands of eps from any actual
    face). The 23-eps figure is not evidence of near-tangency and is dropped.

For Havoc, exact arithmetic still sidesteps the amplification: it resolves the
~50-eps input offset with certainty, never manufactures the ~1.7-eps twin, and so
never faces the sub-eps collapse decision. Both sheets are kept distinct, the
right-only sliver (the lump) survives, the {w_S >= 1} boundary is the single outer
sheet. The eps-collapse of the manufactured twin is what deletes it.

ARCHITECTURAL CONSEQUENCE (the narrowing STRENGTHENS C). The four intrinsic
carriers carry exact-coincident vertices, a vertex exactly on a face, or genuine
sub-eps surface contacts AT THE INPUT. That is precisely what a local exact
arrangement handles and what no placement shortcut ever could - you cannot "place"
your way out of two vertices that are exactly equal, or a vertex that is exactly
on a face. So the pivot's Havoc-only rationale does NOT weaken the arc; it removes
the emission-manufactured hope as a GENERAL argument (which is what kills variant
E below) and leaves the local exact arrangement as the only mechanism that covers
both regimes - manufactured (Havoc) and intrinsic (the other four).


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
    cannot disagree about w_S. Exactness provides CONSISTENT CLASSIFICATION.
  - Around every arrangement edge the winding jumps sum to zero exactly
    (conservation of the signed multiplicities of the input faces meeting there),
    so the number of boundary sectors around each edge is EVEN. Exactness
    provides PARITY.

What exactness does NOT provide (rev-1 overclaim, corrected). The rev-1 text said
the halfedges "pair by construction... the cross-edge coupling the DCEL provided
for free, now provided by exactness." That oversells it. The June memo's durable
result is that manifold output needs a GLOBAL DCEL with RADIAL ORDERING of the
face-sides around each arrangement edge. Exactness gives consistent classification
and parity; it does NOT give the radial ORDER itself. You still must sort the
face-sides cyclically around each edge to PAIR the in/out transitions and to
DETECT genuine non-manifold edges (more than two boundary sectors). The radial-
ordering DCEL work REMAINS; exactness makes its predicates tie-free, and that
tie-freedom is exactly what removes the June FP-disagreement kill.

REOPEN VERDICT: the FP-disagreement kill does NOT survive exact winding with a
global convention. It does not resurrect a shortcut. "Manifold by construction"
still requires the COMPLETE EXACT ARRANGEMENT (all intersection curves
constructed, all faces split, radial order sorted). Corrected option-3 therefore
collapses into a COMPLETE exact arrangement - global (variant B) or local (variant
C) depending on whether a clean local boundary exists (see the B/C merge below).
The kill's mechanism dies; its deeper wisdom - you cannot skip the arrangement -
survives.

SCOPE. Only the FP-raycast-disagreement kill is reopened. The memo's other kills
stand and are NOT FP-disagreement kills: reusing boolean2's per-plane 2D winding
(the per-plane 2D winding is not the 3D solid winding - wrong quantity, not an FP
tie), and the serial seam-assembly matcher (unspecified triple-point rule). Every
measured kill in the campaign (the seven wall-A kills, the escape kill's lump
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
not a free lunch. (Round-1 grep confirmed: no adaptive/exact kernel in tree.)


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
      deciding site you must exactify emission, which is the exact arrangement.
      Subsumed.

  (B/C) EXACT ARRANGEMENT WITH A SCOPE DIAL. One mechanism - exact intersection
      points, complete arrangement, radial-order DCEL, exact per-cell winding
      classification, extract the {w_S >= 1} boundary, round outputs to double
      behind a validity gate - run at one of two SCOPES. These are NOT rival
      architectures (rev-1 presented them as two variants; that overstated the
      distinction). The scope dial:
        - C = flagged-region PRODUCTION mode. The eps pipeline runs; where it
          fails closed, the exact arrangement runs on ONLY the faces incident to
          the flagged region and welds to the eps exterior at a clean boundary.
          Cost scales with DEGENERACY. This is the escape kill's shape with the
          information supplied: the escape kill failed because it was positions-
          only (it guessed collapse-or-not on emitted positions, an invariant
          topology-selection error the positions could not arbitrate); exact w_S
          arbitrates it (probe: the lump is material, the sheets stratify).
        - B = whole-mesh / offline-oracle mode. The same mechanism over the whole
          soup. Cost scales with INPUT size and arrangement complexity, not with
          degeneracy. Role: the correctness backstop (run offline to validate C),
          and the unbounded LIMIT of C when no clean boundary exists.
      VERDICT: ADVANCE - front-runner is the C scope, with the B scope as oracle
      and fallback. GT7081's anatomy is delocalized (over 90% of critical-x gaps
      are sub-eps, one contiguous run of ~34k near-coincident planes packed into a
      sub-eps band), so it lands at the WHOLE-MESH end of the dial: C = B there,
      no speedup. The relocation test passes for TOPOLOGY at both scopes; the only
      eps decision is the output-coordinate snap (gated by Is2Manifold / self-
      intersection). Integration: the C scope is MAX v3 reuse (bulk runs eps); the
      B scope dies most of v3. See the C preconditions and staged sketch.

  (D) EXACT 2D SECTIONS. Per-slab sections computed as exact functions of input
      for flagged slabs.
      VERDICT: KILL - insufficient. Exactifies the in-plane arrangement but leaves
      the CROSS-PLANE emission (junction self-location across cap planes) in eps.
      The twin is a cross-plane phenomenon - two cap planes at distinct criticals
      emit one junction twice; exact-per-slab does not unify the cap images. The
      twins SURVIVE exactness-per-slab. It does not reach the deciding decision.

  (E) EXACT VERTEX PLACEMENT (no arrangement). PROPOSED-AND-KILLED (round-1, bvc
      lane). Proposal: at flagged junctions only, compute the emitted vertex
      position from INPUT data exactly, round once, place both cap images there -
      no local arrangement, no patch - so the twins never diverge because both
      derive from one exact point.
      VERDICT: KILL. Decisive argument: exact PLACEMENT of a vertex is not exact
      CLASSIFICATION of a cell. The lump's survival is a per-cell WINDING fact
      (w_S = 2, independently re-derived; column 0 -> 2 -> 0), proven INVARIANT
      across representative choice by the escape kill - an input-exact point is
      merely a THIRD representative, covered by that invariance. So on the images E
      can touch, making the two images bitwise-equal is a MERGE = the killed weld-
      bump / probe-3 collapse, which loses the lump; leaving them distinct reopens
      the fan (E built no arrangement to pair them). And on Havoc's ACTUAL blockers
      - weld images that name no canonical arr.vert (bare 2D track crossings with
      no junction identity) and far-track interpolations - placement has nothing to
      anchor to; E cannot even place them (this is R2's already-measured never-
      bound wall). It also does nothing for GT7863's bracketed-macro dead zone,
      which is a cell/winding question. E is the naked snap with cleaner
      coordinates on the images it can touch and the never-bound wall on the images
      that block. Placement is not classification. E does NOT demote C - it
      sharpens why C's local arrangement + per-cell winding is necessary.

RECOMMENDATION: build the C scope (flagged-region production: bounded exact
arrangement + clean-boundary weld). The B scope is the offline oracle that
validates it and the fallback for any region C cannot bound with a clean boundary
(percolation). C and B are the same object at two scopes; C is B made to pay only
for degeneracy.


## C preconditions (round-1 review: what must exist before build)

Round-1 (c lane) broke C on three of four attacks. The mechanism survives (exact
w_S arbitrates the lump; exact-arrangement cost is face-count-bounded corpus-wide;
no silent wrong-resolve found), but three preconditions must be built and one
scope claim corrected before C is buildable.

1. LOCALIZER (Stage 2 has no seed today). C's flag-and-grow needs to know WHERE it
   failed and which INPUT faces are incident. Neither exists. Today's flags are
   GLOBAL: NonManifoldEmission is `!IsManifold() || !Is2Manifold()` on the whole
   mesh (overlap3.cpp:1189) with no offending edge/vertex; ArrangementBudget is a
   cumulative piece counter (overlap3_sweep.cpp:246) with no region;
   SplitTouchingSheets returns a bare bool. And emitted OutTri3D is `vec3 v[3]`
   (overlap3.cpp:603) with NO provenance to input/canonical faces. So there is
   neither a SEED (which region failed) nor a MAP (emitted geometry -> input
   faces). The rev-1 line "the flagging machinery already exists as FatalReason
   returns" conflated knowing-it-failed with knowing-WHERE. Building the localizer
   = instrument emission provenance + make each gate report the offending halfedge
   + walk it back to the incident input faces. Named design work, stage-0 below.

2. SELECTIVE WELD (Stage 4 is unsound as sketched). There is ONE assembly weld:
   BuildImpl getVertIdx (overlap3.cpp:1119-1141), a uniform hash grid returning
   the min-index prior vert within Euclidean eps, through which EVERY emitted
   triangle passes - exterior AND any C patch. A faithful toy re-impl shows a
   single uniform radius CANNOT serve an exact patch:
     - two exact-distinct patch verts ~0.4 eps apart are MERGED - the arrangement
       resolved them, the weld re-manufactures the twin;
     - a patch-INTERIOR constructed vert ~0.3 eps from a shared boundary vert is
       MERGED - the interior collapses onto the boundary, exact topology corrupted.
   The two requirements collide: one radius must simultaneously MERGE patch-
   boundary onto the exterior vert (to stitch) and NOT merge sub-eps-distinct
   patch-interior verts (to preserve topology). The corpus DOES emit sub-eps-
   distinct arrangement verts (GT7081's sub-eps critical structure is direct
   evidence), so the collision is real, not hypothetical. The patch needs weld
   EXEMPTION with boundary-conformance constraints: the patch boundary adopts the
   eps exterior's vert positions (so the stitch conforms even where the exterior
   emitted constructed, eps-noisy positions along an input edge); interior verts
   are weld-exempt. A bounded local weld distinct from the global eps grid.

3. DOWNSTREAM CARE (Stage 4 walk-order). The assembly walk is: getVertIdx eps-weld
   -> degenerate/exact-dup-tri drop -> SplitTouchingSheets -> CreateHalfedges ->
   IsManifold/Is2Manifold gate -> SortGeometry/SetNormals/SetEpsilon. Two entries
   need care:
     - getVertIdx eps-weld: NEEDS-BYPASS (precondition 2). The dup-tri drop runs on
       WELDED indices, so it is coupled to the bypass.
     - SplitTouchingSheets: NEEDS-CARE. Its radial pairing returns false ->
       NonManifoldEmission on radial ties, zero-length slivers, or non-alternation.
       A legitimate exact patch with GENUINE tangent sheets (R6 touching contacts)
       is REJECTED, not composed - exact-correct patches can fail closed here.
     - The IsManifold/Is2Manifold gate is patch-safe but GLOBAL: one bad weld
       anywhere aborts the WHOLE result (the fail-closed floor is preserved, but a
       single bad weld loses the model, not just the patch).
     - SortGeometry/SetNormals/SetEpsilon are pure reorder/recompute, no re-weld -
       they genuinely compose.

4. PER-CARRIER SCOPE (correct the "localized" claim). Two locality proxies per
   carrier (dirty core = largest near-degenerate face cluster):

     carrier    fatal reason        core faces   core span (frac of diagonal)
     -------    ----------------    ----------   ----------------------------
     Havoc      NonManifoldEmit     ~8.5%        ~0.96   (model-spanning)
     GT7863     NonManifoldEmit     <1%          ~0.12
     Offset1    NonManifoldEmit     <1%          ~0.03
     openscad   NonManifoldEmit     ~2.6%        ~0.41   (model-spanning)
     GT7081     ArrangementBudget   <1%          ~0.56   (>90% sub-eps crit gaps)

   READING: COST (face count) IS bounded corpus-wide - every dirty core is under
   9% of faces, so C's exact-arrangement region does not swallow the model by face
   count. The "scales with degeneracy" cost story SURVIVES on face count. BUT the
   core is SPATIALLY MODEL-SPANNING on the single-feature carriers (Havoc ~0.96 of
   the bbox diagonal, GT7081 ~0.56, openscad ~0.41). The rev-1 "near-coincidences
   are localized, a clean boundary exists" conflated "one contact region" with
   "spatially small"; they differ. The near-degenerate contact is a model-spanning
   SHEET (few tessellated faces, huge reach), so the clean-boundary weld seam
   threads across the whole model - R2 weld exposure scales with MODEL size, not
   degeneracy. And GT7081 is delocalized at the ARRANGEMENT level (the sub-eps
   critical run), so C grows to most of the arrangement = B; today it fails
   ArrangementBudget, so C would run full B to upgrade it, at B's cost.


## Zero silent wrong-resolves (recorded finding)

Round-1 (c lane, ATTACK 4) spot-checked C's foundational premise - that the eps
pipeline never resolves WRONG silently, it only fails closed. For three RESOLVING
corpus cases (Offset2, Offset3, Cray), exact rational w_S(input soup) was compared
against the eps pipeline's output containment at sampled interior/exterior points.
Everywhere measured, they AGREE - zero disagreements. The zero-oracle-wrong record
holds on this spot check; refusal (fail-closed) is the only observed failure mode.
This supports C's premise where tested and means C only ever UPGRADES a fail-closed
region to correct, never corrects a silent error (there were none to find). SCOPE:
all three are union / self-overlap (all-positive multiplicity); the subtraction /
negative-multiplicity regime (R5) is untested here.


## Surviving architecture (C): staged sketch (rev 2)

Rebuilt with the v4-plan discipline: a stage-0 go/no-go measurable WITHOUT the
exact kernel (the second cheap decision point - the probe was the first), an
objectively-checkable gate and stub-proof pin per stage, and an explicit
stage -> v3-machinery kept/deleted map so "MAX v3 reuse" is a checked claim, not
an adjective.

STAGE 0 - GO/NO-GO WITHOUT THE KERNEL (build this before Stage 1). C's own load-
  bearing assumption is R1 (clean-boundary growth is bounded on the corpus), and
  R1 is measurable with NO exact arithmetic: run the eps pipeline, build the
  localizer (precondition 1), grow each flagged region along guard-clean input
  edges, and measure the size distribution. Also produce the selective-weld design
  note (precondition 2). This is the deciding experiment the plan must not defer.
  GATE: bounded clean-boundary growth -> build C (Stage 1). Unbounded
  (percolation) -> only the whole-mesh B scope survives, re-decide before spending
  on the kernel.
  INITIAL ANSWER (c lane, folded): face-count bounded (< 9% everywhere) = GO on
  cost; but spatially MODEL-SPANNING on Havoc/openscad and DELOCALIZED on GT7081
  (C = B there). So the gate is a PARTIAL go: build C for the face-count-bounded,
  spatially-compact carriers; dial GT7081 to the whole-mesh B scope explicitly.
  The localizer and the selective weld are the two pieces of design work this stage
  must finish before any kernel spend.
  v3 map: reuses the FatalReason returns at BuildSlabs/EmitCaps/EmitStrips as the
  failure SEED, plus the NEW provenance channel; the eps pipeline is untouched.

STAGE 1 - EXACT PREDICATE KERNEL (the reusable primitive; cost scales with
  degeneracy, not size).
  - Adaptive orient3d (FP filter first; exact expansion only when the FP
    determinant is not sign-certain). Cheap in the common non-degenerate case.
  - Exact w_S(p): segment-crossing count on INPUT vertices (orient3d signs), one
    fixed global symbolic-perturbation convention (the Shadows pattern lifted to
    3D orient) so a ray grazing a shared edge/vertex resolves identically for
    every probe. This is the load-bearing consistency guarantee of the reopen.
  PIN (stub-proof): differential sign-agreement vs a throwaway rational driver on
  the probe's query points, bit-for-bit. The probe lanes already pass a partial
  instance (two independent algorithms agree bit-for-bit).
  v3 map: NEW surface area alongside the bulk path; DELETES nothing (the tree has
  no exact kernel to reuse).

STAGE 2 - FLAG + BOUND. The eps pipeline runs unchanged. On a guard fire, the
  localizer (Stage 0) supplies the seed and the incident INPUT faces; GROW the
  region until its boundary is guard-clean: every boundary edge is an INPUT edge
  the eps pipeline resolved without a guard. The region is a submesh whose
  interface with the exterior is shared input geometry (bitwise equal on both
  sides).
  GATE: the grown boundary is entirely guard-clean input edges (checkable).
  v3 map: reuses the localizer's seed + the FatalReason returns; adds the grower.

STAGE 3 - LOCAL EXACT ARRANGEMENT. Compute the exact arrangement of the submesh:
  rational intersection points among its faces, cell decomposition, radial-order
  DCEL around each internal edge (tie-free predicates - the reopen), classify each
  half-face by exact w_S (Stage 1), extract the oriented {w_S >= 1} boundary. The
  topology is exact; constructed intersection coordinates round to double only at
  emission, behind a LOCAL validity gate (2-manifold, no self-intersection within
  the patch).
  GATE: local 2-manifold + no self-intersection (checkable).
  v3 map: NEW.
  COST (bvc lane, folded): level-1 predicates on INPUT are a few machine words
  (a real Havoc coord is an exact m/2^k with m ~53-bit; orient3d exact fallback
  ~165-bit). Building the arrangement reaches level-2 predicates on CONSTRUCTED
  points (order intersection points along an edge; point-in-face) whose operands
  are ~165-bit rationals, so a degree-3 op yields ~500-bit numerators; feeding
  constructed coordinates into further predicates cascades (~1500-bit at level-3).
  Tractable (GMP handles 500-bit) but NAME it: keep decisions as INPUT predicates
  wherever possible so bit-growth stays bounded.

STAGE 4 - WELD. Stitch the local patch to the eps exterior along the shared
  input-edge boundary using the SELECTIVE weld (precondition 2): patch boundary
  adopts the exterior's vert positions; interior verts are weld-exempt (a bounded
  local weld, NOT the global eps grid). Handle the downstream-care entries
  (precondition 3): SplitTouchingSheets must not reject legitimate tangent patches.
  GATE: boundary vertices are INPUT vertices, shared BITWISE (checkable); no
  patch-interior vert merged.
  v3 map: reuses the assembly weld WITH a bounded local exemption, and
  SplitTouchingSheets WITH a needs-care patch for genuine tangents.

STAGE 5 - VALIDITY GATE + FALLBACK. If the welded result has unpaired halfedges or
  fails validity, fall back to the B scope on that connected component, or fail
  closed. The honest fail-closed contract is the floor; the exact resolver only
  ever UPGRADES a fail-closed region to correct.
  v3 map: reuses the global IsManifold/Is2Manifold gate as the floor.

Second-probe status: the local mechanism's CORE is validated at toy scale on the
real cluster (exact winding stratifies the sheet region cleanly; the lump
survives; the boundary is manifold-by-construction). The UNVALIDATED load-bearing
assumption is Stage 0's clean-boundary growth at corpus scale (partially answered:
face-count GO, spatial-spanning caveat; see Risks).


## Risks (the attack surface for the adversarial round)

R1. CLEAN-BOUNDARY GROWTH (the load-bearing assumption). C's cost and correctness
    both hinge on growing every flagged region to a guard-clean input-edge
    boundary. Round-1 (c lane) gave the initial measurement, now Stage 0's job:
    face-count bounded corpus-wide (< 9%), but SPATIALLY model-spanning on
    Havoc/openscad and DELOCALIZED on GT7081 (C = B). So the risk is REAL and
    partly realized: on single-feature carriers the boundary threads across the
    whole model, and on GT7081 no bounded neighborhood is clean. Open work: finish
    Stage 0's growth measurement with the localizer built.

R2. WELD RELOCATION. A boundary input edge may be crossed by a third face near-
    degenerately, so the eps exterior subdivides it at an eps-placed vertex the
    exact patch places exactly -> an unpaired halfedge at the weld. Round-1
    sharpened this: because the dirty core is spatially model-spanning, the weld
    seam is long, so R2 exposure scales with MODEL size, not degeneracy. The
    selective weld (precondition 2) is the mechanism-level mitigation (boundary
    conformance to exterior positions); the residual attack is a boundary edge with
    a hidden near-degenerate crossing.

R3. OUTPUT SNAP-ROUNDING. Rounding constructed intersection coordinates to double
    can create NEW near-coincidences at output (two exact-distinct patch vertices
    round together), silently breaking 2-manifoldness below the validity gate's
    resolution. Iterated snap-rounding may be required. ADDED (bvc lane): to bound
    bit-growth (Stage 3), the cascade may need INTERMEDIATE snapping, and
    intermediate snapping can change TOPOLOGY, not just output coordinates - a
    stronger failure than the output-only case. Attack: a patch whose exact
    vertices are closer than eps after the arrangement, so no snap is both valid
    and faithful.

R4. ONE CONVENTION ACROSS ALL PREDICATES (the subtle reopen residue). Exactness
    cures the FP-disagreement kill ONLY IF a SINGLE symbolic-perturbation
    convention governs EVERY predicate family (orient3d, the winding segment-
    crossing, the arrangement's intersection ordering). A convention consistent
    within orient3d but different for the winding ray reintroduces cross-boundary
    disagreement - the kill in disguise. Round-1: the runtime LOAD is light
    (ray degeneracy is boundary-only, ~0% on interior cell points), so this is a
    correctness constraint, not a hotspot. Attack: construct a mixed-convention gap.

R5. SUBTRACTION / NEGATIVE-MULT REGIME. Every measurement is union-composed (all
    multiplicities +1), including the zero-silent-wrong-resolve spot check.
    Subtraction inputs carry negative multiplicities; w_S can go negative; the
    {w_S >= 1} rule still holds but the strata and the winding jumps differ.
    Unmeasured. Attack: a subtraction carrier with an inward shell.

R6. COINCIDENT / COPLANAR FACES (measure-zero overlaps). Exact w_S is undefined ON
    a coincident face (measure zero); the fill needs a tie-break there. Round-1
    made this concrete, not hypothetical: GT7863 has a vertex exactly on a face,
    Offset1 and openscad have exactly-coincident verts, and the intrinsic carriers
    carry these AT THE INPUT. The tie-break must match the library's (the
    coincident-wall cleaning the witness lane's Case B relied on). Attack: coplanar
    coincident sheets where w_S jumps by 2 across a doubled wall, or an operand wall
    buried inside another.

R7. KERNEL COST + CERTIFICATION. Adaptive orient3d must be certified correct
    (Shewchuk-class); the exact arrangement's rational constructions reach ~500-bit
    at level-2 even when bounded (Stage 3 cost). Attack: a flagged region that is
    bounded by face count but internally DENSE - GT7081's sub-eps critical run
    (~34k near-coincident planes in a thin band) is the concrete instance, where
    the local arrangement is combinatorially large despite a small face count.


## What is settled (do not re-litigate)

The probe results (now confirmed bit-for-bit by two independent algorithms) and
every measured kill across the campaign; A0; the witness theorem; the escape
kill's lump invariance. Variant E is KILLED (placement is not classification). B
and C are MERGED into one mechanism with a scope dial. The zero-silent-wrong-
resolve spot check holds on the union corpus. The v3/v4 branches stay as they are;
zero-oracle-wrong; correct-or-honest. The only object reopened is the FP-premised
RSI-#3 option-3 kill, which does not survive exact winding with one global
convention - and which, corrected, is the exact arrangement this doc designs C to
pay for only where the eps pipeline fails closed.
