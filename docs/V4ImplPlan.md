# V4 implementation plan: two-pass global identity emission

STATUS: STAGED IMPLEMENTATION PLAN (gates-first contract), branch
explore/sweep-plane-3d-v4. The converged design is docs/MaintainedEmission3D.md
variant (iii) - two-pass global identity. This document turns that design into a
pinnable, attackable build. It is NOT the pipeline; it is the plan the adversarial
round-1 lanes attack next. v3 (docs/SweepEmit3D.md, src/overlap3*) stays intact and
green until a stage DELIBERATELY and PINNEDLY changes it.

GATES-FIRST IS LAW (v3's thrice-confirmed lesson). Build agents rescope fail-closed
machinery and weaken gates while narrating happy-path green; that green is
systematically vacuous. The countermeasures that worked, and that every stage below
inherits: contractual gates enumerated BEFORE code; MECHANISM PINS authored red-first
that CANNOT pass against stubs (white-box - they assert internal invariants, not just
outcomes); narrow single-mechanism stages; independent re-runs of every green claim.

SETTLED, not re-litigated here: the design verdict and the six-candidate kill table
(MaintainedEmission3D.md); the two structural conditions; the absolute
zero-oracle-wrong posture; the fence rule (ALL existing tests are every stage's
fences); no counts in docs; ASCII; house-voice commits.

## 0. What the design settled, and what it did NOT

Variant (iii) DISSOLVES the per-plane-independence coupling (v3's many independent
per-plane cap self-locations) and does NOT relocate, under two NAMED structural
conditions that are the design requirements of the build:

  (C1) the near-degeneracy collapse is re-derived from INPUT coordinates via an
       ORDER-INDEPENDENT event ordering (exact-lex event set + exact predicates),
       never consuming a constructed coordinate in a discovery-order-dependent way;
  (C2) each near-degenerate cluster is CANONICALIZED exactly once.

It did NOT settle that the identity-carrying emission reaches ORACLE-TRUE on the
hardest carriers. The escape kill (probe 3) proved a bounded POSITIONS-ONLY
composition closes Havocglass8 to a VALID genus-0 winding-correct 2-manifold that is
oracle-WRONG by a macro volume (~14000x the eps bound), via a macro-region flip
inseparable from closure. Variant (iii)'s claim is that the UPSTREAM arrangement
identity emission discards - which faces and seams actually meet at a junction, not
where its images landed - is exactly what disambiguates that flip. That claim is the
RESEARCH FRONTIER the pipeline exists to reach; it is UNPROVEN on Havoc. Therefore:

  the project's DONE criterion is NOT "the corpus resolves". It is "the identity
  table is built to C1/C2, emission consumes it, the isolated case-1 twin closes
  oracle-true, and the corpus is re-adjudicated with ZERO oracle-wrong resolves".

Carriers that the identity-carrying emission cannot close oracle-true stay RECORDED
FAIL-CLOSED contracts, never silent oracle-wrong resolves. The gates protect against
regression and oracle-wrongness regardless of how far the research frontier moves.

## 1. Architecture delta

### 1.1 The two passes

PASS 1 - GLOBAL IDENTITY COLLECT. Runs once, after FindSeams, before BuildSlabs.
CONSUMES the arrangement's junction-defining geometry: the canonical merged verts
(input-quantized, order-invariant), the seams and their endpoints, the criticals
(seam-seam crossing x's, degenerate-contact x's), and the crossing bundles (the
density set). PRODUCES the IDENTITY TABLE: the named 3D junctions, each junction's
cluster COLLAPSED ONCE from input coordinates, and the density-side thinned critical
set. It is a BATCH collect (not incremental - that is the retired variant i); the
emission pass reads it.

PASS 2 - EMISSION. v3's stages, now reading the identity table where they self-located
before. Every cap and every strip that images a junction takes that junction's
canonical position from the table, so no per-plane self-location exists and the
assembly weld never reconciles a junction twin (it cannot arise).

### 1.2 New data structures and their invariants

```
struct Junction {
  vec3 pos;                  // canonical position, coordinate-determined
  std::vector<int> members;  // the arr.vert cluster this junction collapses
  bool spread;               // classifier: true = one junction, false = distinct
};
struct IdentityTable {
  std::vector<Junction> junctions;
  std::vector<int> vert2junction;   // arr.vert id -> junction id
  std::vector<double> thinnedCriticalXs;  // density-side, reps >= eps apart
};
```

INVARIANTS, each a pin (section 3):
- `Junction::pos` is a PURE FUNCTION of the member INPUT coordinates (exact
  lexicographic-min member, or an exact-predicate representative), hence ORDER-FREE.
  It is NOT a merge first-match representative - R1 measured that order-dependent at
  the eps level (36/40 nonidentical under vertex permutation). This is C1's teeth.
- The whole `IdentityTable` is BITWISE-IDENTICAL under input permutation (C1).
- Each junction's `pos` is chosen EXACTLY ONCE during construction (C2) - a
  per-junction counter, not a per-consumer/per-member re-derivation.
- The classifier `spread` FIRES on every case-1 twin and EXCLUDES the case-2 distinct
  pair (predicate: maxDx <= dxThresh*eps AND maxTr >= maxDx; the transverse-dominated
  test is load-bearing - eps-ball alone merges the distinct pair).

### 1.3 Adjudications the design demanded

PROVENANCE CHAINS are NOT subsumed; they COMPOSE, orthogonally. The provenance
channel (classId / edgeSubdiv, landed) threads per-input-EDGE SUBDIVISION identity -
the ordered arrangement vertices ALONG one edge (a 1D object). The identity table
threads per-JUNCTION identity - which emitted images are the same 3D point (a 0D
object). Different axes. Pass-2 emission reads BOTH: junction positions from the
identity table (replacing the SlabResolver's self-location), edge subdivisions from
the provenance chains (unchanged). The killed "provenance-label twin merge" was
positions-only end-to-end and named NO canonical for the weld images; the identity
table supplies the canonical from the GLOBAL collect, which names every junction from
INPUT geometry (including what v3 produces as per-slab weld points). That is the
distinction that makes the table not the killed candidate. VERDICT: provenance chains
SURVIVE as the subdivision foundation; the identity table is a NEW orthogonal channel
that SUBSUMES only the SlabResolver's positions-only self-location role.

DENSITY-SIDE BUNDLE HANDLING lives in PASS 1. Crossing-bundle thinning (VALIDATED-SAFE
in the wall-A arc: it touches only crossing x's, which SEAMS consumes for x alone, and
keeps representatives >= eps apart so no macro feature's bounding vertex moves) is the
density-side of the collect. It produces `thinnedCriticalXs`. It is landed only under
the zero-oracle-wrong posture (section 3, pin 5): thinning that flips a single mesh to
an un-oracle-checkable resolve is FORBIDDEN.

ARRANGEMENTBUDGET STAYS. It is the honest refusal for genuinely-degenerate-VERTEX
density: GT7081's endpoint clusters stay dense after crossing thinning (skeleton:
maxCluster=77 near-degenerate verts survive), so it still refuses. Its firing DOMAIN
SHRINKS: cases whose density is pure crossing OVER-INCLUSION (self_intersect A/B, hull)
section after thinning and land on their emission-stage outcome instead. The budget is
not recalibrated; it is the resource backstop that converts a de-facto hang into a
recorded refusal.

### 1.4 v3 stage fates in pass 2

| v3 stage            | fate under variant (iii)                                        |
|---------------------|-----------------------------------------------------------------|
| Canonicalize        | SURVIVES; feeds pass 1.                                          |
| FindSeams           | SURVIVES; its verts/seams/criticals are pass-1 INPUT; its        |
|                     | criticalXs are thinned in pass 1.                               |
| BuildSlabs          | SURVIVES as the per-critical collect KERNEL; its section verts    |
|                     | BIND into the global identity (a junction-imaging vert takes     |
|                     | the junction's canonical, not a per-slab position).             |
| SlabResolver        | Self-location (Extend via InterpolateSafe) DELETES for           |
|                     | junction-imaging verts - the wall-A locus. Degenerates to a      |
|                     | table read; non-junction endpoints may keep track subdivision.  |
| ComputeCap          | CONSUMES identity (stages B, C).                                 |
| EmitStrips          | CONSUMES identity; provenance chains unchanged (stage C).       |
| BuildImpl / split   | Assembly weld is a no-op for junction images (they share the     |
|                     | canonical), so the twin cannot arise. Is2Manifold gate KEPT.    |

DIES: per-plane independence of junction self-location; track extension AS the
junction LOCATOR; assembly-weld reconciliation of junction twins. BORN: the identity
table (pass 1); the junction-identity read in caps/strips/assembly.

## 2. Stage sequence

Each stage is narrow (single mechanism), additive-first (new channels default-off /
unconsumed until the consuming stage lands - the negEdges/SweepCapture precedent),
independently green, and independently re-run. For each: mechanisms, files, v3 touch,
CONTRACTUAL GATES (its fences = ALL existing tests, plus its new white-box pins), its
red-first anchor, and its DONE criterion.

The A/B/C/D shape follows the design's decomposition, AMENDED at the B/C boundary:
B and C split along the transverse-vs-x-residual axis, which is the GT7863
bound-and-diverged anatomy (the transverse yz divergence collapses at B; the cap-plane
x-residual closes at C). That split is empirical, not fictitious (section 4, R-BC).

### Stage A: identity table construction (pass 1)

MECHANISM: a new `BuildIdentityTable(arr, eps)` between FindSeams and BuildSlabs.
Cluster arr.verts by an eps-ball union-find (from-input positions, exact predicates);
per cluster apply the classifier; choose the canonical ONCE (coordinate-determined).
SCOPE: vertex junctions (arr.verts collapsed). Resolving seam-seam crossings to 3D
triple-point junctions is a FLAGGED heavier sub-mechanism (section 4, R-TRIPLE) NOT in
stage A - the skeleton shows the vertex-junction table is bounded and buildable, and
the case-1 twin's underlying junction V is itself an arr.vert.
FILES: overlap3.cpp (new function), overlap3.h (IdentityTable, ArrangementGeometry
gains an `IdentityTable identity` field default-empty), test hooks expose it.
V3 TOUCH: ADDITIVE. ArrangementGeometry gains a field; emission IGNORES it. Default
behavior BITWISE unchanged.
GATES: ALL existing tests (fence) + Pin_IdentityTable_OrderInvariant,
Pin_IdentityTable_OnceOnly, Pin_Classifier_Fence (section 3, pins 1-3).
RED-FIRST ANCHOR: a constructed fixture with a known near-degenerate cluster (or
Havoc's) whose junction the table MUST name; the pin reds against an empty /
per-vertex-identity stub.
DONE: table built + three pins green (mutation-verified) + suite green + emission
bitwise-unchanged (a diff of emitted geometry against HEAD is empty).

### Stage B: caps consume junction identity (transverse)

MECHANISM: ComputeCap / BuildCapEdgeSet reads a junction-imaging section vert's
canonical (y,z) from the identity table instead of the SlabResolver's self-located
extension. Both caps of one junction then place it at the SAME (y,z), dissolving the
transverse divergence. A lookup with fallback to the resolver for non-junction verts
(additive-shaped).
FILES: overlap3.cpp SlabResolver::Extend (identity-first path), BuildCapEdgeSet.
V3 TOUCH: the SlabResolver gains an identity-first branch; ComputeCap reads the table.
GATES: ALL existing + Pin_CapImage_CanonicalTransverse (section 3, pin 4a): a
junction's two cap images have BIT-IDENTICAL (y,z), mutation-verified against
self-location; plus no fence regression.
RED-FIRST ANCHOR: the isolated case-1 twin fixture - the transverse bind reds (yz
diverges) without the table, greens (yz collapses bitwise) with it, reproducing
GT7863's R2 bound trace.
HONESTY / the split from C: stage B alone does NOT close the twin. The cap-plane
x-residual survives (GT7863 bound-and-still-diverged: yz bitwise, ~1.1-eps x-gap). B's
DONE is the transverse-identity pin + no regression, NOT twin closure.
DONE: cap images at canonical yz pinned (mutation-verified) + suite green +
emitted-geometry diff confined to the transverse bind (no new resolves yet).

### Stage C: strips + assembly join by name (coordinated re-emission)

MECHANISM: the junction reaches the EMITTED vertex. A junction is emitted at ONE 3D
point - including its canonical x, not each cap plane's own x - so the two adjacent
cap planes imaging one junction UNIFY into one plane; strips join by junction NAME
(provenance chains supply the subdivision between junctions); the assembly weld is a
no-op for junction images. This is the coordinated cap+strip re-emission around a
near-degenerate cluster the wall-A memo named as the real fix.
FILES: overlap3.cpp EmitCaps/ComputeCap (cap-plane unification for junction-imaging
pairs), EmitStrips, BuildImpl.
V3 TOUCH: pair-canonical cap emission; strip chains; assembly weld.
GATES: ALL existing + Pin_IsolatedTwin_ClosesOracleTrue (section 3, pin 4b) +
Pin_ZeroOracleWrong meta-gate (pin 5) + the per-carrier corpus contracts (pin 4),
tightened so any NEWLY-resolving carrier is oracle-checked.
RED-FIRST ANCHOR: the isolated case-1 twin fixture reds (NonManifoldEmission) without
C, resolves ORACLE-TRUE (volume + genus + winding) with C - the escape's proven
bounded partial win, now landed as the mechanism.
HONESTY / the research frontier: on real ENTANGLED carriers (Havoc) the escape kill
proved positions-only closure is oracle-WRONG. Stage C carries the upstream junction
identity (incidence), which is what disambiguates spread from distinct BEYOND
positions - but whether that reaches oracle-true on Havoc is UNPROVEN. C's gate is
zero-oracle-wrong: a carrier C cannot close oracle-true STAYS a recorded fail-closed.
DONE: isolated twin oracle-true + corpus re-adjudicated with zero oracle-wrong (some
may newly resolve oracle-checked; others stay recorded contracts) + suite green.

### Stage D: guard dissolution + corpus re-adjudication + density + perf

MECHANISM: delete the now-structurally-dead machinery and re-adjudicate the guard
family; land crossing-bundle thinning (density-side of pass 1) gated by the oracle
posture; re-adjudicate ArrangementBudget's domain; perf.
PER GUARD (delete at the stage that kills it, or keep the fail-closed posture):
- SlabResolver junction self-location: DEAD at C (junctions read the table). Delete,
  pin-dead (a mutation restoring self-location reds the twin pin).
- Chain-plane wide-run rule + SubEpsFeature wide-run cap backstop: RE-ADJUDICATE.
  Under identity emission the run-spanning interpolation is replaced by junction
  identity, so the wide-run backstop MAY become dead. KEEP until pinned dead; do not
  delete on narration.
- SubEpsFeature single-face coverage guard (BuildSlabs): KEEP. A macro face wholly in
  a merged ulp-run is a genuine merged-run hazard the identity table does not touch.
- ArrangementBudget: KEEP; domain shrinks after thinning (section 5).
- Is2Manifold gate (BuildImpl): KEEP - the final honesty check.
- EngineIdConflict / CoplanarOverlap / EdgeInPlane: already retired (v3).
FILES: overlap3.cpp (guard deletions + thinning in pass 1), overlap3_sweep.cpp
(budget domain), tests (guard-dead pins + corpus re-adjudication + tol-invariance).
GATES: ALL existing + guard-dead pins + the corpus contracts + perf checkpoints
(section 5) + Pin_ZeroOracleWrong.
RED-FIRST ANCHOR: each guard-dead pin reds if the guard's mechanism is restored;
thinning's tol-invariance pin reds if a single-mesh resolve moves under tolerance.
DONE: dead guards deleted and pinned-dead; thinning landed+oracle-fenced OR deferred
with its tol-invariance evidence; ArrangementBudget domain re-adjudicated; perf within
bound; corpus zero-oracle-wrong.

## 3. The pin list (authored as specifications now)

Each pin is white-box (asserts an internal invariant), red-first, and carries a
MUTATION-VERIFICATION recipe: a named stub the pin MUST red against, verified by the
author independently re-running the red (not trusting a build agent's green).

PIN 1 - Pin_IdentityTable_OrderInvariant (C1). Permute the input (vertex + triangle
reorder, several seeds), rebuild the table, assert the sorted junction set is
BITWISE-identical (exact double compare on canonical positions) and the collapse
partition matches. WHY A STUB CANNOT PASS: a merge first-match representative is
order-dependent at the eps level (R1: 36/40 nonidentical under vertex permutation);
only a coordinate-determined choice (exact lex-min / exact-predicate rep) passes.
MUTATION: replace the canonical choice with `members[0]` (discovery-order first) ->
reds under permutation. SKELETON EVIDENCE: bitwise on control + Havoc (section 4).

PIN 2 - Pin_IdentityTable_OnceOnly (C2, counter). Instrument construction with a
per-junction "canonical assigned" counter; assert it equals the junction count (each
canonical computed exactly once), and that a per-consumer path does NOT re-derive it.
WHY A STUB CANNOT PASS: a stub that canonicalizes inside the cap/strip loop
re-derives the representative per consumer -> counter > junction count. MUTATION: move
the canonical choice into the emission loop -> counter reds.

PIN 3 - Pin_Classifier_Fence (both directions). On a fixture carrying BOTH the case-1
spread twin AND the case-2 DISTINCT pair (order tens of eps apart in x, sub-eps in
yz), assert the classifier FIRES on the twin (collapses to one) and EXCLUDES the
distinct pair (never clustered nor collapsed). WHY A STUB CANNOT PASS: an eps-ball
collapse without the transverse-dominated test, or a dxThresh wide enough to reach the
distinct pair's x-gap, MERGES the distinct pair (R4's kill). MUTATION: drop the
`maxTr >= maxDx` conjunct OR widen dxThresh past the distinct pair's dx -> the pair
merges, pin reds. This is the R4 verdict-flipper's fence made a pin.

PIN 4 - per-carrier targets.
  (4a) Pin_CapImage_CanonicalTransverse (stage B): a junction's two cap-input images
       have bit-identical (y,z). MUTATION: restore SlabResolver self-location -> the
       two images diverge, pin reds.
  (4b) Pin_IsolatedTwin_ClosesOracleTrue (stage C): a hand-built isolated case-1 twin
       resolves to a clean 2-manifold, volume-neutral, oracle-true (a+b Boolean).
       MUTATION: disable cap-plane unification -> NonManifoldEmission, pin reds.
  (4c) PAIRS get EXACT Boolean oracles (volume + genus + winding): Havoc, GT7863,
       Cray, GT7081, hull. A resolve MUST match; a fail-closed MUST be a named guard.
  (4d) SINGLES get TOL-INVARIANCE contracts (resolve volume identical across a
       tolerance sweep, verified by EXACT bit comparison, not a fixed-precision
       print) + 2-manifold validity: Offset1-4, openscad, self_intersect A/B. A single
       mesh is NEVER asserted resolved-correct without an oracle.

PIN 5 - Pin_ZeroOracleWrong (meta-gate, absolute). Across the WHOLE corpus + synthetic
suite, in every mode, ZERO oracle-off resolves. Any resolve is oracle-checked (pairs)
or tol-invariant-validity-checked (singles); any carrier that cannot be checked stays
FAIL-CLOSED. WHY IT EXISTS: the escape kill closed Havoc to a VALID manifold that was
oracle-WRONG by a macro volume. Validity is NOT correctness. This meta-gate forbids
exactly that outcome and is the gate the whole project is measured against.

MUTATION-VERIFICATION DISCIPLINE (every pin): the author re-runs the red against the
named stub independently. A green claim from a build agent is re-run before it is
believed. This is v3's thrice-confirmed countermeasure, not optional ceremony.

## 4. Risks (the round-1 attack surface)

Named for a skeptic to break, with the skeleton evidence where it exists.

R-DENSITY (the identity table's size on GT7081's million-critical bundles). EVIDENCE
(skeleton, section below): the vertex-junction table is BOUNDED at the vertex count
(GT7081: ~18.5k junctions from ~18.6k arr.verts), NOT the crossing-bundle count
(~1.37M criticals); table build is sub-second (0.02s clustering). So the table's SIZE
is not the density risk; the ~1.37M criticals feeding BuildSlabs are (v3's existing
ArrangementBudget/thinning problem, ORTHOGONAL to the table). ATTACK: this holds ONLY
for the vertex-junction scoping (R-TRIPLE).

R-TRIPLE (does the vertex-junction scope of stage A suffice, or must pass 1 resolve
triple points?). The skeleton table is built on v3 arr.verts, which does NOT resolve
seam-seam crossings to 3D junctions. A FULL arrangement completion (the RSI-#3 object
the design names) that names triple-point junctions could add up to ~1.37M junctions
on GT7081 - the table then is NOT bounded. ATTACK: construct a carrier where the
spread-vs-distinct disambiguation needs triple-point incidence ABSENT from the
vertex-junction table (the dense endpoint cluster, skeleton maxCluster=77, is where to
look). The plan's "vertex junctions suffice for the case-1 twin" is an UNPROVED
property; the twin's V is an arr.vert, but the disambiguation on dense clusters may not
be.

R-ONCE-VS-STRICTFP (does canonicalize-once conflict with v3's strict-FP per-slab
sectioning?). v3 builds every slab whose xMid is strictly interior and runs an
INDEPENDENT 2D arrangement per slab. If a junction is canonicalized once globally but
each slab still sections independently, the per-slab section and the global canonical
are two structures that can DIVERGE - the wall relocated into the bind. ATTACK: the
collect KERNEL (per-slab 2D arrangement) must BIND its section verts into the global
identity; show a case where a slab's independent strict-FP section disagrees with the
junction's global canonical and the bind does not reconcile it.

R-BC (is stage B green without C, or is the split fictitious?). ARGUMENT: the split is
the GT7863 bound-and-diverged anatomy - the transverse yz divergence collapses at B
(pinnable: bit-identical cap-image yz), the cap-plane x-residual survives to C. The B
pin (4a) is a real internal invariant a C-less implementation genuinely satisfies and
that is not vacuous. ATTACK: show the B pin is either vacuous (passes against a stub)
or requires C (the split is fictitious).

R-CLOSURE (the plan assumes a closure it has not earned). The escape kill proved
positions-only closes Havoc oracle-WRONG; variant (iii)'s oracle-true claim rests on
upstream incidence identity that is UNPROVEN on Havoc. GUARD: the plan does NOT gate on
"Havoc resolves"; stage C's DONE is isolated-twin-closes + zero-oracle-wrong, and the
research frontier is named as a frontier (section 0). ATTACK: find where a stage's DONE
criterion QUIETLY assumes carrier closure instead of zero-oracle-wrong.

R-THINNING-ORACLE (crossing-bundle thinning flips single meshes with no oracle).
Thinning turns self_intersect A/B from ArrangementBudget refusal to a resolve that has
NO a+b oracle. The plan defers landing thinning until an oracle exists and substitutes
a tol-invariance contract. ATTACK: is tol-invariance (stability under tolerance) a
sufficient stand-in for correctness? It proves the resolve is not a thinning artifact;
it does NOT prove it is geometrically right. The plan treats it as validity-plus, not
oracle - attack whether that is honest.

R-ORDER-UNMEASURED-AT-DENSITY (C1 bitwise-ness is measured only on small fixtures).
The skeleton confirmed bitwise order-invariance on the control and Havoc, but GT7081
order-invariance is unmeasured (too heavy: quadratic compare x permutations x
multi-second hooks). Denser input may make seam-vert POSITIONS order-dependent, which
even a lex-min choice cannot fix - only the exact-predicate re-derivation from input
(C1 proper) can. ATTACK: the plan leans on C1 holding at density; the evidence is
small-scale.

## 5. Perf checkpoints

BOUND: the full suite within ~2x of the current ~107s, and no single test > ~120s,
checked PER STAGE (a stage that regresses perf is not done).

- Stage A: table build is O(nVerts) eps-ball clustering - skeleton measured 0.02s on
  ~18.6k verts, negligible. The order-invariance PIN runs several permutations x hooks;
  keep it on SMALL fixtures only (control + Havoc) and do NOT permute GT7081 (its hooks
  are ~30s each). CHECKPOINT: stage A adds a low single-digit percentage to the suite.
- Stages B/C: the identity read is an O(1) map/hash lookup per section vert; it should
  not regress. CHECKPOINT: no test > ~120s; suite within 2x.
- Stage D / the GT7081-hull-budget interplay once bundles collapse: crossing-bundle
  thinning REDUCES GT7081's ~1.37M criticals by roughly two orders of magnitude
  (wall-A evidence: self_intersect A/B flip from budget refusal to resolve after
  thinning), which SPEEDS BuildSlabs and lets the crossing-over-inclusion cases (hull,
  self_intersect) section. GT7081's ENDPOINT clusters stay dense (maxCluster ~ tens),
  so ArrangementBudget still fires there - the budget remains the resource backstop and
  is NOT recalibrated up (relaxing it runs the dense case to the emission blowup and
  STILL fails, strictly worse). CHECKPOINT: after thinning, hull sections and lands on
  its cap-plane outcome (resolve if C closed it, else recorded fail-closed);
  self_intersect resolves tol-invariant (oracle-fenced or deferred); GT7081 still
  refuses ArrangementBudget within seconds and bounded memory; the suite stays within
  bound with the budget as backstop.

## 6. Walking-skeleton probe evidence

A thin real version of the load-bearing joint - the identity table - built on v3's
existing arrangement (arr.verts from RemoveOverlaps3D_TestHooks) and run on one small
fixture (Havocglass8) and one clean control (Gate4a kWedges8), plus a density
measurement on GT7081. Env-gated (OV3_V4SKELETON[_GT7081]), default-skip, REVERTED
after capture (machinery delta zero); the recipe and full numbers live in the lab
notebook (.claude/lane-reports/v4impl-stage0-*). The table clusters arr.verts by an
8*eps ball, applies the probe-3 classifier (maxDx <= 0.7*eps AND maxTr >= maxDx), and
chooses the canonical as the lexicographic-min member (coordinate-determined).

C1 (order-invariance), permute input x8 and compare junction sets bitwise:
- Control (kWedges8): BITWISE identical (worst nearest-deviation 0 eps, no count
  difference, no exact miss).
- Havocglass8: BITWISE identical (worst deviation 0 eps). One collapse fired (a 2-vert
  cluster, dx ~ 0.45 eps - sub-eps-x, classifier clean, zero fence violations).
So a coordinate-determined canonical yields a bitwise-identical table under
permutation at skeleton scale, on both a clean control and the mapped carrier - direct
evidence for C1 (with R-ORDER-UNMEASURED-AT-DENSITY the honest limit).

C2 (once-only), Havoc: the one collapse computed ONE representative (vs ~196 a naive
per-member canonicalization would compute) - the collapse-once shape, concretely.

DENSITY, GT7081 (under ulimit -v 4GB): ~18.6k arr.verts and ~1.37M criticals produce a
~18.5k-junction table (spread ~17.6k, distinct ~0.3k, maxCluster 77); table build
0.02s; hooks 30s (BuildSlabs hits ArrangementBudget, arr still populated). The junction
count tracks the VERTEX count, not the crossing count - the R-DENSITY evidence, and the
R-TRIPLE caveat's setup (this is the vertex-junction table; triple-point completion is
the unmeasured heavier object).

NOTE ON LEVELS: the skeleton table lives UPSTREAM (arr.verts). The escape kill operated
DOWNSTREAM (emitted cap-plane images at ~1.7 eps). That level gap IS the point: the
identity emission discards is present, bounded, and order-invariant UPSTREAM; the
escape failed because it worked downstream without it. The skeleton does not close any
carrier and does not claim to - it pressures the joint and hands its numbers to the
risks above.

## 7. Landing history and fence

- Grounding + skeleton evidence: lab notebook .claude/lane-reports/v4impl-stage0-*.
- Skeleton probe: env-gated, run, REVERTED (git diff src/ test/ empty; canonical test
  binary rebuilt from reverted source). Baseline fence: Overlap3.* 54 pass + 1 skip
  (Gate4c_HullMask_MustResolve), unchanged.
- This document is the artifact. Machinery delta from stage 0: ZERO (docs + notebook).
- NEXT: adversarial round-1 lanes attack this plan (the risks in section 4 are the
  entry points), then stage A implementation under the section-3 pins.
