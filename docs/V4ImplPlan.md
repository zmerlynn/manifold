# V4 implementation plan: two-pass global identity emission

STATUS: STAGED IMPLEMENTATION PLAN (gates-first contract), branch
explore/sweep-plane-3d-v4. REVISION 2 - folds three adversarial round-1 lanes
(.claude/lane-reports/v4plan-verify-closure, -gates, -density; all NEED-CHANGE
except the two density SURVIVEs). The converged design is docs/MaintainedEmission3D.md
variant (iii) - two-pass global identity. This document turns that design into a
pinnable, attackable build. It is NOT the pipeline; it is the plan the build lanes
implement under its pins. v3 (docs/SweepEmit3D.md, src/overlap3*) stays intact and
green until a stage DELIBERATELY and PINNEDLY changes it.

REVISION-2 HEADLINE: the project's central unproven premise - does the upstream
identity CLOSE the hardest carriers oracle-true, or only DIAGNOSE the tension and
refuse correctly? - is no longer deferred to the most expensive stage. It is decided
FIRST, cheaply, on the existing skeleton, as STAGE A0 (the SUFFICIENCY PROBE, the
project's go/no-go milestone). Both reasoning lanes converged on this independently.

GATES-FIRST IS LAW (v3's thrice-confirmed lesson). Build agents rescope fail-closed
machinery and weaken gates while narrating happy-path green; that green is
systematically vacuous. The countermeasures that worked, and that every stage below
inherits: contractual gates enumerated BEFORE code; MECHANISM PINS authored red-first
that CANNOT pass against stubs (white-box - they assert internal invariants, not just
outcomes); narrow single-mechanism stages; independent re-runs of every green claim; a
PROGRESS RATCHET (section 3) so "zero oracle-wrong" is not satisfiable by fail-closing
everything.

SETTLED, not re-litigated here: the design verdict and the six-candidate kill table
(MaintainedEmission3D.md); the two structural conditions; the absolute
zero-oracle-wrong posture; the three round-1 lanes' findings; the fence rule (ALL
existing tests are every stage's fences); no counts in docs (magnitudes and by-kind;
exact numbers cite the lab notebooks); ASCII; house-voice commits.

## 0. What the design settled, and what it did NOT

Variant (iii) DISSOLVES the per-plane-independence coupling (v3's many independent
per-plane cap self-locations) and does NOT relocate, under two NAMED structural
conditions that are the design requirements of the build:

  (C1) the near-degeneracy collapse is re-derived from INPUT coordinates via an
       ORDER-INDEPENDENT event ordering (exact-lex event set + exact predicates),
       never consuming a constructed coordinate in a discovery-order-dependent way;
  (C2) each near-degenerate cluster is CANONICALIZED exactly once.

NECESSITY IS PROVEN; SUFFICIENCY IS NOT. The design (MaintainedEmission3D.md
L471-472, L582, L596-598) proves a NECESSITY-ONLY claim: WITHOUT the upstream
arr.vert identities that emission discards, the oracle-true topology is UNREACHABLE.
It does NOT prove the converse - that WITH them, oracle-true IS reachable on the
hardest carriers. The escape kill (probe 3) proved a bounded POSITIONS-ONLY
composition closes Havocglass8 to a VALID genus-0 winding-correct 2-manifold that is
oracle-WRONG by a macro volume (~14000x the eps bound), via a macro-region flip
INSEPARABLE from closure. Variant (iii)'s claim is that the UPSTREAM arrangement
identity emission discards - which faces and seams actually meet at a junction, not
where its images landed - is what disambiguates that flip. That claim is the RESEARCH
FRONTIER; it is UNPROVEN on Havoc, and both round-1 reasoning lanes established that
the deciding information (incidence + signed mult) is NECESSARY-BUT-POSSIBLY-
INSUFFICIENT: the escape-kill inseparability is a TOPOLOGICAL tension incidence may
only let you diagnose, not resolve.

THEREFORE the sufficiency question is answered EARLY, not built around. Stage A0 (the
SUFFICIENCY PROBE) runs against the existing skeleton BEFORE stages A-D are built and
returns one of two outcomes that RE-SCOPE the whole plan:

  - SEPARABLE-AND-CLOSABLE (a computable incidence signature both distinguishes the
    must-survive cluster from the must-collapse twin AND admits a re-emission that
    closes the fan while preserving the lump): stage C acquires its DECISION RULE, and
    the fan carriers (Havoc, GT7863, ...) return to the target list as RED-FIRST
    targets that MUST change state.
  - NOT-SEPARABLE (the escape kill's measured inseparability is the null hypothesis
    and it holds): the Havoc-class is adjudicated HONESTLY FAIL-CLOSED at A0 - cheaply,
    not after building A-D. The target list becomes the SEPARABLE SUBSET, and the
    project's value statement is re-scoped to what remains real: the density wins
    (stage A subsumes the thinning decision, section 1.3/5), match-consistency (the
    per-slab match made combinatorial, PIN 6), any separable carriers A0 identifies,
    and the machinery simplification (dead-guard deletion, stage D).

So the project's DONE criterion is NOT "the corpus resolves", and is NOT the old
loophole "some may resolve". It is: the identity table is built to C1/C2 and carries
the incidence A0 needs; A0's per-carrier verdict is recorded; emission consumes the
table; the PROGRESS RATCHET's named must-move set changes state per A0's outcome; and
the corpus is re-adjudicated with ZERO oracle-wrong resolves. Carriers the
identity-carrying emission cannot close oracle-true stay RECORDED FAIL-CLOSED
contracts, never silent oracle-wrong resolves. The gates protect against regression
and oracle-wrongness regardless of how far the research frontier moves.

## 1. Architecture delta

### 1.1 The two passes

PASS 1 - GLOBAL IDENTITY COLLECT. Runs once, after FindSeams, before BuildSlabs.
CONSUMES the arrangement's junction-defining geometry: the canonical merged verts
(input-quantized, order-invariant), the seams and their endpoints, the criticals
(degenerate-contact x's, seam-seam crossing x's), and the density set. PRODUCES the
IDENTITY TABLE: the named 3D junctions, each junction's cluster COLLAPSED ONCE from
input coordinates WITH its arrangement INCIDENCE attached (section 1.2), and the
density-side thinned critical set. It is a BATCH collect (not incremental - that is
the retired variant i); the emission pass reads it.

PASS 2 - EMISSION. v3's stages, now reading the identity table where they self-located
before. Every cap and every strip that images a junction takes that junction's
canonical position from the table (matched COMBINATORIALLY, PIN 6, never by
position), so no per-plane self-location exists and the assembly weld never reconciles
a junction twin (it cannot arise).

### 1.2 New data structures and their invariants

Round-1 closure lane FACT: the plan's original IdentityTable carried NONE of the
incidence its own closure prose invoked. Operand-of-origin (A vs B) does NOT exist to
carry - Canonicalize (overlap3.cpp:218-305) folds a PRE-MERGED self-intersecting Impl
by SIGNED MULTIPLICITY, discarding the operand label upstream of overlap3. What DOES
remain, and what the table must carry, is the arrangement INCIDENCE:
member arr.verts -> incident seams (Seam.vertId0/1 == member) -> incident faces
(Seam.faceId0/1) -> each face's signed mult (CanonicalFace.mult).

```
struct Junction {
  vec3 pos;                       // canonical position, INPUT-derived (lex-min member)
  std::vector<int> members;       // the arr.vert cluster this junction collapses
  std::vector<int> incidentSeams; // seams with an endpoint in members
  std::vector<int> incidentFaces; // the faceId0/1 those seams join
  std::vector<int> faceMult;      // signed mult of each incident face
  bool spread;                    // classifier: true = one junction, false = distinct
};
struct IdentityTable {
  std::vector<Junction> junctions;
  std::vector<int> vert2junction;         // arr.vert id -> junction id
  std::vector<double> thinnedCriticalXs;  // density-side, reps >= eps apart
};
```

HOW IDENTITY IS CONSTRUCTED (the density lane's measured correction). Identity is NOT
a pure input-side combinatorial key: the order-free face-pair/triple key OVER-COUNTS
the bundle (density lane: 687,332 combinatorial keys on GT7081 where the geometry has
16,889 junctions - it does NOT reach bundle scale). Identity is:
  - constructed-point GROUPING: cluster the arr.verts / contact points by an eps-ball
    union-find over their CONSTRUCTED coordinates. This is measured ORDER-SAFE at
    density - the partition is BITWISE order-invariant (density lane: same 16,889
    clusters under input permutation on GT7081; bitwise on selfA/selfB/Havoc/GT7863);
  - input-derived POSITION: each junction's `pos` is the exact lexicographic-min
    MEMBER (a from-input coordinate), order-free by construction and measured
    order-invariant to within 0.000689 eps worst-case on GT7081 (one junction of
    16,889 wobbles sub-milli-eps), bitwise on selfA/B (density lane).

INVARIANTS, each a pin (section 3):
- `Junction::pos` is a PURE FUNCTION of the member INPUT coordinates (exact
  lexicographic-min member, or an exact-predicate representative), hence ORDER-FREE.
  It is NOT a merge first-match representative - R1 measured that order-dependent at
  the eps level (36/40 nonidentical under vertex permutation, skeleton notebook). This
  is C1's teeth.
- The whole `IdentityTable` - members, incidence, AND positions - is BITWISE-IDENTICAL
  under input permutation (C1). The pin asserts the FULL table, not a summary hash.
- Each junction's `pos` is chosen EXACTLY ONCE during construction (C2) - a
  per-junction counter at the canonicalization site, not a per-consumer/per-member
  re-derivation.
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

DENSITY-SIDE HANDLING lives in PASS 1, and the density lane CORRECTED its mechanism.
The dense critical set is DEGENERATE CONTACTS (near-tangent face pairs, seamLen <=
eps), NOT seam-seam crossing bundles: the density lane measured ZERO M1 seam-seam
triple points on GT7081 (M1=0; a negligible skeleton tail of 15) and on selfA/selfB
(M1=0). The mechanisms still transfer verbatim - degenerate-contact x's are consumed
x-ONLY by SEAMS, exactly as crossing x's were, so the thinning safety argument
(touches only x-consumed criticals; keeps representatives >= eps apart; no macro
feature's bounding vertex moves) holds unchanged. But R-TRIPLE's specific fear
(triple-point completion adds ~1.37M junctions) is MOOT: there are no triple points to
complete on the heavy carriers.

Degenerate-contact clustering produces `thinnedCriticalXs`. It is landed only under
the zero-oracle-wrong posture (section 3, pin 5): thinning that flips a single mesh to
an un-oracle-checkable resolve is FORBIDDEN, and STAGE A SUBSUMES THE THINNING
DECISION (the same eps-ball clustering that names junctions collapses the bundle-
interior criticals to the junction's canonical x - density lane, section 5). This
ABSORBS the open thinning owner-decision left by the v3 wall-A arc (walla-arc
candidate C2, "validated-safe, recorded, not landed"): it is no longer a separate
owner call bolted onto FindSeams, it is a byproduct of the stage-A collect, oracle-
fenced by pin 5.

ARRANGEMENTBUDGET STAYS. It is the honest refusal for genuinely-degenerate-VERTEX
density: GT7081's endpoint clusters stay dense after contact thinning (skeleton:
maxCluster=77 near-degenerate verts survive), so it still refuses. Its firing DOMAIN
SHRINKS: cases whose density is pure contact OVER-INCLUSION (self_intersect A/B, hull)
section after thinning and land on their emission-stage outcome instead. The budget is
not recalibrated; it is the resource backstop that converts a de-facto hang into a
recorded refusal.

### 1.4 v3 stage fates in pass 2 (COMPLETE - every mechanism has a row)

Round-1 gates lane required an explicit fate for EVERY mechanism in
src/overlap3*.{h,cpp} and overlap3_sweep.cpp, so a stage-D agent deleting "dead"
machinery cannot improvise. The two previously-unaccounted primaries
(SplitTouchingSheets, the section-vert->junction match) and two loose ends (coplanar
family, chain-plane loX/hiX) are now rows.

| v3 mechanism (file:site)                | fate under variant (iii)                 |
|-----------------------------------------|------------------------------------------|
| Canonicalize (overlap3.cpp:218)         | SURVIVES; feeds pass 1. Signed-mult fold  |
|                                         | is why operand-of-origin is absent (1.2).|
| FindSeams (:323)                        | SURVIVES; its verts/seams/criticals are   |
|                                         | pass-1 INPUT; contact criticals thinned.  |
| BuildSlabs (overlap3_sweep.cpp:95)      | SURVIVES as the per-critical collect       |
|                                         | KERNEL; its section verts BIND into the   |
|                                         | global identity by the COMBINATORIAL match|
|                                         | (PIN 6), taking the junction's canonical. |
| SlabResolver::Extend self-location      | Self-location (Extend via InterpolateSafe)|
| (overlap3.cpp:519)                      | DELETES for junction-imaging verts - the  |
|                                         | wall-A locus. Degenerates to a table read;|
|                                         | non-junction endpoints keep track subdiv. |
| section-vert -> junction MATCH (NEW,    | pass-2 mechanism: which section vert is    |
| pass 2)                                 | junction J's image. MUST be COMBINATORIAL |
|                                         | (edge/seam incidence to a member of J),   |
|                                         | once-only, cross-slab consistent - PIN 6. |
|                                         | A positional match RELOCATES the wall.    |
| ComputeCap / BuildCapEdgeSet (:785/821) | CONSUMES identity (stages B, C).          |
| EmitCaps pair driver (ci == li+1, :918) | CONSUMES identity; pair-canonical cap      |
|                                         | emission (stage C). Post-mod form pinned. |
| EmitStrips + provenance chains          | CONSUMES identity; provenance chains       |
| (StripChains, :662)                     | (edgeSubdiv) UNCHANGED (stage C).         |
| chain-plane loX/hiX placement           | constructional-closure cap-plane placement|
| (StripChains, :664/700/943)             | (Coplanar_PerpFacesSubEps money fixture). |
|                                         | RE-ADJUDICATE at C/D: name-joined strips  |
|                                         | may unify loX/hiX; KEEP until pinned dead.|
| chain-plane wide-run rule + SubEps       | RE-ADJUDICATE at D. Junction identity      |
| wide-run cap backstop (:860/871)        | replaces run-spanning interpolation; the  |
|                                         | wide-run backstop MAY die. KEEP until      |
|                                         | pinned dead; do not delete on narration.  |
| SubEpsFeature single-face coverage       | KEEP (stage D). A macro face wholly in a  |
| guard (BuildSlabs)                      | merged ulp-run is a genuine hazard the    |
|                                         | table does not touch.                    |
| SplitTouchingSheets (overlap3.cpp:985,  | KEEP. Called by BuildImpl (:1175) AFTER   |
| called at :1175 in BuildImpl)           | weld/dedup, BEFORE Is2Manifold. It is the |
|                                         | wall's PROXIMATE failure site (unpaired   |
|                                         | twin cap sheet + micro-edge). Under        |
|                                         | identity-aware assembly its INPUTS change:|
|                                         | junction images share the canonical, so   |
|                                         | the twin sheet no longer ARISES to be     |
|                                         | paired. It still runs UNCHANGED for       |
|                                         | genuine touching contacts (Touch_EdgeEdge,|
|                                         | Touch_VertexOnly fences). NOT deleted; the|
|                                         | twin dissolves UPSTREAM at C, verified by |
|                                         | pin 4b, not by editing this function.     |
| coplanar-group family (face2Group,      | SURVIVES under "BuildSlabs survives", but |
| GroupedEdge, CollectGroupedMemberEdges, | with an EXPLICIT contract: a coplanar-    |
| SlabResolver class-iii planar,          | group junction enters the identity table  |
| overlap3.cpp:159/355/480/522)           | like any other; class-iii planar          |
|                                         | resolution reads the canonical, not a     |
|                                         | self-located planar image. Fences the     |
|                                         | ~10 Coplanar_*_Oracle tests. Pinned at C. |
| BuildImpl / assembly weld hash grid     | Assembly weld is a no-op for junction      |
| (:1107, cell=eps)                       | images (they share the canonical), so the |
|                                         | twin cannot arise. KEPT.                  |
| Is2Manifold gate (BuildImpl :1189)      | KEPT - the final honesty check.           |
| ArrangementBudget                       | KEEP; domain shrinks after thinning (5).  |
| EngineIdConflict / CoplanarOverlap /    | already RETIRED (v3).                     |
| EdgeInPlane                             |                                          |

DIES: per-plane independence of junction self-location; track extension AS the
junction LOCATOR; assembly-weld reconciliation of junction twins. BORN: the identity
table WITH incidence (pass 1); the combinatorial junction-identity match/read in
caps/strips/assembly (PIN 6).

## 2. Stage sequence

Each stage is narrow (single mechanism), additive-first (new channels default-off /
unconsumed until the consuming stage lands - the negEdges/SweepCapture precedent),
independently green, and independently re-run. For each: mechanisms, files, v3 touch,
CONTRACTUAL GATES (its fences = ALL existing tests, plus its new white-box pins), its
red-first anchor, and its DONE criterion. Every stage DONE below is a NAMED-SET
observable criterion checkable by a non-author (no judgment words, no permissive
parentheticals).

### Stage A0: the SUFFICIENCY PROBE (project go/no-go, FIRST milestone)

WHY FIRST. The project's central premise - upstream incidence CLOSES the hardest
carriers, not merely DIAGNOSES the tension - is genuinely open, and both round-1
reasoning lanes proved it is decidable CHEAPLY, NOW, against the existing stage-0
skeleton, without building A-D. Deferring it to stage C (the old plan) builds the
whole architecture to learn something a throwaway probe settles. A0 is env-gated,
throwaway (skeleton-style, reverted after capture), and gates the rest of the plan.

MECHANISM. On the stage-0 skeleton (RemoveOverlaps3D_TestHooks arr.verts), build the
identity table AND attach the incidence of section 1.2 to each junction (for each
member: incident seams via Seam.vertId0/1 == member; the faces those seams join via
Seam.faceId0/1; each face's signed CanonicalFace.mult). Run it on two calibrated
clusters on Havocglass8, plus GT7863 as a second-anatomy data point:

  - THE MUST-SURVIVE cluster (J_lump): the junction whose member sits at the
    right-operand vertex the escape kill lost the ~336-unit lump at - the Havocglass8
    x ~ 12940-12942 cluster whose two cap images are the ~1.7-eps twin (walla-arc
    step 3: v4/v509, v5/v510, 3D dist 3.83e-8 = 1.70 eps at eps=2.25e-8, same y to
    1e-8, split in z ~1.6 eps and x ~2.4e-10).
  - THE MUST-COLLAPSE control (twin): the sub-eps-x 2-vert cluster the skeleton
    already collapses clean (dx = 0.445 eps, classifier fires, zero fence violations -
    v4impl-stage0 notebook).
  - GT7863 as a SECOND A0 data point: its in-run macro content is a DIFFERENT anatomy
    - the M4 cancelled-edge macro content bracketed between the twin's two cap planes
    at a ~1.12-eps x-residual (closure lane attack 4). Separability may differ per
    carrier, so A0 measures each carrier it will later target.

THE THREE TESTS A0 ANSWERS (recorded in its notebook, decisive):
  (i) SEPARABILITY: is there a COMPUTABLE incidence signature (over incidentSeams /
      incidentFaces / faceMult, order-free from input) that DISTINGUISHES J_lump
      (must survive) from the collapsible control twin (safe to collapse)? Operand-of-
      origin does NOT exist in the input (closure lane FACT 1); the signature is
      incidence+mult INFERENCE - does J_lump's member bound real faces with nonzero
      mult that the control's does not? The escape kill's measured inseparability is
      the NULL HYPOTHESIS.
  (ii) IF SEPARABLE, the DECISION-RULE SKETCH: can a re-emission keyed on that
      signature close the fan AND preserve the lump (a paper walk plus the smallest
      runnable probe)? The escape-kill inseparability (every config that keeps a
      lump-region cluster distinct REOPENS the fan) predicts this is FAIL-CLOSED, not
      oracle-true closure. A0 refutes or confirms that prediction before C is built.
  (iii) THE BLANK-TABLE MUTATION for the isolated-twin pin (pin 4b, section 3): with
      the incidence table BLANKED, does the isolated twin still close positions-only?
      The escape kill's positive result says YES - which is exactly why pin 4b must
      carry the blank-table mutation, and A0 measures the baseline the mutation must
      flip.

OUTCOMES (written into this plan by the A0 lane, re-scoping stages C/D):
  - SEPARABLE-AND-CLOSABLE: stage C acquires the decision rule A0 sketched; the fan
    carriers A0 cleared return to the PROGRESS RATCHET as red-first must-move TARGETS.
  - NOT-SEPARABLE (the predicted outcome): the Havoc-class is adjudicated HONESTLY
    FAIL-CLOSED at A0. The target list becomes the separable subset (if any); the
    plan's value statement is re-scoped to density wins + match-consistency + any
    separable carriers + machinery simplification (section 0). C's scope for a
    NOT-SEPARABLE carrier is the honest fail-closed contract, and B/C are NOT built to
    chase a closure A0 already ruled out.

FILES: env-gated hooks on the existing skeleton path (test/overlap3_test.cpp), NO
production machinery. V3 TOUCH: none (reverted after capture, like stage 0).
GATES: the skeleton fence (Overlap3.* baseline) unchanged; A0 changes zero production
behavior. DONE: the three tests answered and RECORDED in the A0 notebook; the
SEPARABLE / NOT-SEPARABLE verdict written back into stages C/D and the ratchet; probe
reverted (git diff src/ test/ empty). A0 is the go/no-go: stages B/C are scoped by its
verdict.

### Stage A: identity table construction (pass 1)

MECHANISM: a new `BuildIdentityTable(arr, eps)` between FindSeams and BuildSlabs.
Cluster arr.verts by an eps-ball union-find (from-input positions, exact predicates);
per cluster apply the classifier; choose the canonical ONCE (coordinate-determined);
attach the incidence (section 1.2). SCOPE: vertex junctions (arr.verts collapsed),
each carrying its arrangement incidence. Resolving seam-seam crossings to 3D triple-
point junctions is NOT needed - the density lane measured ZERO triple points on the
heavy carriers (M1=0), closing R-TRIPLE (section 4). The case-1 twin's underlying
junction V is itself an arr.vert.
FILES: overlap3.cpp (new function), overlap3.h (Junction + IdentityTable,
ArrangementGeometry gains an `IdentityTable identity` field default-empty), test hooks
expose it.
V3 TOUCH: ADDITIVE. ArrangementGeometry gains a field; emission IGNORES it. Default
behavior BITWISE unchanged.
GATES: ALL existing tests (fence) + Pin_IdentityTable_OrderInvariant (full table
bitwise), Pin_IdentityTable_OnceOnly (global once-only counter), Pin_Classifier_Fence
(section 3, pins 1-3).
RED-FIRST ANCHOR: a constructed fixture with a known near-degenerate cluster (or
Havoc's) whose junction AND incidence the table MUST name; the pin reds against an
empty / per-vertex-identity / incidence-blank stub.
DONE: table built with incidence + three pins green (mutation-verified) + suite green +
emission bitwise-unchanged (a diff of emitted geometry against HEAD is empty).

### Stage B: caps consume junction identity (transverse)

MECHANISM: ComputeCap / BuildCapEdgeSet reads a junction-imaging section vert's
canonical (y,z) from the identity table instead of the SlabResolver's self-located
extension. The section vert is matched to its junction COMBINATORIALLY (PIN 6), not by
proximity. Both caps of one junction then place it at the SAME (y,z), dissolving the
transverse divergence. A lookup with fallback to the resolver for non-junction verts
(additive-shaped).
FILES: overlap3.cpp SlabResolver::Extend (identity-first path), BuildCapEdgeSet.
V3 TOUCH: the SlabResolver gains an identity-first branch; ComputeCap reads the table.
GATES: ALL existing + Pin_CapImage_CanonicalTransverse (pin 4a): a junction's two cap
images have BIT-IDENTICAL (y,z) EQUAL to the table's canonical yz (white-box table
read), mutation-verified against self-location AND against a blanked table; + PIN 6
(combinatorial match, cross-slab consistent); + the STRENGTHENED B FENCE below.
RED-FIRST ANCHOR: the isolated case-1 twin fixture (exact coordinates pinned in pin
4b) - the transverse bind reds (yz diverges) without the table, greens (yz collapses
bitwise to the canonical) with it, reproducing GT7863's R2 bound trace.
STRENGTHENED B FENCE (closure lane required this): B moves EVERY junction-imaging vert
corpus-wide from self-located to canonical yz, and the design's own residual-1
(MaintainedEmission3D:205-215) warns that perturbing a cap arrangement can flip a
retained region. So B's fence is NOT "tests still pass + diff confined": it is
BITWISE / ORACLE-INVARIANCE on EVERY currently-resolving carrier (the named ratchet
set: Cray, Offset2/3/4, all currently-green synthetic + Coplanar gates), matching
stage A's discipline. The emitted-geometry diff must be EMPTY except at the enumerated
junction-imaging verts (a NAMED set), bitwise - not "confined to the transverse bind"
by author judgment.
HONESTY / the split from C: stage B alone does NOT close the twin. The cap-plane
x-residual survives (GT7863 bound-and-still-diverged: yz bitwise, ~1.1-eps x-gap). B's
DONE is the transverse-identity pin + PIN 6 + the strengthened fence, NOT twin closure.
DONE: cap images bit-identical to the table's canonical yz (mutation-verified, incl.
blank-table) + PIN 6 green + every ratchet-set carrier bitwise/oracle-invariant + the
emitted-geometry diff empty outside the named junction-imaging vert set + suite green.

### Stage C: strips + assembly join by name (coordinated re-emission)

MECHANISM: the junction reaches the EMITTED vertex. A junction is emitted at ONE 3D
point - including its canonical x, not each cap plane's own x - so the two adjacent
cap planes imaging one junction UNIFY into one plane; strips join by junction NAME
(provenance chains supply the subdivision between junctions); the assembly weld is a
no-op for junction images. This is the coordinated cap+strip re-emission around a
near-degenerate cluster the wall-A memo named as the real fix.
THE DIFFERENTIATION PARAGRAPH (what C does that killed probe 4 / the escape kill did
NOT - stated as data + decision rule, not vocabulary). Round-1 confirmed stages B+C
COMPOSED are the killed probe-4 composition (OV3_IDBIND + OV3_CAPMERGE on GT7863,
outcome (iii) dead zone) and that stage C's collapse+re-mesh mechanism IS the escape
kill's. Two things and only two make C's OUTCOME differ:
  (1) the ZERO-ORACLE-WRONG GATE (pin 5): where the escape kill closed oracle-WRONG,
      C REFUSES (recorded fail-closed). This is a GATE, legitimately fencing the
      failure mode - not a new closing mechanism. It is why C is honest, not why C
      closes anything.
  (2) A0's DECISION RULE, IF A0 returned SEPARABLE: the incidence signature that marks
      J_lump "distinct/preserve" vs the control "collapse", plus the re-emission A0
      proved closes-the-fan-while-preserving-the-lump. This is the ONLY thing that
      can make C CLOSE an entangled carrier rather than refuse it, and it EXISTS only
      if A0 delivered it. The load-bearing obstruction probe 4 hit - the fate of the
      REAL macro content bracketed between a twin's two cap planes (collapse it ->
      SubEpsFeature backstop; keep it -> inconsistent macro face) - is resolved by
      A0's rule or it is not resolved at all. If A0 returned NOT-SEPARABLE for a
      carrier, C's scope for that carrier is the honest fail-closed contract, and C
      does NOT reconstruct probe 4 hoping for a different outcome.
FILES: overlap3.cpp EmitCaps/ComputeCap (cap-plane unification for junction-imaging
pairs), EmitStrips, BuildImpl.
V3 TOUCH: pair-canonical cap emission; strip chains; assembly weld.
GATES: ALL existing + Pin_IsolatedTwin_ClosesOracleTrue (pin 4b, exact coords + blank-
table mutation) + Pin_ZeroOracleWrong meta-gate (pin 5, with the CorpusSingleGate
rewrite it requires) + the PROGRESS RATCHET (section 3): the named must-move set
changes state per A0's verdict; the named resolving set stays oracle-true/bitwise.
RED-FIRST ANCHOR: the isolated case-1 twin fixture reds (NonManifoldEmission) without
C, resolves ORACLE-TRUE (volume + genus + winding) with C, and REOPENS under the
blank-table mutation - the escape's proven bounded partial win, now landed as the
mechanism AND pinned load-bearing on the table.
DONE (named-set, not permissive): isolated twin oracle-true (and blank-table-reopens)
+ the ratchet's A0-named must-move carriers changed state as A0 predicted (or, if A0
returned NOT-SEPARABLE, the must-move set is empty by A0's recorded verdict and the
Havoc-class carriers are recorded fail-closed) + zero oracle-wrong across the corpus +
suite green. "Some may resolve" is REMOVED; the carrier state-changes are named by A0.

### Stage D: guard dissolution + corpus re-adjudication + density + perf

MECHANISM: delete the now-structurally-dead machinery and re-adjudicate the guard
family; land degenerate-contact thinning (density-side of pass 1) gated by the oracle
posture; re-adjudicate ArrangementBudget's domain; perf.
PER GUARD (delete at the stage that kills it, or keep the fail-closed posture):
- SlabResolver junction self-location: DEAD at C (junctions read the table). Delete,
  pin-dead (a mutation restoring self-location reds the twin pin).
- Chain-plane wide-run rule + SubEpsFeature wide-run cap backstop: RE-ADJUDICATE.
  Under identity emission the run-spanning interpolation is replaced by junction
  identity, so the wide-run backstop MAY become dead. KEEP until pinned dead; do not
  delete on narration.
- Chain-plane loX/hiX constructional-closure placement: RE-ADJUDICATE (may unify under
  name-joined strips). KEEP until pinned dead; Coplanar_PerpFacesSubEps is its fence.
- SubEpsFeature single-face coverage guard (BuildSlabs): KEEP. A macro face wholly in
  a merged ulp-run is a genuine merged-run hazard the identity table does not touch.
- ArrangementBudget: KEEP; domain shrinks after thinning (section 5).
- Is2Manifold gate (BuildImpl): KEEP - the final honesty check.
- SplitTouchingSheets: KEEP (touching-contact fences); verify junction twins no longer
  reach it (pin 4b covers the dissolution upstream).
- EngineIdConflict / CoplanarOverlap / EdgeInPlane: already retired (v3).
FILES: overlap3.cpp (guard deletions + thinning in pass 1), overlap3_sweep.cpp
(budget domain), tests (guard-dead pins + corpus re-adjudication + tol-invariance).
GATES: ALL existing + guard-dead pins + the corpus contracts + the progress ratchet +
perf checkpoints (section 5) + Pin_ZeroOracleWrong.
RED-FIRST ANCHOR: each guard-dead pin reds if the guard's mechanism is restored;
thinning's tol-invariance pin reds if a single-mesh resolve moves under tolerance.
DONE: dead guards deleted and pinned-dead; thinning landed+oracle-fenced (a single-
mesh resolve stays PROVISIONAL, never counted as a win - pin 4d) OR deferred with its
tol-invariance evidence recorded; ArrangementBudget domain re-adjudicated with a named
per-carrier before/after (not the word "re-adjudicated"); perf within bound; corpus
zero-oracle-wrong with the ratchet's named set intact.

## 3. The pin list (authored as specifications now)

Each pin is white-box (asserts an internal invariant), red-first, and carries a
MUTATION-VERIFICATION recipe: a named stub the pin MUST red against, verified by the
author independently re-running the red (not trusting a build agent's green). Round-1
gates lane found 4/5 original pins stub-passable and a fail-closed loophole; the
hardened specs below close each hole named pin-by-pin.

PIN 1 - Pin_IdentityTable_OrderInvariant (C1). Permute the input (vertex + triangle
reorder, several seeds), rebuild the table, assert the WHOLE table is BITWISE-identical
- members, INCIDENCE (incidentSeams / incidentFaces / faceMult), AND positions (exact
double compare) - NOT a summary hash and NOT only the sorted junction set. STUB HOLE
CLOSED: the original pin tested only the OUTCOME (bitwise on small fixtures) and passed
for LexMin-over-CONSTRUCTED-coords; tighten to assert the canonical is recomputable
from INPUT-vertex ids, independent of the arr build. WHY A STUB CANNOT PASS: a merge
first-match representative is order-dependent at the eps level (R1: 36/40 nonidentical
under vertex permutation). MUTATION: replace the canonical choice with `members[0]`
(discovery-order first) -> reds under permutation. SKELETON + DENSITY EVIDENCE: bitwise
on control + Havoc (skeleton); bitwise partition + position within 0.000689 eps on
GT7081, bitwise on selfA/B (density lane) - the R-ORDER-UNMEASURED-AT-DENSITY gap now
CLOSED for the grouping and near-closed for the position.

PIN 2 - Pin_IdentityTable_OnceOnly (C2, counter). A GLOBAL once-only counter on the
REPRESENTATIVE-SELECTION PRIMITIVE, asserted == junctionCount across the WHOLE pipeline
run (construction AND emission), at the canonicalization SITE, wrapper-proof. STUB HOLE
CLOSED: the original per-junction construction counter caught "moved the only
canonicalization out" but NOT "ADDED a redundant re-derivation in a consumer" - the
global primitive counter catches both. WHY A STUB CANNOT PASS: a stub that
canonicalizes inside the cap/strip loop re-derives the representative per consumer ->
global counter > junction count. MUTATION: move the canonical choice into the emission
loop, OR add a per-consumer re-derivation -> counter reds either way.

PIN 3 - Pin_Classifier_Fence (both directions). On a fixture carrying BOTH the case-1
spread twin AND the case-2 DISTINCT pair, assert the classifier FIRES on the twin
(collapses to one) and EXCLUDES the distinct pair (never clustered nor collapsed). STUB
HOLE CLOSED: the original fixture's distinct pair was "tens of eps apart in x" -
excluded by the dx threshold ALONE, so a dx-threshold-only stub (no transverse test)
passed. The fixture MUST ALSO carry a SUB-EPS-X, TRANSVERSE-SEPARATED pair (maxDx
small, maxTr < maxDx) that an eps-ball merges but the `maxTr >= maxDx` conjunct
excludes. MUTATION: drop the `maxTr >= maxDx` conjunct -> the sub-eps-x transverse pair
merges, pin reds (genuinely exercising the load-bearing conjunct); OR widen dxThresh
past the far distinct pair's dx -> that pair merges, pin reds. This is the R4
verdict-flipper's fence made a pin, now with a fixture that reds BOTH mutations.

PIN 4 - per-carrier targets.
  (4a) Pin_CapImage_CanonicalTransverse (stage B): a junction's two cap-input images
       have (y,z) BIT-IDENTICAL to identityTable.junction(v).pos.yz (WHITE-BOX table
       read, not merely equal to each other). STUB HOLE CLOSED: "both images equal" is
       passed by a global-snap stub that never reads the table; asserting equality to
       the TABLE'S canonical + the blank-table mutation closes it. MUTATION: restore
       SlabResolver self-location -> the two images diverge, pin reds; AND blank the
       table -> the pin cannot find the canonical, reds.
  (4b) Pin_IsolatedTwin_ClosesOracleTrue (stage C): a hand-built isolated case-1 twin
       resolves to a clean 2-manifold, volume-neutral, oracle-true (a+b Boolean).
       EXACT FIXTURE COORDINATES PINNED IN THIS PLAN (not author-chosen, to defeat
       detect-and-dodge): two wedge faces whose section produces ONE junction V at
       x0=12940.0, sub-eps in x from two bracketing criticals at x0 +/- 0.30 eps, on a
       steep (~100x) extension track that amplifies the sub-eps x-gap to a ~1.7-eps
       transverse image split (z-dominated, matching walla-arc's Havoc twin: dist ~1.70
       eps, same y, split in z), eps = 2.25e-8. The two cap images land at
       (x0, y, z) and (x0, y, z + 1.7*eps). MUTATIONS (BOTH required): (1) disable
       cap-plane unification -> NonManifoldEmission, pin reds; (2) BLANK THE IDENTITY
       TABLE -> the twin REOPENS. Mutation (2) is load-bearing: without it the pin
       validates positions-only closure the escape kill ALREADY proved works WITHOUT
       the table (probe 3's positive result), certifying nothing about the table.
  (4c) PAIRS get EXACT Boolean oracles (volume + genus + winding): Havoc, GT7863,
       Cray, GT7081, hull. A resolve MUST match; a fail-closed MUST be a named guard.
  (4d) SINGLES get TOL-INVARIANCE contracts (resolve volume identical across a
       tolerance sweep, verified by EXACT bit comparison, not a fixed-precision
       print) + 2-manifold validity: Offset1-4, openscad, self_intersect A/B. A single
       mesh is NEVER asserted resolved-correct without an oracle; a tol-invariant
       single-mesh resolve with no a+b oracle is labeled PROVISIONAL and NEVER counted
       as a win (pin 5's rule; tol-invariance proves not-a-thinning-artifact, not
       geometric correctness).

PIN 5 - Pin_ZeroOracleWrong (meta-gate, absolute). Across the WHOLE corpus + synthetic
suite, in every mode, ZERO oracle-off resolves. Any resolve is oracle-checked (pairs)
or tol-invariant-validity-checked+PROVISIONAL (singles); any carrier that cannot be
checked stays FAIL-CLOSED. WHY IT EXISTS: the escape kill closed Havoc to a VALID
manifold that was oracle-WRONG by a macro volume. Validity is NOT correctness. This
meta-gate forbids exactly that outcome. HARNESS CHANGE NAMED (gates lane): the current
CorpusSingleGate (test/overlap3_test.cpp ~L1619-1637) checks single-mesh resolves on
VALIDITY ONLY (ours.Status()==NoError) - NO oracle, NO tol-invariance. Pin 5 / pin 4d
are NOT enforced until CorpusSingleGate is REWRITTEN to run the tolerance sweep + exact
bit compare and to label the result PROVISIONAL. A build agent claiming "pin 5 green"
WITHOUT that rewrite is vacuous; the orchestrator verifies the rewrite landed.

PIN 6 - Pin_JunctionMatch_Combinatorial (NEW; all three round-1 lanes converged on
this as the plan's most improvised-around gap). The per-slab section-vert -> junction
MATCH - "is this section vert junction J's image?" - is a NEW per-slab, per-consumer
decision pass 2 introduces (density lane R-ONCE-VS-STRICTFP). It MUST be COMBINATORIAL
(the section vert lies on a track/edge whose defining edge or seam is INCIDENT to an
arr.vert that is a MEMBER of J - read from the table's incidence), NEVER positional
(within eps of J's canonical). WHY: a POSITIONAL match can DIVERGE across two adjacent
slabs (one matches J, the other falls through to track-extension), RELOCATING the
variant-(ii) "which section vertex is the image of V" kill (MaintainedEmission probe 1)
into the match. The pin asserts (1) the match predicate is combinatorial (reads
incidence, not position) and (2) CROSS-SLAB CONSISTENCY: the two adjacent slabs imaging
one junction assign the SAME junction id to their shared-junction section verts, by
construction. It is ONCE-ONLY (folded into pin 2's global counter). The gates lane
counted this bind ASSERTED six times and MECHANIZED zero; PIN 6 mechanizes it.
MUTATION: replace the combinatorial predicate with `la::length(pt - J.pos) <= eps`
(the v3 positional resolver match) -> cross-slab consistency reds on the twin fixture
(one slab matches, the other extends).

THE PROGRESS RATCHET (replaces the fail-closed loophole; gates + closure lanes). "Zero
oracle-wrong" is a REGRESSION guard, not a progress gate - it is satisfiable by
FAIL-CLOSING everything (the current corpus already fail-closes Havoc/GT7863 ->
NonManifoldEmission, GT7081 -> ArrangementBudget). So the DONE criteria pair it with a
ratchet:
  - THE NAMED RESOLVING SET stays oracle-true / bitwise EVERY STAGE: Cray (a+b oracle
    rel=0), Offset2/3/4 (bitwise volume), all currently-green synthetic + Coplanar
    gates. A stage that regresses any of these is NOT done. This is the floor.
  - THE NAMED MUST-MOVE SET is written by A0: stage C/D's DONE names the EXACT fixture
    set that MUST change state, per A0's verdict. If A0 returned SEPARABLE for a
    carrier, that carrier is a must-move target and C is not done until it moves
    oracle-true. If A0 returned NOT-SEPARABLE, the must-move set is EMPTY by A0's
    recorded verdict, the Havoc-class carriers are recorded fail-closed, and the
    plan's honest deliverable is density wins + match-consistency + machinery
    simplification (section 0) - NOT "some may resolve". Either way the DONE names the
    carrier state-changes; "some may resolve" is removed.

MUTATION-VERIFICATION DISCIPLINE (every pin): the author re-runs the red against the
named stub independently. A green claim from a build agent is re-run before it is
believed. Every "mutation-verified" DONE depends on the orchestrator INDEPENDENTLY
re-running the mutation - it lives outside the build agent's gate text. This is v3's
thrice-confirmed countermeasure, not optional ceremony.

## 4. Risks (round-1 outcomes folded)

Named for a skeptic to break, with the round-1 verdict where the lanes settled one.

R-DENSITY (the identity table's size on GT7081's dense criticals). SETTLED SURVIVE
(density lane): folding the FULL degenerate-contact set into junction identities keeps
the table bounded at ~vertex scale (GT7081: 1,374,649 criticals -> 16,889 junctions,
tracking nVerts ~18.6k, 9.2s / 134MB; selfA 124,142 -> 8,546). The table's SIZE is not
the density risk; the criticals feeding BuildSlabs are (v3's existing
ArrangementBudget/thinning problem, addressed by the stage-A thinning subsumption,
section 1.3/5).

R-TRIPLE (must pass 1 resolve triple points?). SETTLED SURVIVE, the feared object is
ABSENT (density lane): the heavy carriers carry ZERO M1 seam-seam triple points
(GT7081 M1=0, selfA/B M1=0; a negligible skeleton tail of 15 on GT7081). The ~1.37M
"criticals" are 99.99%+ DEGENERATE CONTACTS (near-tangent face pairs), not crossings.
The triple-point completion that would explode the table does not exist to be resolved;
the plan's vertex-junction scoping is VINDICATED. The escape-kill disambiguation on the
dense endpoint cluster (skeleton maxCluster=77) is what A0 tests - the residual risk is
sufficiency (R-CLOSURE), not table size.

R-ONCE-VS-STRICTFP (does canonicalize-once conflict with strict-FP per-slab
sectioning?). SETTLED NEED-CHANGE -> ADDRESSED by PIN 6 (density lane Q3). C2 covers
the junction POSITION (canonicalized once); the GAP is the CONSUMPTION-side per-slab
MATCH pass 2 adds. If positional it can diverge cross-slab and relocate the wall; PIN 6
requires it combinatorial + cross-slab-consistent + once-only. Mild (a spec + one pin),
not a break: the table's C2 is sound, the match is the residual.

R-BC (is stage B green without C, or the split fictitious?). SETTLED (closure lane):
the DECOMPOSITION is anatomically REAL (transverse yz collapses at B, cap-plane
x-residual at C = the GT7863 trace) - NOT fictitious. But B's original fence was WEAKER
than A's against the design's residual-1 (moving V's image to canonical yz can flip a
neighboring retained region). ADDRESSED: B's fence is strengthened to bitwise/oracle-
invariance on every ratchet-set carrier (stage B, section 2), matching A's discipline.

R-CLOSURE (the plan assumes a closure it has not earned). SETTLED (both reasoning
lanes): necessity proven, sufficiency UNPROVEN, and the deciding experiment is now
STAGE A0 (moved from stage C). The escape kill proved positions-only closes Havoc
oracle-WRONG; whether upstream incidence reaches oracle-true is what A0 measures BEFORE
B/C are built. The plan does NOT gate on "Havoc resolves"; the DONE is A0's verdict +
zero-oracle-wrong + the ratchet's named set. The incidence is NECESSARY-BUT-POSSIBLY-
INSUFFICIENT (the tension is topological, not informational); A0 settles which cheaply.

R-THINNING-ORACLE (thinning flips single meshes with no oracle). Thinning turns
self_intersect A/B from ArrangementBudget refusal to a resolve with NO a+b oracle. The
plan defers landing thinning until an oracle exists and substitutes a tol-invariance
contract - LABELED PROVISIONAL, never a win (pin 4d/5). Tol-invariance proves the
resolve is not a thinning artifact; it does NOT prove geometric correctness, so it is
validity-plus, not oracle. Honest by the PROVISIONAL labeling.

R-ORDER-UNMEASURED-AT-DENSITY (C1 bitwise-ness at density). SETTLED near-closed
(density lane): the GROUPING is BITWISE order-invariant at density (GT7081 same 16,889
clusters under permutation; bitwise on selfA/B/Havoc/GT7863); the POSITION is
order-invariant to within 0.000689 eps worst-case on GT7081 (one junction of 16,889),
bitwise elsewhere. The tiny position residual is exactly what C1's exact-predicate
re-derivation from input zeroes. The skeleton-only gap the original plan flagged is now
measured at density for the contact set.

R-MATCH-BIND-WELD (closure lane attack 2, the weld-image pole). Havoc's blocking
junction images include WELD images (naming NO canonical arr.vert) and FAR-track
interpolations (neither is an arr.vert, so vert2junction cannot map them). The
combinatorial match (PIN 6) reaches section verts with arr.vert/seam provenance;
a pure 2D-crossing weld image with NO arr.vert provenance is the never-bound case. If
A0 returns NOT-SEPARABLE this is subsumed (Havoc-class fail-closed); if A0 returns
SEPARABLE, C's decision rule must state what junction id (if any) such a vert takes, or
record it as the structural fail-closed boundary. This is the section-level face of the
spread-vs-distinct ambiguity; A0 is where it is first measured.

## 5. Perf checkpoints

BOUND: the full suite within ~2x of the current baseline, and no single test beyond the
~120s single-test ceiling, checked PER STAGE (a stage that regresses perf is not done).

- Stage A0: throwaway skeleton probe, env-gated, adds nothing to the default suite.
- Stage A: table build is O(nVerts) eps-ball clustering with incidence attach -
  skeleton measured 0.02s on ~18.6k verts; the density lane's contact-folding clustering
  cost 9.2s / 134MB on GT7081's 1.37M contacts, sub-second on selfA/B - all one-time,
  bounded. The order-invariance PIN runs several permutations x hooks; keep it on SMALL
  fixtures only (control + Havoc) and do NOT permute GT7081 (its hooks are ~30s each).
  CHECKPOINT: stage A adds a low single-digit percentage to the suite.
- Stages B/C: the identity read is an O(1) map/hash lookup per section vert; the
  combinatorial match (PIN 6) is an incidence lookup. Should not regress. CHECKPOINT: no
  test beyond the single-test ceiling; suite within 2x.
- Stage D / the density interplay once contacts collapse: degenerate-contact thinning
  REDUCES the critical set an ORDER OF MAGNITUDE (density lane, MEASURED on the critical
  set: GT7081 slabs 231,616 -> 19,530 under eps-merge, ~12x; selfA 63,089 -> 8,942,
  ~7x), which SPEEDS BuildSlabs and lets the contact-over-inclusion cases (hull,
  self_intersect) section. The PROJECTED budget relief (thinned crits could flip GT7081
  from budget-refusal to sectionable) is a NAMED STAGE-A FOLLOW-UP MEASUREMENT (the
  load-bearing BuildSlabs-on-thinned-crits run), NOT assumed - the density lane
  projected but did not run it. GT7081's ENDPOINT clusters stay dense (skeleton
  maxCluster=77), so ArrangementBudget still fires there and is NOT recalibrated up
  (relaxing it runs the dense case to the emission blowup and STILL fails, strictly
  worse). CHECKPOINT: after thinning, hull sections and lands on its cap-plane outcome
  (resolve if C closed it per A0, else recorded fail-closed); self_intersect resolves
  tol-invariant PROVISIONAL (oracle-fenced or deferred); GT7081 still refuses
  ArrangementBudget within seconds and bounded memory; the suite stays within bound with
  the budget as backstop.

## 6. Walking-skeleton probe evidence

A thin real version of the load-bearing joint - the identity table - built on v3's
existing arrangement (arr.verts from RemoveOverlaps3D_TestHooks) and run on one small
fixture (Havocglass8) and one clean control (Gate4a kWedges8), plus density
measurements on GT7081 and the self-intersection pair. Env-gated
(OV3_V4SKELETON[_GT7081]), default-skip, REVERTED after capture (machinery delta zero);
the recipes and full numbers live in the lab notebooks (.claude/lane-reports/
v4impl-stage0-* for the skeleton; v4plan-verify-density-* for the contact-set density).
The table clusters arr.verts by an 8*eps ball, applies the probe-3 classifier (maxDx <=
0.7*eps AND maxTr >= maxDx), and chooses the canonical as the lexicographic-min member
(coordinate-determined).

C1 (order-invariance), skeleton, permute input x8 and compare junction sets bitwise:
- Control (kWedges8): BITWISE identical (worst nearest-deviation 0 eps, no count
  difference, no exact miss).
- Havocglass8: BITWISE identical (worst deviation 0 eps). One collapse fired (a 2-vert
  cluster, dx ~ 0.45 eps - sub-eps-x, classifier clean, zero fence violations).

C1 AT DENSITY, folding the FULL contact set (density lane):
- Partition BITWISE order-invariant under permutation: GT7081 same 16,889 clusters
  (countDiff 0); selfA/selfB/Havoc/GT7863 exactly identical.
- Position order-invariant to within 0.000689 eps worst-case on GT7081 (one junction of
  16,889; exactMiss 1/16889), bitwise on the others. This CLOSES
  R-ORDER-UNMEASURED-AT-DENSITY for the grouping and near-closes it for the position -
  the near-tangent contacts are clearly within eps or clearly not, so the sub-milli-eps
  FP wobble flips no partition membership.
- Input-side COMBINATORIAL keys do NOT reach bundle scale: 687,332 distinct face-pair/
  triple keys on GT7081 where the geometry has 16,889 junctions (~40x over-count). So
  identity = GROUPING (order-safe, measured) + input-derived POSITION, NOT a pure
  input-combinatorial key (section 1.2).

C2 (once-only), Havoc: the one collapse computed ONE representative (vs ~196 a naive
per-member canonicalization would compute) - the collapse-once shape, concretely.

DENSITY MECHANISM (density lane, the R-TRIPLE + R-DENSITY closers): GT7081's 1,374,649
criticals are 1,374,634 DEGENERATE CONTACTS / 0 M1 triples / 15 skeleton; they fold to
16,889 junctions (table build 9.2s / 134MB). selfA 124,142 -> 8,546; selfB 122,324 ->
8,555; Havoc 1,431 -> 101; GT7863 4,476 -> 330. The junction count tracks the VERTEX
count, not the crossing count. The critical set itself shrinks ~12x (GT7081) / ~7x
(selfA/B) under bundle-clustering to the junction's canonical x - the stage-A thinning
subsumption.

NOTE ON LEVELS: the skeleton table lives UPSTREAM (arr.verts). The escape kill operated
DOWNSTREAM (emitted cap-plane images at ~1.7 eps). That level gap IS the point: the
identity emission discards is present, bounded, and order-invariant UPSTREAM; the
escape failed because it worked downstream without it. The skeleton does not close any
carrier and does not claim to - it pressures the joint and hands its numbers to the
risks above. WHETHER the upstream identity CLOSES a carrier is exactly Stage A0's
question, unanswered by any probe run so far.

## 7. Landing history and fence

- Grounding + skeleton evidence: lab notebook .claude/lane-reports/v4impl-stage0-*.
- Skeleton probe: env-gated, run, REVERTED (git diff src/ test/ empty; canonical test
  binary rebuilt from reverted source). Baseline fence: Overlap3.* 54 pass + 1 skip
  (Gate4c_HullMask_MustResolve), unchanged.
- Round-1 adversarial review: three lanes (.claude/lane-reports/v4plan-verify-closure,
  -gates, -density), all folded into this REVISION 2 - the A0 sufficiency probe,
  incidence in the table, hardened pins + PIN 6 + the progress ratchet, completed v3
  fates, restructured B/C, and the density reframe (degenerate contacts, R-TRIPLE
  closed). Fold notebook: .claude/lane-reports/v4plan-fold-*.
- This document is the artifact. Machinery delta from stage 0 and from this fold: ZERO
  (docs + notebooks).
- NEXT: STAGE A0 (the sufficiency probe) - commissioned as a separate lane after this
  fold lands. Its verdict (SEPARABLE / NOT-SEPARABLE, per carrier) re-scopes stages
  C/D and the progress ratchet before stage A implementation begins.
</content>
