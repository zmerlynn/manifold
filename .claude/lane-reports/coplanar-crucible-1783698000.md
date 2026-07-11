# Coplanar-family crucible - main-agent notebook

Branch explore/sweep-plane-3d-v3 from bd2a330e. User directive:
resolve coplanar FIRST (huge CAD class), ahead of the corpus
measurement.

GOAL: coplanar face overlap + edge-in-plane IN SCOPE (resolved, not
fail-closed). DONE GATE: existing 35+1 gates + fences + suite green
PLUS coplanar oracle gates; Gate4c hulls unskipped if EdgeInPlane
lands. Cap 3 rounds.

## Framing analysis (before any code)

Why coplanar is fatal today, decomposed:
(a) Shared plane PERPENDICULAR to sweep (normal ~ +-x): both faces
    are section-invisible (parallel-to-section skip); their content
    is pure CAP arithmetic - the cap at that critical computes
    region differences of adjacent limits, which IS the coplanar
    resolution. HYPOTHESIS: free once the fatal is lifted.
(b) Shared plane PARALLEL/OBLIQUE to sweep: sections cut pi into
    COLLINEAR OVERLAPPING segments from different faces. The engine
    merges/cancels these natively (PolySetAdd, MergeVerticals1D,
    Smith overlap handling); multiplicities sum (+2 same-oriented,
    0 anti-oriented stacking - both algebraically right). What
    breaks: ATTRIBUTION (MergeSrcId -> -1 -> EngineIdConflict
    fatal) and TRACKS (a merged piece's endpoint is a breakpoint
    created by ONE member's boundary edge - which member's track?).
    HYPOTHESIS: per-breakpoint attribution is determinate (the
    breakpoint IS some member's segment endpoint = that member's
    cutting edge = class-i for that member); the piece needs a
    GROUP to search, not a single face id.
(c) EdgeInPlane (T-contact, hulls): F transversal, one edge in G's
    plane. Sections get T-junctions (incidence pre-split handles);
    the F-G seam is a real segment usable as a class-ii track.
    HYPOTHESIS: the fatal is a June-PSLG-architecture carryover;
    the current architecture may mostly handle it.

RISKS -> lanes (after probes sharpen them):
R1 track-resolution rule for merged pieces (load-bearing)
R2 perpendicular-case-is-free claim
R3 EdgeInPlane relaxation vs hulls
R4 winding semantics: |m|=2 boundaries, mixed-orientation groups,
   inverted shells, the flap (same-winding double cover)
R5 engine attribution change must be null-default (2D fences)
R6 criticals completeness: in-plane edge-edge crossing x's
   (M1-analog) or slabs change combinatorially mid-slab
R7 near-coplanar (within-eps vs just-outside-eps) razor seams
R8 stage-A interplay: partially-overlapping coplanar faces sharing
   verts (the flap/touching family)

## Step 1: probes (lift fatals, observe what ACTUALLY breaks)

P1: stacked-along-x boxes (perpendicular shared plane) - expect cap
    arithmetic to just work.
P2: boxes sharing a z-plane with overlapping area (parallel case) -
    expect EngineIdConflict in slabs.
P3: Gate4c hulls with EdgeInPlane lifted - observe.
Probe method: temporary env-var bypass of the two fatals + a scratch
test binary run; no committed changes.

### Probe results (scratch driver, all fatals bypassed)

BASELINE (no bypass): identical-face stacking ALREADY WORKS dv=0
(stage-A annihilates coincident faces; caps do the rest). Every
partial-overlap case dies on EdgeInPlane FIRST (the small box's
side-face edges lie in the big face's plane) - EIP is the gating
detector for the entire stacking family, not CoplanarOverlap.

ALL BYPASSED:
- P1 stackedX smallOnBig (plane PERP to sweep): manifold, dv=0.
  Hypothesis (a) CONFIRMED - pure cap arithmetic, free.
- P2 stackedZ smallOnBig + sideBySide (plane PARALLEL to sweep,
  ANTI-oriented faces): manifold, dv=0, conflicts=0! Cancellation:
  +1/-1 collinear segments annihilate in PolySetAdd BEFORE
  MergeSrcId - no conflict is ever recorded. The dominant CAD class
  (stacking/walls) is conflict-free by cancellation.
- P2c sameTop overlap (SAME-oriented +1+1 = +2): THROWS "retained
  piece has no face track" - srcId=-1 from the merge, no track.
  THE one mechanism gap.

### Design shape (post-probe, much smaller than framed)

1. SEAMS: CoplanarOverlap fatal -> COPLANAR GROUPS (union-find by
   plane-within-eps). Coplanar pairs skip seam creation but add the
   M1-ANALOG criticals: in-plane boundary-edge crossing x's (for
   axis-aligned these coincide with vert x's - P2c reached the
   track crash with correct geometry - but rotated coplanar pairs
   need them).
2. ATTRIBUTION WITHOUT ENGINE CHANGE: seed section edges of grouped
   faces with the GROUP id (ids >= nFaces) instead of the face id.
   Coincident group segments then merge with NO MergeSrcId conflict
   (same id). Genuine non-coplanar conflicts still fatal. The
   engine is untouched.
3. TRACKS: per-ENDPOINT member selection for group-attributed
   pieces - an endpoint matching a member's section-segment end is
   class-i via THAT member's track; else seam track; else weld.
   Two endpoints may track different members; both trajectories lie
   in the shared plane, so the strip quad stays in-plane.
4. EdgeInPlane fatal dies on this evidence (T-junctions are
   incidence-handled); lanes hunt residual hazards.
5. THE FLAP re-adjudicates: same-winding shared-edge overlap
   becomes +2 content and RESOLVES (regularized single cover) -
   the correct overlap-removal semantics, superseding the
   fails-closed-downstream story.

Next probe: Gate4c hulls with EIP bypassed.

### Gate4c hull probe (all bypasses)

Dies DEEPER: "retained directed edges must form closed walks"
(boolean2.cpp OutEdgesToPolygons - the checked-extraction guard
from the M4 audit round, throwing in debug). Real geometry breaks
beyond the mechanism probes for an undiagnosed reason. Calibration:
boxes prove the mechanism shape; hulls remain the EMPIRICAL
acceptance fixture. D1's highest-value assignment.

## Step 2: design written + lanes launched (round 1 of 3)

Design = "COPLANAR (extension design, under crucible review)"
section appended to docs/SweepEmit3D.md (uncommitted; probe
bypasses also uncommitted in tree, synced to both review copies).
Four mechanisms: (1) coplanar groups via union-find at SEAMS,
CoplanarOverlap retires; (2) group-id seeding (ids >= nFaces), zero
engine change, genuine conflicts stay fatal; (3) per-endpoint
member track selection; (4) in-plane crossing criticals (M1
analog, axis-aligned fixtures mask it). EdgeInPlane dies; flap
re-adjudicates to resolve-as-+2; Gate4d evolves to
resolve-or-named-guard.

Lanes: D1 (Codex, review tree) = engine reality: R3 seeding trace,
R1 endpoint matching, R7 DIAGNOSE THE HULL OPEN WALK empirically.
D2 (Codex, review-reuse) = winding semantics: R2 |m|=2 end-to-end,
R6 mixed/inverted groups, flap semantics, Gate4d wrongly-resolving
attack. D3 (Sonnet, canonical read-only) = R4 criticals
completeness enumeration, R5 near-coplanar eps boundary +
group-chain diameter question, test-surface gaps.

## Step 3: D1/D2 verdicts (D3 pending)

D1 BREAK (+3 NC): the four mechanisms do NOT close Gate4c -
instrumented the open cap walk to x=-107.465, open path
32->33->35->45, endpoint gap 2.98e-8 ~ 85 eps (hull scale eps
~3.4e-10). Root shape: BuildCapEdgeSet extends adjacent pieces
INDEPENDENTLY; two pieces sharing a section vertex resolve tracks
independently and land apart. Its fifth-mechanism proposal:
per-critical endpoint BINDING - resolve the extension ONCE PER
SECTION VERTEX (3D incidence identity), all incident pieces reuse
it. My working hypothesis for the 85-eps magnitude: extrapolation
amplification (a track with sub-eps dx and macro dy/dz) OR
weld-vs-track cascade mismatch between neighbors - needs my own
validation probe before the fold commits to a mechanism shape.
D1-3: grouping must PRECEDE the shared-edge exemption (flap pairs
must join groups before seam-skip). D1-4: non-coplanar faces
cannot section-overlap collinearly (planes meet in a line -> point
crossing only) - I verified this geometric claim myself: residual
conflicts = missed grouping (near-coplanar), reframe as guard.

D2 NEED-CHANGE/4: (1) +2 content REGULARIZES to unit boundary
inside the engine (probe: duplicate same-id square -> out=4 not 8,
area 1) - R2 mostly evaporates, design text needs the precise
invariant; (2) strips are positive-region boundaries only
(negEdges is caps-only) - need pins for mixed-sign and inverted
coplanar groups; (3) converges with D1-2: group id cannot flow
through the one-id-one-track table - the resolver must be defined
with tie rules, fail closed when unbounded; (4) flap needs a REAL
closed-shell folded fixture (its naive one never reached the
pipeline). Gate4d attack FAILED to break: near-parallel resolved
dv=3.49e-15 vs oracle - evolution to oracle-or-named-guard stands.

CONVERGENCE: D1-2 + D2-3 (resolver underspecified) + D1-1 (hull
break) all point at ONE design revision: per-section-vertex
extension binding replacing per-piece-endpoint matching.

## Step 4: my Phase-4 validation of the BREAK (EXTLOG probe)

Instrumented ExtendPtWithSeams at the failing cap. The 85-eps gap:
- LEFT: pt f=9461 CLASS-I -> (127.196017974, 153.402142203) =
  vert 32 EXACTLY (the 3D vertex at the critical - class-i is
  exact by construction, edges terminate at criticals).
- RIGHT: pt=(127.196017997, 153.402142184) f=9347 WELD (interior,
  atP0=6.57, seam search MISSED) -> constant-frozen at xMid =
  vert 45. The junction moved 2.3e-8 between xMid and the critical
  and the weld cannot track it.
WHY the seam search missed: the junction's bounding structure on
the right is a coplanar/degenerate contact - under bypass (and
UNDER THE DRAFT DESIGN TOO) those contribute criticalXs but NO
seam/track object. The 9347x9461 seam ended AT the critical (9461
died there), so the right slab legitimately has no such track; what
bounds 9347's piece there needed an IN-PLANE track that does not
exist. D1's BREAK is real and sharper than stated: mechanism 4
gives in-plane junctions their CRITICALS but nothing to EXTEND
ALONG.

FOLD SHAPE (mechanism 5, replacing draft mechanism 3):
(a) resolve extensions PER SECTION VERTEX once (all incident
    pieces reuse - same-side closure by construction);
(b) track candidates grow a third class: in-plane boundary edges
    of coplanar-group members (ordinary 3D edge interpolations);
(c) 3D-identity preference: when the governing track terminates at
    an arr.verts vertex at the target critical, the extension IS
    that vertex (exact, not eps);
(d) weld shrinks to genuine block-rule artifacts, cap coverage
    argument unchanged.
EXTLOG instrumentation to remove before the fold commit.

## Step 5: round-2 fold + relaunch

Design section REWRITTEN in docs/SweepEmit3D.md ("round 2 - post
D1/D2 fold"): five mechanisms - groups-as-pre-pass incl.
shared-edge pairs [D1-3]; group-id seeding + unit-boundary
regularization invariant [D2-1] + strips-positive-only [D2-2] +
near-coplanar-guard reframing [D1-4]; PER-VERTEX extension
resolution with class-iii in-plane member edges + 3D-identity
preference + deterministic ties + fail-closed [D1-1/D1-2/D2-3];
in-plane crossing criticals; checked closure as the Gate4c
boundary. Flap re-adjudication + closed-shell fixture requirement
[D2-4]. Gate4d evolution recorded with D2's failed-attack evidence.

Lane housekeeping: D3 (Sonnet) DIED at the 32k final-message cap
with an EMPTY notebook log - it ignored the incremental convention;
only its plan survived. Relaunch D3b was POLICY-BLOCKED instantly
(false positive, 0 tool uses - suspect the "reuse the dead lane's
plan" phrasing). D3c relaunched with rephrased prompt + hard
output discipline (findings INTO the notebook as derived, final
message < 10 lines), attacking the round-2 design directly.
D1-round2 convergence lane launched in parallel: walk ITS OWN hull
junction through the revised mechanism 3, adjudicate class-iii
enumerability (per-slab member edges vs hidden in-plane
arrangement), hunt new holes (per-vertex table vs pair-canonical
binding; caps vs raw extensions).

## Step 6: round-2 verdicts

D1-R2: SURVIVE - all four round-1 findings RESOLVED, including
walking its own hull junction through revised mechanism 3
(class-iii -> 3D-identity -> same arr.verts vertex both sides) and
adjudicating class-iii enumerability (member-edge scan suffices, no
hidden in-plane arrangement). No new findings.

D3c: DIED at the 32k cap AGAIN (10 tool uses, 42 min) - notebook
got plan + grounding notes (progress over D3's empty log) but the
attack derivations burst in the final message and were lost. Two
deaths on the same angles -> executing D3's attacks MYSELF
(main-agent-for-failed-lanes, the confirmed pattern).

## Step 7: D3 attacks, main-agent execution

### R4 criticals completeness - enumeration

Events changing a coplanar group's section combinatorics along x
(section of group = 1D interval arrangement on the line L(x) =
plane cut, breakpoints = L(x) crossing in-plane 1D structure):
1. member vertex crosses sweep plane: vert x's, COVERED.
2. member-edge x member-edge in-plane crossing: mechanism 4,
   COVERED.
3. in-plane tangency (endpoint-on-edge): the endpoint is a vertex,
   COVERED by vert x's.
4. collinear-overlapping member edges in-plane: combinatorics
   change only at edge ENDPOINTS (verts, covered); between them
   the coincident breakpoints are one 3D line - class-iii
   interpolations agree; tie rule picks deterministically. BENIGN.
5. seam(T, member) x member-edge in-plane crossing, T transversal:
   the seam of a transversal face with a member LIES IN THE PLANE;
   where it crosses ANOTHER member's boundary edge, a section
   junction transfers (T's endpoint slides across a group
   breakpoint). NOT covered: M1 covers seam x seam sharing a face;
   mechanism 4 covers edge x edge. GAP - REAL (fixture: two
   partially-overlapping coplanar rectangles + a rotated
   transversal triangle whose seam with member 1 sweeps across
   member 2's boundary edge at a non-vertex x).
   FIX: define the group's IN-PLANE SKELETON = member boundary
   edges + in-plane seams (transversal x member); mechanism 4
   becomes pairwise skeleton-crossing x's (subsumes the edge x edge
   case; M1 remains for the non-coplanar seam x seam class).
6. group membership changing along x: membership is per-face and
   x-independent (a face is in the plane for its whole extent);
   entry/exit happens at vert x's. NOT AN EVENT CLASS. COVERED.
7. near-perpendicular planes (L(x) sweeping fast): not a
   combinatorial class; x-extent <= eps falls to SubEpsFeature;
   noted under R5 as an amplification band, not a critical gap.

### R5 near-coplanar eps boundary

(a) Separation just OUTSIDE eps (razor band, ~2-10 eps): not
    grouped; sections are parallel segments a few eps apart - they
    do NOT exact-merge (no conflict, the guard never fires) and
    instead arrange as thin slivers -> BuildImpl weld partially
    collapses them. UNTESTED BAND, honest status: add a band
    fixture to the test surface; behavior decides whether the
    grouping threshold widens or a band guard lands. NEED-CHANGE
    (test surface), not BREAK.
(b) eps-CHAINS (pairwise-within-eps spanning >> eps): the group is
    one id, but geometric coherence is delegated to the ENGINE's
    vert merge - which is ITSELF an eps-union-find over the same
    geometry, so section verts chain-weld exactly where faces
    chain-group. Consistent by construction; no diameter guard
    needed for the SECTION story (the old design's diameter guard
    served triple-point unification, a dead mechanism). PIN IT
    (3-face eps-chain fixture) and STATE the delegation in the
    design. Class-iii tracks are per-member TRUE 3D edges (no
    common-plane projection), so chained groups do not distort
    tracks.
(c) predicate symmetry: union-find membership must be the
    symmetric OR of the two direction tests (the current detector
    short-circuits after B-in-A); the grouping pre-pass must
    restructure this. PROSE/IMPL note.

### R6/R7 test surface + resolver

- coplanar x sub-eps critical pair: no new mechanism (pair-
  canonical binding + per-critical caps own it) but ADD a fixture:
  a group overlap corner within eps of a vertex plane.
- group at first/last critical: exterior-cap path already covers
  (P1 fixture ends at extremes). COVERED.
- coplanar + M1: subsumed by the skeleton rule (R4 fix).
- R7 resolver mechanics: adjudicated by D1-R2 (member-edge scan);
  candidate matching = evaluate candidate tracks at xMid, eps-match
  to the vertex, deterministic ties - same posture as today's seam
  matching. NO sixth mechanism.

D3 VERDICT (mine): NEED-CHANGE - one real R4 gap (skeleton
crossings), three test-surface additions (razor band, eps-chain,
coplanar x sub-eps pair), one predicate-symmetry note. No BREAK.

## Step 8: IMPLEMENTATION (user: "push on")

Order per brief: core gates RED first (P1/P2/P2c/Gate4d-evolution +
skeleton pin), mechanisms 1 -> 2 -> 4 -> 3, EIP retirement with 1,
then extended fixtures (three-face, mixed, inverted, chain, razor,
corner, flap), then Gate4c acceptance run.

Design notes settled while planning:
- Mechanism 4 REUSES SeamSeamCrossX verbatim (skeleton segments are
  coplanar 3D segments - exactly its contract). Skeleton = ALL
  member edges (incl. internal diagonals - over-inclusion harmless,
  finer slabs) + in-plane seams (seam with a group-member face),
  deduped by sorted vert pair.
- Mechanism 3 resolver: per-slab map from EXACT (y,z) (capture
  coords are exact section-vert coords) -> Track{vec3 a,b} | Weld.
  Candidates: face-track segment ENDS (class-i), seam yzMid
  (class-ii), group-member edge crossings of the xMid plane
  (class-iii), weld fallback. Priority = class then lowest id.
  3D-identity falls out of InterpolateSafe evaluated at a
  terminating endpoint's exact x (t=0 -> bitwise endpoint); crits
  ARE vert x's bitwise, so no separate snap mechanism is needed -
  choosing a TERMINATING track over weld is what class priority
  does. Resolver stores tracks (not precomputed points) because
  pair-canonical binding extends across unbuilt runs to foreign x.
- Mechanism 1 kills coplanarInteriorOverlap + ClipPolyByHalfplane +
  SegTriInteriorLen2D (EIP helper) - grouping is by PLANE not by
  overlap, non-overlapping same-plane faces grouping is benign.

## Step 9: implementation landed

Mechanisms 1/2/4 (FindSeams grouping pre-pass + detector deletion,
BuildSlabs group seeding, skeleton criticals) + mechanism 3
(SlabResolver: per-vertex, class i/ii/iii/weld, built per slab,
consumed by BuildCapEdgeSet; ExtendPt/ExtendPtWithSeams deleted;
per-piece attribution fatal deleted with it).

First full run: 35/40 - ALL FOUR core oracles green immediately
(incl. the +2 probe-crasher). The five reds each taught something:
- Skeleton pin: crossings ARE added but are ALSO seam endpoints
  (side faces rising from boundary edges) -> already verts. Pin
  reworked to the property form; mechanism 4 = belt for
  degenerate-adjacent cases; red-first-vs-stub not achievable on
  closed solids - recorded honestly.
- P4/P4b: edge-on-face TOUCHING = genuinely non-manifold union
  (4 faces at the welded contact line); was reaching SortGeometry's
  debug assert THROUGH BuildImpl's weaker IsManifold gate ->
  BuildImpl now gates Is2Manifold -> clean NonManifoldEmission;
  pins evolved + renamed (EdgeOnFace_Touching).
- Gate4d: RESOLVES (as designed) -> evolved to
  resolve-with-oracle-or-named-guard, oracle added.
- Gate4c: near-coplanar guard fires (shallow-crossing planes
  coincide in section locally; global grouping test misses) ->
  skip narrowed to {EngineIdConflict, NonManifoldEmission}; hulls
  remain the next arc's opening problem.
- Corner fixture accidentally built the in-run macro-change dead
  zone (M4-close (a)) - now REACHABLE with EIP retired; converted
  to recorded-contract (PerpFacesSubEpsApart_Recorded).

Extended fixtures: ThreeFaceGroup (+3), MixedOrientation, EpsChain,
RazorBand (recorded), InvertedStacking - 9/10 green on first run.
Flap closed-shell fixture DEFERRED (owed; +2 semantics covered).
Enum retirement: CoplanarOverlap + EdgeInPlane deleted; 4d/4e/4f
guard sets updated; helper deletions (clip, seg-tri-interior).

FULL SUITE: 591 + 1 skip (592 total; was 581). Overlap3 45 + 1.

## Step 10: implementation audit round + close

FIDELITY: all five mechanisms FAITHFUL with evidence (gdb-verified
Gate4d's resolve branch; verified the M4-belt claim is recorded
honestly, not quietly vacuous; Extend's miss arm unreachable -
resolver keys exactly the piece endpoints caps extend). 2 LOW
stale comments -> fixed (suite header skip claim; P4b clip-story).

SIMPLICITY: NEED-CHANGE/3 -> (1)+(2) FOLDED: shared
CollectGroupedMemberEdges (one concept, one collection - skeleton
partitions by group, resolver filters by x-span) +
AppendPairwiseCrossXs for the skeleton's pairwise loop (M1's
seam-pair loop keeps its two-line body - predicate indirection for
two differently-filtered sites is heavier than the duplication;
adjudicated). (3) skeleton pin doesn't isolate mechanism 4:
DEFENDED as recorded-owed - the doc close section already records
the degenerate-adjacent fixture as future test surface and the
fidelity lane called that recording honest.

CLOSE: design 2 rounds + impl 1 round (cap respected). Full suite
591 + 1 skip. The user's directive is delivered: the coplanar
family resolves - stacking, walls, same-oriented overlap, triple
groups, mixed orientation, near-parallel-within-eps, eps-chains,
inverted stacking - with the remaining boundaries NAMED and pinned
(near-coplanar shallow-crossing hulls at EngineIdConflict;
edge-on-face touching at NonManifoldEmission; in-run macro change
recorded; razor band recorded; flap closed-shell fixture owed).

## Step 11: touching contacts (user: "take it on now")

Boolean3 probe settled the target: union of edge-touching cubes =
16 verts, Decompose()=2 - coincident geometry, separate topology
(epsilon-valid). Our BuildImpl weld is geometric, so it fuses
touching sheets into 4-fan edges + non-disk vertex links.

MECHANISM (SplitTouchingSheets, in BuildImpl before CreateHalfedges):
1. RADIAL EDGE PAIRING: for every multi-fan edge (kF + kB
   halfedges), sort incident faces by angle around the edge axis.
   Orientation algebra (derived, then sanity-checked on P4 + the
   edge-edge cube case): a forward halfedge (lo->hi) has material
   just BELOW its angle, a backward just ABOVE -> material wedges
   alternate with empty ones -> each backward pairs with the NEXT
   forward CCW. Non-alternating pattern = overlapping material ->
   fail closed. Angle ties (tangent sheets) -> fail closed
   (pairing would be a coin flip; wrong pairing keeps volume but
   garbles topology - which is also why the tests must check
   component COUNT, not just oracle volume).
2. VERTEX SPLIT: union-find corners per vert, connected through
   PAIRED halfedges only; each component gets its own vert copy.
   Covers edge-on-face, edge-edge, vertex-only touch, and
   self-touch (C-arms) uniformly.
|F| != |B| or slivers or ties -> false -> NonManifoldEmission
(the in-run macro dead zone keeps its guard).

Tests: P4/P4b flip from pinning the fatal to pinning RESOLUTION
(oracle vs Boolean3 a+b, Decompose()==2); new edge-edge and
vertex-only cube fixtures with component-count checks.

### Step 11 outcome

All four touching fixtures green on FIRST build of the splitter
(P4, P4b, EdgeEdge cubes, VertexOnly cubes - each with
Decompose()==2 + oracle). Full gates 47+1 skip; full suite 593+1.
Discovery before implementing: vertex-only contact ALREADY passed
(shared vert, disjoint fans, no shared edge is topologically
tolerated) - the red class was exactly the 4-fan edges. The
component-count assertion is the load-bearing check: wrong radial
pairing preserves volume/SA and would pass a pure oracle test.
Doc: implementation-close note superseded; new "Touching contacts
resolved" section. P2c-perp recorded fixture unaffected (its
failure is upstream of assembly).

### Step 11 audit: SURVIVE / 0 findings

The lane re-derived the material-side convention independently
(confirms: backward-pairs-next-forward-CCW IS the non-crossing
pairing), verified against a cross-pairing counterfactual, built an
asymmetric 6-fan (three rotated corner cubes -> 3 components,
vol 3), and swept the tangent band empirically: 5e-10 rad fails in
the cap walk, 2e-9..1e-7 fails closed as unresolvable sheet
contact (FP noise breaks alternation above the tie guard - still
the correct posture), 1e-4+ resolves. Pinched-vertex single-mesh
scratch case: valid 2-component output. Touching-contacts arc
CLOSED with zero folds.

## Step 12: corpus histogram (user: "get that histogram now")

Corpus = test/models/ (the RSI-era real meshes): 4 left/right
operand pairs (Cray, 2x Generic_Twin, Havocglass8) + hull
body/mask + 4 Offset singles + openscad-nonmanifold-crash +
self_intersectA/B. Driver: one case per subprocess (crash
isolation, timeout 120s), modes pair (compose, resolve, oracle =
boolean union) and single (resolve, oracle = input if it round
trips as a clean Manifold). Outcome classes: RESOLVED-ORACLE-OK /
RESOLVED-ORACLE-OFF / RESOLVED-NO-ORACLE / FATAL(reason) / THROW /
TIMEOUT / INPUT-REJECTED.

## Step 12 results: the histogram (12 cases)

RESOLVED-ORACLE-OK: 4 - Cray PAIR (first real operand pair
end-to-end!), Offset2/3/4.
FATAL-EngineIdConflict (near-coplanar): 2 - hull pair,
Generic_Twin_7081.
FATAL-NonManifoldEmission (unbalanced fan = emission-closure
defect): 2 - Generic_Twin_7863 (1F/0B), Havocglass8 (1F/2B).
Instrumented: the splitter's unbalanced arm - the emitted
triangulation has an unpaired edge BEFORE the splitter; upstream
strip/cap closure bug on real geometry, NEW class.
TIMEOUT: 4 - Offset1, openscad, self_intersectA/B. Diagnosed via
stage probe: openscad = 717 verts -> 3.6s through slabs -> 528,307
pieces / 8794 slabs / 24k criticalXs -> assembly's quadratic weld
is the wall (~1M emitted verts, linear-scan getVertIdx). PERF
class, deferred by design; first targets known (hash the weld,
bound ChainSplitVerts, criticals thinning).
ZERO ORACLE-OFF: nothing resolved WRONG anywhere in the corpus.

Next-arc priorities from the data: (1) minimize + fix the
emission-closure defect (Generic_Twin_7863 is the smaller repro),
(2) hulls' near-coplanar guard, (3) perf campaign (weld hash
first). The 900s openscad run left going for slow-vs-hung
confirmation; result to check next session.

## Step 13: emission-closure crucible (user directive)

GOAL: the unbalanced-fan class (1F/0B Generic_Twin_7863, 1F/2B
Havocglass8) - root cause + fix. Suspects, ranked before evidence:
(a) BuildImpl exact-duplicate drop removing a NEEDED copy (odd
multiplicities on dense real geometry - the R3 reasoning assumed
exactly-two); (b) TriangulateCap silently under-covering a
degenerate/self-touching cap loop; (c) chain/zipper mismatch
(strip edge subdivided differently from its cap partner); (d)
PushSimpleLoops splitting a figure-eight and dropping an edge.
NOTE: index-degenerate triangle drops are balance-PRESERVING
(their two live edges cancel internally) - not a suspect.
Method: TEMP provenance instrumentation in BuildImpl (dump
unpaired edges with positions + which mechanism emitted the
incident triangles), then minimize (7863_left is 11KB - tiny).

### Step 13 close (5 instrumented iterations)

Layer-by-layer: (1) doubled in-run caps -> canonical-only emission
[KEPT]; (2) sub-eps section twins resolving to different tracks ->
per-cluster resolution [KEPT]; (3) capEps=8 noise floor [KEPT,
scoped honestly - floor not bound]; (4) triangulator-eps experiment
[REVERTED - the mismatch theory was wrong; the polygon itself
carried the hair]; (5) final diagnosis: L/R limits disagree about
junction position by 10-100 eps (retained loops split hairs across
a 10-eps micro-edge; "equal limit locally" fails at sub-cluster
scale). Named requirement: provenance-exact strip/cap closure
(per-input-edge retained chains + partial-retention semantics).
Corpus gates -> Recorded contracts; Pin_OneArrangementPerCritical
updated to the canonical rule. Full suite 595+1 skip (596 total).

### Step 13 audit: NEED-CHANGE/1 (prose) -> folded -> CLOSED

The lane's attacks on all three landed changes SURVIVED (cluster
root determinism verified against DisjointSets' tie rules; no
mixed-scale consumer; exterior runs resolve; right-end run variant
fails closed as recorded; Pin_PerCriticalCaps oracle-true; fences
111/111). It INDEPENDENTLY confirmed the no-bounded-provenance-fix
read: srcId seeding dies at cancellations - per-input-edge chains
need partial-retention semantics = the design arc. The one finding:
stale R3 doc text + pin banner still claimed every-critical caps
with identity dedup - folded (R3 AMENDED in doc, pin banner
rewritten). Emission-closure crucible CLOSED: 2 structural fixes +
1 scoped floor landed, remainder diagnosed/named/recorded, corpus
fixtures as recorded contracts.
