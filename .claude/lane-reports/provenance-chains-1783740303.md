# Provenance-exact chains arc

Branch: explore/sweep-plane-3d-v3, base d07e8ca6. Solo editor.

## Goal / acceptance gate

Flip test/overlap3_test.cpp's `Corpus_Havocglass8_Recorded` and
`Corpus_GenericTwin7863_Recorded` from recorded-contract back to
MustResolve (no fatal + OracleCompare). Do the flip FIRST, confirm both
red, then build the mechanism (provenance-exact strip chains).

## The inherited diagnosis (settled, do NOT re-derive)

Strips reconstruct their cap-side subdivision GEOMETRICALLY
(ChainSplitVerts: eps-projection of arrangement verts onto the snapped
chord). On steep-track geometry junction clusters spread 10-100 eps, L
and R limits disagree about junction positions at that scale, retained
loops split hairs across micro-edges, and geometric reconstruction
cannot match the retained graph -> unpaired edges reach assembly and
fail closed (NonManifoldEmission "unresolvable sheet contact"). Three
hardening fixes landed (canonical-only caps, per-cluster resolver,
capEps=8eps floor) - correct but insufficient. Closure requires strips
to consume the retained cap graph's subdivision by PROVENANCE, not
geometric matching. Also audited: naively seeding cap-input edges with
piece-slot srcIds is INSUFFICIENT - PolySetAdd cancellations + MergeSrcId
conflicts erase provenance where coincident L/R content cancels (anti-
oriented stacking cancels by design, ubiquitous).

## Design questions the mechanism must answer

(a) Partial retention: strip spans whole chord, retained graph covers
    only parts; chain must subdivide consistently w/ retained boundary
    where it exists AND w/ coincident neighbor piece where content
    cancelled (cancelled segments -> L,R strips coincide, pair w/ each
    other, chains must agree exactly).
(b) Coincident chords (L + R anti-oriented): one arrangement sub-edge,
    two owners; both pieces need same subdivision of shared segment.
(c) Pieces vanishing at cap (extension collapses): single-vert chains
    today; must stay consistent.
(d) Determinism + fence: same input->same output; all pins green
    (Pin_OneArrangementPerCritical, Pin_PerCriticalCaps,
    Pin_StripSubdivFromCap, coplanar+touching suites).

## Starting hypotheses (to use or refute)

- The incidence pre-split in RemoveOverlaps2D computes per-input-edge
  vertex lists BEFORE the sweep (BuildIncidenceLists) - a per-input-edge
  out-channel from the engine (input edge idx -> ordered sub-edge vertex
  ids, surviving further sweep splits) may be the honest primitive.
- Winding pass forced-through splits add verts after incidence -> channel
  must capture sweep-time splits too (SweepEdge carries srcId through
  splits; a parallel per-ORIGINAL-edge id could ride the same plumbing,
  null-default like negEdges/capture).
- 2D fences (CrossSection*:Boolean2*, 111 tests) must stay green; default
  paths zero-cost.

## Fence rule (absolute)

Never weaken an existing test. Existing gates/pins/fences + full suite
stay green. Honest partial > weakened gate. Stop at clean commit if the
mechanism is too large; notebook records where + why.

## PLAN (steps; append observations as I go)

1. Read grounding: coplanar-crucible-1783698000.md Step 13; overlap3.cpp
   key fns; boolean2_sweep.cpp / boolean2.h/cpp key fns. [in progress]
2. Flip the two corpus gates to MustResolve; confirm BOTH red; capture
   exact failure mode.
3. Design the provenance out-channel; write it up; red-first pins.
4. Implement engine channel (null-default, 2D-fence-safe).
5. Restructure ComputeCap/chains/EmitStrips to consume provenance.
6. Iterate on Overlap3.* until gates green; guard determinism + pins.
7. Full suite (~4min) before each commit.
8. Doc section PROVENANCE CHAINS in SweepEmit3D.md.

Stop discipline: <=3 major design iterations, then commit honest partial
+ record wall.

---

## LOG

### Step 1 - grounding (DONE)
Read: SweepEmit3D.md (full), Step 13 diagnosis, overlap3.cpp
STRIPS/CAPS (ChainSplitVerts, ZipperEmit, EmitStrips, ComputeCap,
BuildCapEdgeSet, EmitCaps, SlabResolver), boolean2.h, boolean2.cpp
(RemoveOverlaps2D driver, MergeVerts, RemapAndCollapse,
BuildIncidenceLists, incidence pre-split), boolean2_sweep.cpp
(SweepWinding, SweepEdge, SweepPass, PolySetAdd, MergeVerticals1D,
CollectArrangement, CollectThenMeasure).

KEY LOCUS: overlap3.cpp buildChains (in ComputeCap, ~L876-890) calls
ChainSplitVerts (L604) = GEOMETRIC projection of ALL r.verts onto the
straight chord [r.verts[m0], r.verts[m1]] within capEps=8eps. This is
the geometric reconstruction the diagnosis condemns.

### Step 2 - DESIGN: provenance-exact chains by class-keyed subdivision

INSIGHT: the cap boundary is ALREADY provenance-exact - it is built
from arrangement sub-edges (OutEdgesToPolygons over the winding
measure's retained edges, which are true-arrangement sub-edges split
at every arrangement vert; OutEdgesToPolygons keeps collinear verts).
ONLY the strips reconstruct geometrically (ChainSplitVerts). So the
fix is localized: replace ChainSplitVerts with the ENGINE's actual
per-input-edge subdivision. Cap side unchanged.

MECHANISM: per-input-edge subdivision channel keyed by CLASS =
unordered merged-endpoint pair (mLo,mHi). Coincident edges (L piece +
R piece anti-oriented = SAME merged endpoints) share a class ->
IDENTICAL subdivision by construction (requirement (b) exact). Twins
(near the chord but on a DIFFERENT line) have a different class ->
excluded (kills the ChainSplitVerts twin-ambiguity). A real crossing
on THIS edge is included by provenance.

WHERE THE SUBDIV VERTS COME FROM, per class c:
  - incidence pre-split verts (vertex-on-edge, known before sweep) ->
    seed classSubdiv[c] in RemoveOverlaps2D at sub-edge build.
  - sweep-constructed crossings + block-rule forced-through points ->
    recorded during the ARRANGEMENT pass (SplitAt, ProcessEvent
    forced-through) against the SweepEdge's class.
Cancelled edges (L+R annihilate at seed): get only incidence verts;
correct because an absent segment gets no sweep crossings (a proper
crossing needs both edges present) and its coincident partner shares
the same incidence set - self-consistent, verified by reasoning.

CANCELLATION HANDLING (the audited srcId-seeding failure): classId is
a LABEL in PolyVal, never affecting m/erase/conflict. On PolySetAdd/
PendingAdd merge, keep existing classId (coincident -> same class ->
no loss). Different-class SAME-oriented collinear overlap (+2, a
coplanar sub-case) would drop incoming class's future splits -> flagged
as a residual to INSTRUMENT; corpus is non-coplanar steep-track so
likely absent. Worst case = fail closed (honest), not garbage.

PLUMBING (all null-default, 2D zero-cost, arrangement geometry BITWISE
unchanged - classId is a pure label):
  - PolyVal gains int32 classId=-1 (rides like srcId).
  - SweepEdge gains int32 classId=-1.
  - SweepWinding gains const vector<int>* edgeClass (per-input-edge)
    + map<int,vector<vec2>>* classSubdiv (in/out; ADDS crossings).
  - RemoveOverlaps2D gains vector<vector<vec2>>* edgeSubdiv out;
    assigns classIds by (min,max) merged-endpoint dedup, seeds
    incidence verts, threads to SweepWinding, assembles per-INPUT-edge
    ordered chains (endpoints + interior filtered to param in (0,1),
    sorted, oriented input v0->v1).
  - ComputeCap: pass &edgeSubdiv; buildChains -> leftChains[k] =
    edgeSubdiv[k], rightChains[k] = edgeSubdiv[nL+k]. Delete
    ChainSplitVerts geometric path.

Positions are bitwise in r.verts (endpoints = merged verts; crossings
flow through getId), so strip corners == cap corners bitwise -> pair.

Determinism: std::map ordering + deterministic seed order.

### Step 2a - about to flip gates red-first (mandate 4)
Flip Corpus_Havocglass8_Recorded + Corpus_GenericTwin7863_Recorded
from recorded-contract (fatal-OR-resolve) to MustResolve (no fatal +
OracleCompare). Confirm BOTH red + capture failure mode BEFORE
building.

### Step 2a RESULT - both gates RED (confirmed)
Havocglass8_MustResolve: fatal=3 NonManifoldEmission "unresolvable
sheet contact" (244 ms - my iteration fixture).
GenericTwin7863_MustResolve: fatal=3 same (8.6 s - slow, use
Havocglass8 to iterate).
fatal=3 = the sheet splitter's unpaired-fan arm = the diagnosed
geometric-reconstruction failure. Red-first satisfied.

### Step 4 - MECHANISM LANDED, but insufficient alone (diagnosis)

Provenance chains implemented (engine channel + ComputeCap consume).
Build clean, 47/50 Overlap3 pass (1 skip Gate4c), NO regressions in
pins/coplanar/touching. Corpus gates STILL red but the FAILURE MODE
CHANGED: old = UNBALANCED fan (1F/0B, unpaired edges); new = OVERLAP
(k=4, in-plane cap edge, two same-facing adjacent). So chains FIXED
the unpaired-edge class (their designed target) and exposed a deeper
blocker.

INSTRUMENTED ROOT CAUSE (Havocglass8, cap x=13494.10429, eps~2.25e-8):
- The cap has micro-edge twins in the RAW INPUT (extended piece
  endpoints), 1.8e-7 (~8eps=capEps, borderline) to 4.7e-6 (~208eps!)
  apart. capEps=8eps does NOT merge them.
- The 208eps twin is WELD-vs-TRACK: one piece endpoint extends along
  a real 3D track (lands at the true junction at c); its twin is a
  forced-through WELD (frozen at its section position ~208eps off).
- Source: the SECTION already carries these as twins ~62eps apart
  (winding pass constructs T-junction verts without merging). The
  SlabResolver clusters section verts within 1x eps, so 62eps twins
  are NOT clustered -> resolved independently -> one finds a track,
  the other welds -> extension amplifies 62eps to 208eps.
- The cap arrangement keeps both twins (>capEps) -> micro-edge ->
  TriangulateCap makes ~zero-area SLIVER cap triangles (third vert ==
  the micro-edge endpoint, >eps so not weld-dropped) -> SplitTouching
  Sheets sees k=4 with two same-facing adjacent -> OVERLAP fatal.

This is EXACTLY the diagnosed "irreducible remainder": L/R limits
disagree about a junction by 10-100eps (here 208eps), amplification
unbounded, no constant radius fixes it. The provenance CHAINS (the
briefed mechanism) address SUBDIVISION agreement, not endpoint
POSITION divergence. The chains are correct + necessary (unpaired
class gone) but the twin-position divergence is a SEPARATE mechanism
(the SlabResolver/extension, coplanar mechanism-3 3D-IDENTITY, which
is NOT implemented in Extend - it interpolates, never snaps to a
terminating/incident 3D vertex).

### Step 5 - decision: attempt provenance-based cap-input unification
The chains alone can't close the gate; the twins must be unified by
PROVENANCE (same 3D junction -> same cap-input position). Checking
what 3D-vertex data is available at a critical to snap endpoints to
3D-identity. If it doesn't converge in-budget, commit the honest
partial (chains + doc + notebook), restore gates to recorded-contract
(failure is still NonManifoldEmission = a named guard, improved from
unbalanced to a documented sliver-overlap residual).

### Step 6 - BOTH gates share the twin blocker; DECISION: honest partial

GenericTwin7863 also fails on a twin micro-edge in a cap plane (dx=0):
UNBAL 1F/0B at A,B ~2.8e-6 apart (~62eps at scale ~31165). Its
manifestation differs (an unpaired cap edge whose micro-edge no strip
covers - two DIFFERENT pieces' twin endpoints joined by a cap edge)
but the ROOT is identical: L/R extensions of one junction land 60-200
eps apart, cap keeps both as a micro-edge, and either slivers
(Havocglass8) or a strip-less cap edge (GenericTwin) results.

WHY CHAINS CANNOT CLOSE IT (the precise wall): design question (b)
assumes coincident L/R pieces map to "one arrangement sub-edge, two
owners." That premise REQUIRES the endpoints to merge. They don't
(twins 60-200 eps > capEps=8eps), so L/R become TWO sub-edges with a
micro-edge between - and no provenance SUBDIVISION scheme repairs a
POSITION divergence. Closing needs the cap INPUT twins UNIFIED by
provenance (same 3D junction -> same position), i.e. either a
topological/provenance junction merge or a true 3D-IDENTITY extension
(snap terminating/incident endpoints to the 3D arrangement vertex).
The prior lane (Step 13 layer 3) already established a CONSTANT radius
does NOT close the class (capEps is a floor, not a bound; amplification
unbounded) - so bumping the merge is a known dead end. A safe
provenance junction-merge is research-grade (the diagnosis's named
"irreducible remainder", companion of the hulls' near-coplanar arc).

DECISION (stop at 1 major design iteration, per discipline): commit
the HONEST PARTIAL. The provenance chains are correct, 2D-fence-safe,
null-defaulted, and NECESSARY (they fixed the unpaired-edge class that
was Havocglass8's first-hit failure; both gates are now blocked SOLELY
by the isolated twin-position divergence). Keep them as the mechanism
the next arc builds on. Restore the two corpus gates to their pre-arc
recorded-contract (NOT a weakening - the pre-existing honest contract;
failure is still NonManifoldEmission = a named guard), with the comment
updated to name the twin residual. Doc section + this notebook record
the wall precisely. Remove all TEMP DEBUG; full suite green before
commit.

### Step 7 - FINAL STATE (committed honest partial)

LANDED: the provenance channel (edgeSubdiv), 2D-fence-safe, null-
defaulted, arrangement geometry bitwise unchanged (classId pure label).
- boolean2.h/.cpp: RemoveOverlaps2D edgeSubdiv out-param; classId
  assignment by (min,max) merged-endpoint dedup; incidence seeding;
  AssembleChain; per-input-edge assembly.
- boolean2_sweep.cpp: PolyVal/SweepEdge classId; MergeVerticals1D
  preserves classId (ResolveActiveId helper); SplitAt/ProcessEvent
  record crossings + forced-through; SweepWinding edgeClass/classSubdiv
  threading.
- overlap3.cpp: ComputeCap consumes edgeSubdiv directly; ChainSplitVerts
  DELETED.
- test/boolean2_test.cpp: Pin EdgeSubdivProvenanceChannel (red-first,
  mutation-verified).
- docs/SweepEmit3D.md: PROVENANCE CHAINS section (mechanism +
  adjudications + the twin-divergence wall).
- test/overlap3_test.cpp: corpus gates RESTORED to recorded-contract,
  comment updated to name the twin-position residual.

RESULTS:
- Full suite 597 = 596 pass + 1 skip (Gate4c). +1 vs baseline = the
  new pin. No regressions. 2D fences (Boolean2/CrossSection) 112 green.
- Corpus gates: still fail-closed (NonManifoldEmission) BUT the failure
  mode is now the isolated twin-position divergence, not mis-
  subdivision. Havocglass8: chains eliminated its unpaired-edge
  first-failure (unbalanced -> sliver-overlap from a weld-vs-track
  micro-edge). GenericTwin7863: strip-less cap edge from the same
  micro-edge class.

WALL (precise, for the next arc): closing the gates needs a NON-
constant-radius provenance junction unification = a true 3D-IDENTITY
extension (snap L and R endpoints of ONE 3D junction to the same
point; the coplanar mechanism-3 preference, unimplemented in
SlabResolver::Extend, which interpolates and never snaps). A constant
merge radius is a known dead end (Step 13 established capEps is a
floor). This is the research-grade remainder, companion of the hulls'
near-coplanar arc; the provenance chains are its necessary, landed
foundation. Stopped at 1 major design iteration per discipline.

### (earlier) Step 3 - engine channel
Order: (1) PolyVal+SweepEdge classId label; (2) SweepWinding channel
params + recording in arrangement pass; (3) RemoveOverlaps2D classId
assignment + edgeSubdiv assembly; (4) ComputeCap consume. Build +
Overlap3.* after each. TEMP-DEBUG instrumentation to count
different-class merges (the flagged residual) - env-gated, removed
before commit.
