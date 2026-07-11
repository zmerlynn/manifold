# Performance campaign - overlap3 sweep-native emission

Branch explore/sweep-plane-3d-v3 at c36cd2a1. Solo, clean context.
Suite green at start: 597 tests (596 pass + 1 documented Gate4c skip).
Sole editor of this tree; untracked .claude/ subdirs unrelated, left alone.

## Mandate

Make the four TIMEOUT corpus cases TERMINATE (~120s each) with an honest
outcome (resolve+oracle OR a named fail-closed guard), then add them as
Recorded-contract corpus fixtures. Cap at ~4 landed optimizations;
simplicity beats the last 30%. MEASURE FIRST always: no optimization
lands without a before/after timing on the motivating case, recorded here.
Behavior preservation: full suite green; the 8 completing corpus cases
keep outcome classes (ideally identical volumes); prefer bitwise-identical
transforms. Parallelism OUT OF SCOPE. 2D engine (boolean2*) OUT OF SCOPE.

## The four timeout cases (single-mesh, self-overlap removal)

- Offset1.obj (2287 verts)
- openscad-nonmanifold-crash.obj (717 verts) - the DIAGNOSED case:
  ~3.6s through CANONICALIZE+SEAMS+SLABS -> 8794 slabs / 24k criticalXs
  / 528,307 retained pieces, then drowns in assembly. BuildImpl getVertIdx
  is a linear scan per vert over the growing weld list = O(emittedVerts^2),
  ~1M emitted verts here.
- self_intersectA.obj (8582 verts)
- self_intersectB.obj (8582 verts)

Oracle for single-mesh = the input if it round-trips clean as a Manifold.

## Known quadratic suspects (verify by measurement, ranked)

(a) BuildImpl weld getVertIdx linear scan O(emittedVerts^2) - expected
    primary wall. Fix: uniform hash grid, cell eps, query 3x3x3
    neighborhood, return MINIMUM-INDEX match within eps (exact greedy
    semantics, bitwise-identical output).
(b) SlabResolver ctor calls CollectGroupedMemberEdges(arr) PER SLAB;
    arr-constant -> hoist. Trivial, measure anyway.
(c) Per-cap engine runs: one RemoveOverlaps2D per canonical critical.
(d) BuildSlabs per-slab full-face scan; Canonicalize all-pairs merge;
    FindSeams all-pairs; M1/skeleton pairwise loops.
CRITICALS THINNING IS OUT OF SCOPE (under-inclusion breaks correctness).

## Plan

1. Build manifold_test, confirm suite green baseline.
2. Recreate the corpus + stage-timing probe (documented in Step-12
   journal; /tmp binaries stale). Run the four cases with per-stage
   timing to get a BEFORE profile per case.
3. Attack the measured top wall (expect weld first). Measure after each.
4. Cap at ~4 wins. Add the four cases as Recorded-contract fixtures.
5. Stop when the four terminate and the top remaining bottleneck needs
   engine/parallelism/divergent paths - record the profile honestly.

## Log

## Step 1: BEFORE profile (env-gated stage timing added, // TEMP DEBUG)

Instrumentation: T3LOG timers around Canonicalize/FindSeams (RemoveOverlaps3D)
and BuildSlabs/Resolvers/EmitCaps/EmitStrips/BuildImpl (SweepEmit). Probe:
/tmp/probe_corpus.cpp (single-mesh: ReadOBJ -> Impl -> RemoveOverlaps3D).

openscad-nonmanifold-crash.obj (716 verts, 1440 tris), timeout 150s:
  Canonicalize 0.001s faces=1435
  FindSeams    0.079s seams=3325 criticalXs=24029 numGroups=20
  BuildSlabs   3.346s nSlabs=8794 built=3492 pieces=528307
  Resolvers    2.013s
  EmitCaps     0.938s capArrangements=3493 emittedTris=4658
  EmitStrips   0.128s emittedTris=1055604   <-- over 1M triangles
  BuildImpl    HANGS (killed at 150s)

CONFIRMED suspect (a): BuildImpl getVertIdx linear scan O(emittedVerts^2).
1,055,604 emitted tris -> ~3.1M vertex lookups each scanning the growing
weld list. Everything upstream of BuildImpl completes in ~6.5s.
Second tier once (a) is fixed: BuildSlabs 3.3s + Resolvers 2.0s.

Offset1.obj (2285 verts): Canon .015 Seams .093 Slabs .361 Resolvers 2.213
  Caps .203 Strips .022 (246048 tris) BuildImpl 132.601s -> FATAL
  NonManifoldEmission. TERMINATES at 135s (over budget), BuildImpl dominates.
self_intersectA.obj (8582 verts): Canon .124 Seams .967 Slabs 7.141
  Resolvers 3.760 Caps 4.463 Strips .389 (4,156,680 tris) BuildImpl HANG.
self_intersectB.obj (8582 verts): Canon .138 Seams .923 Slabs 6.653
  Resolvers 3.139 Caps 3.637 Strips .321 (3,881,890 tris) BuildImpl HANG.

ALL FOUR: BuildImpl is the wall. Emitted-tri counts 0.25M-4.16M. Second
tier: BuildSlabs (up to 7s) + Resolvers (2-3.8s) + EmitCaps (up to 4.5s).
Downstream question after weld fix: CreateHalfedges/Is2Manifold/SortGeometry
on multi-million-tri output - measure whether it is O(n)-tractable at 4M.

## Step 2: WIN #1 - assembly weld hash grid (BuildImpl getVertIdx)

Replaced the O(emittedVerts^2) linear-scan weld with a uniform hash grid
(cell = eps, GridCell/GridCellHash, 3x3x3 neighbourhood query, MIN-INDEX
match within eps). Bitwise-identical to first-match linear scan: at any
query the vert set is identical in both, a within-eps match lies in the
27-cell neighbourhood, and the min-index tie-break reproduces first-match.
BuildImpl before -> after:
  openscad     HANG(>150s) -> 3.9s   (1.06M tris)
  Offset1      132.6s      -> 0.7s   (0.25M tris)
  self_A       HANG(>200s) -> 17.7s  (4.16M tris)
  self_B       HANG(>200s) -> 16.3s  (3.88M tris)
ALL FOUR now TERMINATE within budget. Outcome each: FATAL
NonManifoldEmission (a named guard - the honest steep-track junction-spread
residual documented in docs/SweepEmit3D.md, now REACHABLE instead of hung).

## Step 3: WIN #2 - hoist CollectGroupedMemberEdges out of the per-slab loop

Suspect (b). SlabResolver ctor recomputed CollectGroupedMemberEdges(arr)
(a set-dedup over all face-edges) once PER BUILT SLAB; it is arr-constant.
Compute once in SweepEmit, pass const ref to each ctor. Pure function of
arr -> bitwise identical. Resolvers stage before -> after:
  openscad  1.909 -> 0.856s
  Offset1   2.213 -> 0.116s   (numGroups=722, biggest win)
  self_A    3.760 -> 3.569s   (numGroups=3, little grouped work)
  self_B    3.139 -> 2.905s
self_intersect's Resolvers residual is the per-slab O(n^2) endpoint
cluster union-find (not CollectGroupedMemberEdges).

## Step 4: DEAD END - seenTris std::set -> std::unordered_set

Hypothesis: the dedup set's 4M per-node allocations were part of the 14s
weld+dedup. Measured (self_A weld+dedup): 14.392 -> 14.004s = NOISE.
The 14s is getVertIdx's 337M (12.5M x 27) hash-grid lookups, not the
dedup. Reverted (no measured benefit; simplicity wins, no custom TriKeyHash).

## Remaining profile after wins #1+#2 (honest, self_intersectA @ 34.5s total)

  BuildSlabs        7.1s  - per-slab full-face straddle scan O(builtSlabs*nFaces)
                            = 8932*17160 = 153M checks. Suspect (d). A fix
                            needs a face-x-range sweep (added structure) - a
                            real algorithmic change, not a simple swap.
  weld+dedup (Impl) 14.0s - getVertIdx 27-cell hash-grid lookups. This is the
                            simple-hash-grid FLOOR at 4M tris; the 27 cells are
                            required for exact min-index semantics (a within-eps
                            match can sit in any neighbour). No bitwise-safe
                            reduction.
  EmitCaps          4.2s  - one RemoveOverlaps2D per built-slab boundary
                            (8933); the engine is out of scope, per-cap setup
                            already ~0.5ms.
  Resolvers cluster 3.5s  - per-slab O(n^2) within-eps endpoint union-find.
                            A 2D hash grid could cut it BUT changes the cluster
                            ROOT (Resolve(pts[root]) is root-sensitive) ->
                            order-changing, load-bearing (5 emission-closure
                            crucible iterations). NOT touched without a proven
                            order-preserving redesign + evidence.
  SplitTouchingSheets 2.2s - radial per-fan-edge sort.

STOP DECISION: the four terminate (worst 34.5s << 120s). The top remaining
bottlenecks each need added structure (BuildSlabs sweep), are at the
simple-hash-grid floor (weld), are engine-scoped (caps), or are
order-changing on load-bearing code (Resolvers cluster). That matches the
mandate's stop condition. Landing 2 clean measured wins; capping here.

## Step 5: landed + close-out

Commits on explore/sweep-plane-3d-v3 (local only, NOT pushed):
  1f07891c  overlap3: hash-grid the assembly weld (O(n^2) -> O(n))   [WIN #1]
  a8e35d15  overlap3: hoist CollectGroupedMemberEdges out of per-slab loop [WIN #2]
  8bf7dda7  overlap3 test: record the four self-overlap corpus singles

All TEMP DEBUG timing instrumentation stripped before committing (verified:
0 occurrences of T3/TEMP DEBUG/OVERLAP3_TIMING in the committed src).
scripts/format.sh equivalent (clang-format) run on both changed files.

ACCEPTANCE MET: the four timeout cases now TERMINATE, each far within the
~120s budget, with an honest named guard:
  Offset1                     135.5s -> 1.6s    NonManifoldEmission
  openscad-nonmanifold-crash  HANG   -> 9.6s    NonManifoldEmission
  self_intersectA             HANG   -> 36.8s   NonManifoldEmission
  self_intersectB             HANG   -> 32.0s   NonManifoldEmission
(times = the new gtest fixtures' wall-clock.)

BEHAVIOR PRESERVATION: both wins are bitwise-identical transforms by
construction (min-index hash-grid weld == first-match scan; hoist of a pure
function of arr).  Evidence: full suite 596 pass + 1 skip at the weld-only
state AND the weld+hoist state (== baseline); 600 pass + 1 skip with the 4
new fixtures.  Overlap3 corpus pair fixtures (Havocglass8, GenericTwin7863)
still NonManifoldEmission; Offset2/3/4 still RESOLVE (probe verts/tris
unchanged: Offset2 16579/32978, Offset3 184/364, Offset4 381/758).  Zero
ORACLE-OFF anywhere.

STOP RATIONALE (2 landed wins, cap was ~4): the four terminate and the top
remaining bottlenecks each hit the mandate's stop condition -
  - weld getVertIdx 14s@4M-tris: the simple-hash-grid FLOOR; the 27-cell
    sweep is required for exact min-index semantics, no bitwise-safe cut.
  - BuildSlabs 7s: per-slab full-face straddle scan; a fix needs a
    face-x-range sweep = added structure / algorithmic change.
  - Resolvers per-slab O(n^2) endpoint cluster 3.5s: a 2D hash grid would
    change the cluster ROOT (Resolve(pts[root]) is root-sensitive) ->
    order-changing on load-bearing code (5 emission-closure iterations);
    needs a proven order-preserving redesign + evidence, out of this cap.
  - EmitCaps 4.2s: one RemoveOverlaps2D per built-slab boundary; the engine
    is out of scope.
Simplicity beats the last 30%; capping here per the mandate.
