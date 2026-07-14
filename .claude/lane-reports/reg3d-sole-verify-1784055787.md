# reg3d SOLE-IMPL landing - ADVERSARIAL VERIFICATION notebook

Branch: explore/sweep-plane-3d-v5  HEAD d0131fbd (4 atop bbde30e2, UNPUSHED).
Verifying: 1bdabfb0 (delete sweep) / 4ce59366 (re-pin corpus) / c650c47a (rename+docs) / d0131fbd (notebook).
Canonical READ-ONLY. This notebook is the only write; no commits.

## Plan
- A1 CLEAN FROM SCRATCH: fresh cmake -B vbuild Release -j8; zero orphans; overlap3_sweep.cpp gone; no dangling CMake.
- A2 BEHAVIOR FROZEN: 27 Regularize pins green + bitwise FNV vs bbde30e2; mutations red-capable (SoS-off, snap-off).
- A3 CORPUS RE-PIN HONESTY: reproduce 4 new pins; GT7081/GT7863 fail-closed WHY; Cray/Havoc/Offsets bitwise.
- A4 RENAME CONSISTENCY: RemoveOverlaps3D=resolver; no stale candidate-B/RegularizeImpl; phase map; SUPERSEDED; ASCII; halfedge.
- A5 ACCOUNTING: 82->31, 55 retirements named; net lines; doc open list.


## FINDINGS

### A1 CLEAN FROM SCRATCH - PASS (one stale comment)
- Fresh vcopy = git archive HEAD (zero working-tree residue); overlap3_sweep.cpp absent, untracked, no CMake ref (canonical clean).
- Fresh cmake -B vbuild Release + cmake --build -j8: BUILD_EXIT=0, ZERO warnings/errors (flags have no -Werror; no surprises). Binary built.
- Orphan grep (src+test) of all deleted sweep symbols + renamed-away names: ZERO code orphans.
- Shared islands survived: OutTri3D, SplitTouchingSheets, GridCell/GridCellHash, BuildImpl.
- NIT (real): src/overlap3.cpp:15-17 file header is STALE - still reads "The 3D sweep-native emission
  prototype: canonicalize -> seams -> slabs -> caps -> strips -> assembly. Slab construction lives in
  overlap3_sweep.cpp. Design: docs/SweepEmit3D.md." Describes the DELETED pipeline + references the
  DELETED file + points design at the SUPERSEDED doc. overlap3.h's header WAS rewritten; the .cpp's was missed.
  (test 92-94 sweep terms are the legitimate record-note; not an orphan.)

### A2 BEHAVIOR FROZEN - PASS (source-identity, stronger than FNV)
- Rename-normalized function-body diff parent bbde30e2 vs HEAD: RecordSeams/WindingAt/RobustWinding/
  FoldCoplanarClusters/EmitSeamedFace/EmitCleanFaces/BuildImpl/DecomposeComponents/ComposeComponents =
  BYTE-IDENTICAL. ResolveComponent/EmitComponentBoundary/SnapNearCoplanarClusters/GateComponent/
  ResolveComponentDirect differ ONLY in comments, whitespace-reflow, and fail-closed detail strings
  ("candidate B:"->"resolver:"). DISPATCH (RegularizeImpl->RemoveOverlaps3D): strict logic diff = IDENTICAL
  (63 lines each, comments/strings/ws stripped). => emitted GEOMETRY provably bitwise-identical; only
  runtime delta is fail-closed detail text (intentional, not asserted by kept tests).
- MUTATION red-capability:
  * snap-off (force SnapNearCoplanarClusters no-op): 3 NearCoplanar tests -> RED. Load-bearing, confirmed.
  * SoS tie-break: disabling the K>0 cascade AND flipping it AND a full sign-flip of SoSOrient3D were all
    NO-OPS on the resolver-level fixtures. Reason (honest): SoS is consumed via straddle (su==sv side-
    agreement) which is invariant to a global sign flip, AND the perturbation cascade isn't reached (these
    integer fixtures decide at e^0); coplanarity routes through DetectCoplanarClusters, not the pierce kernel.
  * exact-tie DETECTION off (Orient3DExactSign never returns 0): resolver fixtures still green, BUT
    Orient3DExactSign_PropertyPin -> RED. So the exact/SoS kernel IS red-tested - by its dedicated property
    pin (the "kernel tripwire"), not by the corpus fixtures. Consistent with the design's narrow-backstop note.
- vcopy reverted to clean; overlap3.cpp md5 == HEAD.

### A3 CORPUS RE-PIN HONESTY - PASS
- Clean pass-through pins assert BITWISE input==output vert multiset + comp/clean/dirty/regularized/failClosed
  counts + Decompose (Cray comp=1, Havoc comp=2, Offsets 39/45/1/1). Reproduction = the green run below.
- GT7081-composed FAIL-CLOSED mechanism: 13 comp / 11 clean / 2 dirty; the 2 dirty shells hit
  FatalReason::DirtyComponentUnresolved ("seam sub-face arrangement not exactly resolvable" - triple-point/
  degenerate/filter-uncertain residue, the RESOLVE/ENUMERATE axis). Any component fail-closed suppresses the
  whole compose (failClosed=2, no partial output). Reason recorded truthfully in the pin.
- GT7863-composed FAIL-CLOSED mechanism (DIFFERENT axis): 4 comp / 2 clean / 2 dirty (1 self-intersecting +
  1 within-component coplanar overlap that IsSelfIntersecting misses, caught by the re-scoped coplanar gate).
  The dirty pieces PASS the exact-zero SoS tie gate then fail NARROWER at EMISSION: coplanar-dominated
  boundary reduces to sub-eps/touching-sheet slivers with no representable double manifold =>
  FatalReason::NonManifoldEmission (stage-7 thin-cell, RE-GATE axis). Recorded truthfully.

### A4 RENAME CONSISTENCY - PASS (2 nits)
- Public entry: RegularizeResult RemoveOverlaps3D defined EXACTLY once (overlap3.cpp:2167); all test callers
  use it; the old sweep RemoveOverlaps3D (returned Overlap3Result) is gone.
- No stale candidate-B / RegularizeImpl / RunCandidateB / RegularizeDirtyDirect in src/test or the LIVE spec
  Regularize3D.md. Renamed functions all exist (Decompose/Gate/ResolveComponent/EmitComponentBoundary/
  ResolveComponentDirect/ComposeComponents). DISPATCH body wires DECOMPOSE->GATE->early-exit/ResolveComponent
  ->RE-GATE->ComposeComponents.
- Phase map (Regularize3D.md lines 90-105) matches code names exactly.
- SUPERSEDED headers present: SweepEmit3D.md, SweepPlane3D.md, MaintainedEmission3D.md. ExactArrangement3D.md
  = "HISTORICAL RECORD" pointing at Regularize3D.
- ASCII: 0 non-ASCII bytes across all changed src/test/docs.
- NIT: "candidate B" survives ONLY in docs/ExactArrangement3D.md (HISTORICAL RECORD, not in the changed set) -
  pre-existing, defensible as history.
- NIT: internal locals still named bRes/bImpl (from "candidate B"); overlap3.cpp:1333 uses hyphenated
  "half-edges". Cosmetic.

### A5 ACCOUNTING - PASS
- Test count parent 82 -> HEAD 31 (verified by counting TEST(Overlap3,...)).
- Real set-diff: naive shows 56 retired / 5 added, BUT DirtySingleComponent_RoutesToFailClosedStub ->
  _RoutesToResolver is a documented RENAME, netting 55 retired / 27 kept (26 same-name + 1 renamed) / 4 added.
  Reconciles the 55/27/4 claim exactly.
- All 55 mechanism-retired pins are NAMED in the in-tree record-note block (comm -23 = empty).
- Spot-check 5 retired: Gate1_EventParity (Overlap3Internals+TestHooks), EmissionAlgebra_SingleCube /
  Coplanar_SharedWall_Oracle / Touch_EdgeEdge (Overlap3Result+sweep RemoveOverlaps3D), SelfIntersectA_Recorded
  (CorpusSingleGate sweep scaffold) - ALL called the DELETED sweep API; none a live-behavior pin. siA/siB
  live coverage preserved (Regularize_SelfIntersect*/BMechanism_*, green).
- Line deltas match STAGE-4 claims exactly: overlap3.cpp 3402->2358 (-1044), overlap3.h 321->170 (-151),
  overlap3_sweep.cpp 290->0, test 3270->1849 (-1421).
- Doc open list (Regularize3D.md "Open list (honest)" + Risks R1) INTACT: openscad negative-winding residue,
  GT7863/PokedCube thin-cell slivers, weld artifact (R1 + SplitTouchingSheets), gate blind spots (non-fusion
  cross-component + within-component coplanar via DetectCoplanarClusters). Consumed into a gap census.

### FULL RUN (close) + resource caps
- Fresh vbuild binary, ulimit -v 6000000, timeout 900: 31/31 PASSED, 0 failed, 0 skipped, exit 0 (60.2s).
- GT7081 under strict ulimit -v 4000000: PASSED, exit 0 (~22s) - stays under 4GB, fail-closes as pinned.

## VERDICT: PASS. Landing is clean, behavior-frozen (source-identity), accounting honest, corpus mechanisms
## recorded truthfully. Non-blocking nits: stale overlap3.cpp:15-17 file header (references deleted
## overlap3_sweep.cpp + SweepEmit3D doc); ExactArrangement3D.md retains historical "candidate B";
## bRes/bImpl locals + one "half-edges". SoS perturbation cascade is exercised only by its property-pin
## tripwire, not the corpus fixtures (by design).
