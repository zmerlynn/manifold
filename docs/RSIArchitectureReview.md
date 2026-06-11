# RSI architecture review: the from-scratch design

Working document for an architectural review pass. Part 1 answers: if
`RemoveSelfIntersections` were built from scratch TODAY as a first-class
manifold feature - by someone fluent in the house architecture, with the
algorithm already proven (the pipeline in `docs/OverlapRemoval.md` is taken
as given and correct) - what would its API and internal architecture look
like? Part 2 reviews the actual branch against that design.

The algorithm itself is NOT in scope: steps 1-13, the gates, the eps
contract, and the known limitations are the proven content this review
holds constant. In scope: API surface, type boundaries, integration with
manifold's internal machinery, reuse vs duplication, threading/cancellation
posture, and file/test organization.

## Part 1: the from-scratch design

### 1.1 House integration patterns (the ground truth)

How manifold features actually integrate, from the tree:

- **Members are thin; Impl does the work.** `Manifold::Simplify`:
  `GetCsgLeafNode().GetImpl()` -> status check (`PropagateStatus` on
  errored input) -> copy the Impl -> Impl-level operations
  (`SimplifyTopology`, `SortGeometry`) -> wrap via the private
  `Manifold(std::shared_ptr<Impl>)` ctor. No serialization, no
  public-API re-entry.
- **Heavy algorithms are Impl-to-Impl.** `Boolean3(const Impl&, const
  Impl&, OpType, ExecutionContext::Impl* ctx)` -> `Result()` -> `Impl`;
  `csg_tree` wraps results with `ImplToLeaf(Impl&&)`. Not a categorical
  house ban on internal `Manifold` use (Minkowski composes through
  public wrappers internally) - the reason RSI is Impl-to-Impl is that
  it rewrites geometry directly and needs STRUCTURAL fallback identity,
  which a wrapper-typed seam cannot give.
- **Cancellation is threaded, and marks status.** ctx-aware members load
  the context with `std::atomic_load(&ctx_)` and pass it to both
  `GetCsgLeafNode(ctx.get())` (which initializes progress counters and
  threads cancel into CSG evaluation) and the Impl-level operation.
  Inside Impl code the idiom is `ADVANCE_PHASE_OR_RETURN(ctx)`:
  cancelled -> `MakeEmpty(Error::Cancelled)` - the result CARRIES the
  cancellation as a status, it is not a silent fallback. Even the
  `MeshGL` ctor sweep is ctx-aware this way.
- **`GetCsgLeafNode` collapses in place.** It replaces `pNode_` with the
  evaluated leaf under the node mutex, so after a member's first line,
  `*this` already wraps the evaluated geometry - returning `*this` later
  is both wrapper identity and geometric identity.
- **Diagnostics live on Impl.** `Impl::IsSelfIntersecting()` (properties.cpp)
  already does BVH broad phase + tri-tri narrow phase over `vertPos_` /
  `halfedge_` directly, using `GetFaceBoxMorton` + the house collider.
- **Properties are first-class.** `boolean_result.cpp` interpolates
  per-vert properties through cuts (barycentric on the source face). New
  geometry-rewriting features are expected to either preserve properties or
  state precisely why not.
- **MeshGL/MeshGL64 is the SERIALIZATION boundary** - the public I/O
  format. Internal passes do not round-trip through it; the public
  `Manifold(MeshGL64)` ctor performs construction-time processing (merge
  hints, zero-volume stripping, fresh originalID) that internal code should
  invoke deliberately or not at all.
- **Tests may include `../src/*.h`** (impl.h, utils.h, execution_impl.h
  precedents) - white-box stage seams in an internal header are
  house-conformant.

### 1.2 Public API

```cpp
/** Removes geometric self-intersection pierces... (contract as today:
 * bit-identical early-exit and fallbacks, pierce-monotonicity,
 * fail-closed gates, winding-faithful welds).
 */
Manifold RemoveSelfIntersections() const;
```

- No parameters in v1. The pipeline's scalar is a WORKING EPSILON - an
  event-identity / coincidence scale - not a Simplify-style "you may
  deviate by this much" budget; surfacing it as `tolerance` would be a
  category error, and surfacing it at all locks a still-settling
  semantic into every binding surface on day one. White-box tests that
  need an explicit eps use the internal seam (house-conformant; same as
  every other stage seam).
- If a lever is ever surfaced, it is `double epsilon = 0` (0 derives
  AlphaBudgetEpsilon from the bbox), documented together with the output
  tolerance floor it implies. That is a deliberate later decision, not a
  v1 default.
- No other public surface. The diagnostic (pierce count/depth) could later
  become public alongside `Status()`-style introspection, but that is a
  separate decision; v1 keeps it internal.
- Bindings expose the same zero-parameter form across every hand-wired
  surface this tree maintains: C (manifoldc.h/.cpp), Python
  (manifold3d.cpp), WASM embind + JS wrapper + TS types, plus the CMake
  source list. No install/export change beyond the existing library
  target; the seam headers stay private.

### 1.3 The member and the internal seam

```cpp
// manifold.cpp
Manifold Manifold::RemoveSelfIntersections() const {
  auto ctx = std::atomic_load(&ctx_);
  auto leafImpl = GetCsgLeafNode(ctx.get()).GetImpl();
  if (leafImpl->status_ != Error::NoError)
    return PropagateStatus(leafImpl->status_);
  std::optional<Impl> result =
      overlap_removal::RemoveOverlaps(*leafImpl, /*eps=*/0.0, ctx.get());
  if (!result) return *this;  // every FALLBACK arm, bit-identical
  return Manifold(std::make_shared<Impl>(std::move(*result)));
}
```

```cpp
// src/overlap_removal/overlap_removal.h (internal, not installed)
namespace overlap_removal {
// The pipeline. Returns nullopt on every FALLBACK arm (early-exit
// included): the CALLER owns "return the input bit-identically", which
// is a Manifold-identity concern, not a pipeline concern. Cancellation
// is NOT a fallback: it returns an Impl made empty with
// Error::Cancelled (the ADVANCE_PHASE_OR_RETURN idiom), which the
// member wraps like any errored result. eps <= 0 derives
// AlphaBudgetEpsilon from the impl's bbox.
std::optional<Manifold::Impl> RemoveOverlaps(const Manifold::Impl& in,
                                             double eps,
                                             ExecutionContext::Impl* ctx);
}
```

Consequences, each removing a whole class of today's residue:

- **No `Manifold` anywhere inside the pipeline.** No friend declaration,
  no `GetCsgLeafNode` plumbing, no capture-less-lambda access trick. The
  member has the access; the pipeline takes an `Impl`.
- **Fallback identity is structural.** `nullopt` -> `return *this`. The
  bit-identity contract cannot be violated by a pipeline arm constructing
  something almost-identical; there is nothing to construct. And because
  `GetCsgLeafNode` collapsed the tree in place on the first line,
  `return *this` is safe even for lazy CSG inputs - the wrapper already
  holds the evaluated leaf.
- **The three outcomes are distinct types, not conventions.** Rebuilt
  `Impl` = success; `nullopt` = fallback (-> input, bit-identical);
  `Impl` with `status_ == Error::Cancelled` = cancelled (-> empty with
  status, observable, sticky - matching Boolean3 and the MeshGL ctor).
  Conflating cancellation with fallback would make a cancelled call
  indistinguishable from "nothing to do", which the house already
  decided against.
- **Status propagation is the member's first line** - identical to every
  other member; nothing to test beyond the observable contract (which is
  already the case today, but here it is self-evident).

### 1.4 Step 1 (eps-merge) as an Impl operation

The pre-migration branch merged by writing MeshGL64 merge hints and
reconstructing a `Manifold` (superseded by D1's landed Impl-level
merge). From scratch: a private pipeline pass over a mutable `Impl`
copy -

```cpp
// Cluster verts within eps (BVH + union-find, centroid positions),
// rewrite tri verts to cluster representatives, then run the FULL
// invariant sweep the MeshGL ctor would otherwise provide:
// CleanupTopology() (pinched verts), SetNormalsAndCoplanar(),
// RemoveDegenerates(), RemoveUnreferencedVerts() (which only MARKS
// unreferenced positions NaN), and SortGeometry() (which rebuilds
// collider_ and bBox_ and compacts) - dropping a stage silently
// violates invariants downstream consumers assume. Positions-only
// intermediate: prop policy is an explicit drop (DedupePropVerts
// no-ops at NumProp()==0). Returns the max applied displacement
// (the tolerance-claim term).
double MergeVertsEps(Manifold::Impl& impl, double eps);
```

- No serialization, no `Manifold(mesh)` construction mid-pipeline, no
  exposure to construction-time semantics (the "Subtract-derived inputs
  lose run-transform info" lossiness today's code documents and dodges
  with a mergedCount==0 special case - that special case dissolves).
- This is NOT a free lunch, and the spec must say so: the merge-hint
  path delegates collapsed-tri removal, pinched-vert splitting, and
  unreferenced-vert cleanup to the well-tested ctor sweep. A direct-Impl
  merge OWNS that sweep - skip a stage and the result silently violates
  invariants every downstream consumer assumes. The hint-based path is
  therefore an architecturally defensible v1; the Impl pass is the right
  end-state because it removes the serialization round-trip and the
  lossiness dodge, not because it is simpler. The collapse machinery
  overlaps `edge_op.cpp`'s short-edge collapse and should share its
  helpers where practical.

### 1.5 The diagnostic as an Impl sibling

`Impl::IsSelfIntersecting()` already exists (bool, conservative,
distance-based via `DistanceTriangleTriangleSquared`, epsilon-distance
shared-vert skip). The pipeline's gate and the tests need a DIFFERENT
predicate - count + max depth, strict Moller-Trumbore interior pierces,
index-based shared-vert skip. The predicate difference is load-bearing
(different false-positive/negative profiles; the pierce counts are
pinned), so they stay two functions. From scratch they are one FAMILY in
properties.cpp:

```cpp
struct SelfIntersectionInfo { int interiorPierces; double maxDepth; };
SelfIntersectionInfo Impl::SelfIntersections(double relTol) const;
// IsSelfIntersecting() stays as-is; both predicates documented
// side-by-side, divergence explained at the definitions.
```

- Both reuse the stored `collider_` as the broad-phase INDEX while
  computing fresh `GetFaceBoxMorton` query boxes per call (the house
  diagnostic's exact shape) - no parallel BVH construction, no
  `GetMeshGL64` read. Valid only on finalized geometry: `collider_` is
  built by `SortGeometry`, so the diagnostic runs after the sweep, never
  mid-rebuild.
- The white-box test helper calls it through an Impl, not a Manifold.

### 1.6 Output construction

The emit builds the result `Impl` directly, running the same invariant
sweep the MeshGL ctor runs (there is no `Impl::Finish()`; the sweep IS
the contract, stated call by call):

```cpp
Manifold::Impl out;
out.vertPos_ = ...;              // ring positions
out.CreateHalfedges(tris);       // triangulated kept cycles
if (!out.IsManifold()) return std::nullopt;  // fail closed, as the
                                 // ctor's NotManifold arm does today
// Metadata FIRST, before any sweep stage: SetNormalsAndCoplanar and
// RemoveDegenerates read/write triRef, and SortGeometry permutes it -
// the MeshGL ctor likewise builds triRef before its sweep. The
// DERIVED posture: one fresh reserved meshID on every triRef
// (identity transform), meshRelation_.originalID = -1 - exactly what
// Manifold(MeshGL64) produces. InitializeOriginal() is the rejected
// alternative: it would flip OriginalID() from -1 to a fresh id, a
// public behavior change with no driver; callers who want a
// provenance root call AsOriginal() themselves.
const int meshID = Impl::ReserveIDs(1);
<triRef = {meshID, meshID, -1, tri}; meshIDtransform[meshID];
 originalID = -1>
out.CalculateBBox();
out.tolerance_ = <the measured formula>;  // BEFORE SetEpsilon: it
out.SetEpsilon();                // floors tolerance_ at the derived
                                 // epsilon_ and never lowers it
out.CleanupTopology();           // pinched verts, pre-normals
out.SetNormalsAndCoplanar();
out.RemoveDegenerates();
out.RemoveUnreferencedVerts();
out.SortGeometry(ctx);           // ctx-aware, like Hull/LevelSet
if (IsCancelled(ctx)) { out.MakeEmpty(Error::Cancelled); return out; }
if (!out.IsFinite()) return std::nullopt;
```

- No `Manifold(MeshGL64)` in the loop: construction-time degenerate
  stripping, merge-hint processing, and metadata policy stop being
  implicit pipeline semantics and become deliberate calls (or deliberate
  absences) - but every stage of the ctor sweep must be accounted for,
  present or deliberately dropped (DedupePropVerts: dropped, no props
  in v1).
- `epsilon_` stays construction-derived (`SetEpsilon` from the bbox);
  the WORKING eps never overwrites it. The measured displacement bound
  lands in `tolerance_` only. Conflating the two would either inflate
  epsilon_ past its bbox meaning or understate the tolerance claim.
- The final gate's pierce check runs `out.SelfIntersections()` on the
  Impl before the member ever wraps it.

### 1.7 Cancellation and parallelism posture

- `RemoveOverlaps` takes `ctx` and polls `IsCancelled(ctx)` at stage
  boundaries (entry, after the input pierce count, after merge, after
  step 6/7, after 9.5, per-face-batch in the partition loop, after the
  cell complex and classification). Cancelled -> an Impl made empty
  with `Error::Cancelled` (section 1.3) - observable, sticky, matching
  the house idiom; NOT a silent input-return, which would make
  cancellation indistinguishable from "nothing to do".
- The finalize sweep is ctx-aware too: `SortGeometry(ctx)` can return
  early on cancel with partial state, so the sweep checks
  `IsCancelled` immediately after and converts to
  `MakeEmpty(Error::Cancelled)` - the Hull/LevelSet pattern. v1 does
  NOT do phase/progress accounting (no `ADVANCE_PHASE_OR_RETURN`
  budget): RSI observes cancellation only, and the public WithContext
  doc says exactly that rather than overclaiming progress reporting.
- The per-face partition loop and the per-face chord passes are
  embarrassingly parallel in their READS (each face touches only its
  own edge/chord/on-edge lists plus read-only shared state); the
  from-scratch design makes the WRITES match: per-face output slots
  pre-sized by face id, written independently, gathered in face order
  afterwards - the house's deterministic-parallel idiom. Accumulating
  through a shared `push_back` is ruled out on day one: it serializes
  the loop and a later policy flip would be a refactor, not a flip. v1
  may still ship with a serial policy, but the loop shape is the
  parallel-safe one from the start.

### 1.8 Properties

v1 drops non-position properties - acceptable ONLY as an explicit,
documented reduction with a designed v2 path: every `MergedPolygon`
carries its source face; emit-time property interpolation is the
`boolean_result.cpp` pattern (barycentric on the source face for new
verts, pass-through for original verts). The data flow required for v2
(polygon -> source face -> barycentric) must not be designed away. The
public comment states the reduction (it does today).

### 1.9 Numeric helpers

- `AlphaBudgetEpsilon` lives in shared.h next to `MaxEpsilon` (as today) -
  it is a general Smith-bound helper, not RSI-private.
- The 1-leaf `Collider` limitation gets fixed IN `collider.h` (or at least
  reported upstream) rather than permanently wrapped; the `SortedBVH`
  convenience wrapper is acceptable as RSI-internal glue, but its
  single-leaf branch is a workaround for a host defect and should say so
  pointing at an upstream issue. Its Collider-facing interface uses
  `Vec`/`VecView` (house types); plain `std::vector` for CPU-only
  algorithm internals is fine.

### 1.10 File and test organization

- `src/overlap_removal/` package: `arrangement.cpp` (steps 1-9.5),
  `cells.cpp` (10-13 + emit), `overlap_removal.cpp` (driver + gates),
  `overlap_removal.h` (the one internal seam), `internal.h` (stage seams
  for tests). The current single file is ~3k lines against a house norm of
  files an order smaller; the seams already exist (the style lanes mapped
  them).
- Tests keep the white-box stage-seam pattern (house-conformant), plus the
  feature tests through the public member.
- `docs/OverlapRemoval.md` unchanged as the algorithm contract; this
  document is the architecture contract.

### 1.11 What the from-scratch design does NOT change

To be explicit about scope: the pipeline stages and their order, the eps
contract, every gate and its math, the determinism pins, the bit-identical
fallback semantics, the known limitations, the no-exceptions-in-release
posture, and the test fixtures all carry over unchanged. This is an
architecture-boundary redesign, not an algorithm redesign.

## Part 2: the branch reviewed against Part 1

Produced by two independent review lanes: a design-critique lane that
attacked Part 1 itself (its corrections are folded into Part 1 above;
section 2.5 records what changed), and a branch-vs-design mapping lane
that walked the branch surface against the design; Codex mirror lanes
fold in when they land. Classification: ARCHITECTURE-DEBT (design
right, branch diverges, migration warranted), DELIBERATE-V1
(defensible documented reduction), FINE (conformant or equally
house-faithful).

Sections 2.1-2.4 are the mapping lane's findings AGAINST THE
PRE-MIGRATION BRANCH, kept as the review record (past tense where the
finding has since been resolved); the EXECUTION STATUS block below is
the current state of the tree. Section 2.7 records the Codex mirror
lanes, which reviewed the POST-migration tree independently (Part 2
stripped from their copy).

EXECUTION STATUS (user decisions, then landed): D1+D2 executed (the
Impl-to-Impl seam, ctx threading, direct-Impl merge/emit, explicit
derived metadata, friend/lambda removal); D4 resolved by fixing
`Collider` single-leaf trees in collider.h itself (user call: reuse
over reimplementation) with a dedicated collider_test.cpp, retiring
SortedBVH's wrapper-side guard; the `OriginalID() == -1` and
cancellation pins exist and are mutation-witnessed; the V2 hook
comment is in place. V3 is superseded by D1's Impl-level merge.
Remaining: D3 (diagnostic family in properties.cpp sharing the stored
collider_), D5 (file split), V1 (write-shape restructure, then
optional parallel policy) - post-landing items per the
recommendation, with this document as the record.

### 2.1 Architecture debt

- **D1 (keystone, LANDED): the internal seam was Manifold-typed.**
  `RunOverlapRemoval(const Manifold&, double)` forced everything
  downstream of it: the friend declaration + capture-less-lambda
  `LeafImplFn` access trick in manifold.h/overlap_removal.cpp; MeshGL64
  as internal currency at all three sites (step-1 merge hints +
  `Manifold(mesh)` reconstruction with its documented mergedCount==0
  lossiness dodge; `CheckSelfIntersection` reading `GetMeshGL64()`;
  emit via `Manifold(outMesh)`); `InferEps(const Manifold&)` for a
  bbox the Impl holds directly; and output metadata policy riding the
  ctor's empty-runOriginalID branch implicitly instead of an explicit
  derived-posture assignment. One refactor (the Part 1 seam:
  `optional<Impl> RemoveOverlaps(const Impl&, double, ctx*)`)
  dissolved the whole set, exactly as predicted. Behavior was pinned
  by the suite through the migration (full-field fallback identity,
  rebuild determinism, tolerance formula); the mergedCount==0 dodge's
  REASON survives as the documented export-reconstruction lossiness
  class (now a test-helper caveat, not a production special case).
- **D2 (LANDED): no ExecutionContext anywhere.** The member loaded no
  ctx; the pipeline took none; no `IsCancelled` poll at any stage
  boundary - a caller attaching a context got silent non-cancellation
  despite RSI being in the heavy class with Boolean3/RefineToLength.
  Landed with D1 (the seam carries the ctx): polls at stage boundaries
  and every 256 partition faces, cancellation as the Cancelled-status
  Impl of Part 1, pinned by cancellation tests witnessed red against
  poll removal.
- **D3 (narrowed by D1): the diagnostic duplicates the broad phase.**
  Pre-migration, `CheckSelfIntersection` ALSO read `GetMeshGL64()`;
  D1 removed the serialization (it reads vertPos_/halfedge_ directly
  now). The REMAINING debt is the RSI-local `SortedBVH` rebuilt per
  call where `Impl::IsSelfIntersecting` queries the stored `collider_`
  with fresh `GetFaceBoxMorton` boxes. The PREDICATE difference is
  load-bearing and stays (strict Moller-Trumbore pierce count/depth vs
  conservative distance). Migration: M - `SelfIntersections` on Impl
  in properties.cpp as the documented sibling, sharing the stored
  collider.
- **D4 (LANDED, beyond the asked fix): SortedBVH's 1-leaf branch
  wrapped a host defect silently.** The finding asked for an upstream
  pointer; the user chose reuse over reimplementation and the defect
  was fixed in `Collider` itself (NumLeaves derivation, one-leaf
  Collisions arm, GetBoundingBox root-as-leaf, UpdateBoxes early-out),
  with collider_test.cpp pinning the degenerate forms and the
  wrapper-side guard retired.
- **D5: one ~3k-line file.** The stage seams already exist; the
  src/cross_section/ package precedent fits RSI's size. Migration: M,
  mechanical, best done after D1 settles the seam signatures.

### 2.2 Deliberate v1 reductions

- **V1: serial driver.** Fine as a v1 policy. The accumulation shape is
  the debt edge: `facePolygons` / `outTris` grow by cross-iteration
  push_back, so enabling `autoPolicy` parallelism is a restructure
  (pre-sized per-face slots + gather), not a flip. Part 1 rules that
  shape in from day one; the branch needs the restructure first.
  Migration: M, last - only after the Impl pipeline is stable and
  re-verified.
- **V2: properties dropped (positions-only output).** Documented in the
  public comment, and the v2 hook EXISTS in the data flow:
  `MergedPolygon.face` carries the first contributor's source face, and
  `EmitTopology::ring2Vert` preserves enough to resolve original-vs-new
  per ring vert at emit time (a post-cancellation polygon can span
  several source faces, so v2 interpolation resolves per vert, not per
  polygon). Gap: nothing in the code names these as the v2 hook -
  `face` reads as dead. S: a comment claiming the hook so a cleanup
  sweep cannot delete it.
- **V3 (SUPERSEDED by D1): step-1 merge via MeshGL64 merge hints.**
  Defensible v1 per Part 1's revised 1.4: it delegated collapse/cleanup
  to the well-tested ctor sweep instead of hand-owning those
  invariants. D1's Impl-level `MergeVertsEps` now owns the sweep
  explicitly (the shared BuildImplFromTris chain) and the lossiness
  dodge is gone.

### 2.3 Conformant

- **F1: zero-parameter public API.** Matches revised Part 1 (the
  working eps is an event-identity scale, not a Simplify-style budget;
  surfacing it as `tolerance` would be a category error). The two tests
  needing explicit eps use the internal seam - house white-box
  precedent; that include loses its reason to exist only if a public
  `epsilon` is ever deliberately surfaced.
- **F2: fallback arms return the input bit-identically** at every site;
  member-side status propagation is `PropagateStatus` first, identical
  to Simplify.
- **F3: AlphaBudgetEpsilon in shared.h** next to MaxEpsilon.
- **F4: overlap_removal_internal.h test seams** (not installed),
  matching the impl.h/execution_impl.h precedent.
- **F5: container choices.** std::vector internals with `Vec`/`VecView`
  at the Collider interface; std::map/set where deterministic order is
  incidentally useful - all within house practice (boolean_result uses
  std::map too).
- **F6: per-face read isolation holds** (each PartitionFace reads only
  its face's inputs plus read-only shared state) - the parallel-ready
  half of 1.7 that is already true; only the write side (V1) lags.

### 2.4 Unpinned behaviors a migration could silently regress

- `OriginalID() == -1` on a rebuilt output (the derived posture). If a
  D1 emit accidentally called `InitializeOriginal()`, OriginalID() would
  flip to a fresh nonnegative id with no test failing. PINNED with the
  migration (witnessed red against exactly that mutation).
- Cancellation semantics: nothing existed to pin until D2; D2 brought
  its own tests (pre-cancelled ctx -> Cancelled status, witnessed red
  against poll removal; uncancelled-ctx passthrough).
- The output-tolerance formula's `unified.maxMove` term has no
  end-to-end discriminating fixture (dropping it passes the suite;
  the `merged.maxMove` term IS pinned by the eps-chain weld). A
  discriminating fixture needs a controlled conditioned-snap remap
  moving a vert past the 10 eps floor - recorded TEST DEBT; note the
  conditioned band is also a documented known limitation whose
  residual sits outside the claim.
- D6 (no Progress() contribution) is only weakly pinnable through the
  public API: Progress() reads 1.0 both for zero scheduled phases and
  for completed accounting, and the phase counters are not reachable
  from tests. The uncancelled-ctx test asserts the 1.0 convention;
  full discrimination is accepted as unpinned.
- The non-entry cancellation polls and the merge/emit Cancelled-status
  arms have no deterministic pin: `ExecutionContext::Impl::cancel` is
  private with no test-reachable writer, no public accessor yields the
  Impl pointer, and every end-to-end path is masked by the entry poll.
  A pin would require new surface added solely for tests - declined;
  the entry-class pins (plain, empty, lazy-CSG inputs) plus the
  finalize sweep's reviewed conversion arm are the coverage.

### 2.5 What the design critique changed in Part 1

For provenance: the from-scratch design itself went through an
adversarial pass before the branch was mapped against it. It corrected
the member ctx idiom to the real `std::atomic_load(&ctx_)` +
`GetCsgLeafNode(ctx.get())` pattern; replaced an ambiguous
nullopt-for-everything seam with the three-outcome form (cancellation
carries `Error::Cancelled` status, matching `ADVANCE_PHASE_OR_RETURN`,
instead of masquerading as fallback); replaced a nonexistent
`Impl::Finish()` with the actual MeshGL-ctor sweep, stage by stage;
split `tolerance_` (measured bound) from `epsilon_` (bbox-derived,
never the working eps); withdrew the public `tolerance` parameter as a
category error; downgraded the direct-Impl merge from "simpler" to
"owns the invariant sweep" (making the branch's hint path a defensible
v1); corrected the diagnostic-family claim to the real collider_
lifecycle; and replaced the false "today's loop shape already
satisfies this" parallelism claim with the pre-sized-slot
prescription. The branch's no-parameter public API and its hint-based
merge were VINDICATED by this pass - two places the original draft of
the design would have steered worse than what is built.

### 2.6 Migration order

1. Immediately (S, no dependencies): the V2 hook comment; the
   `OriginalID() == -1` pin. [DONE]
2. D1, the keystone (L). Everything routes through it. [DONE]
3. D2 ctx threading (M) - the new seam takes ctx on day one of D1, so
   in practice D1+D2 land together. [DONE - landed with D1]
4. D3 diagnostic family (M). [remaining]
5. D5 file split (M) - after seam signatures settle. [remaining]
6. V1 write-shape restructure, then optional parallel policy (M, last).
   [remaining]

(D4 was resolved alongside, ahead of order, by fixing Collider in the
host class.)

Recommendation on timing vs upstreaming - followed: D1 (+D2 riding it)
landed before any upstream PR; the seam is exactly what maintainers
will review, and shipping the friend/lambda + serialization
round-trips would have invited a mandatory rework round. D3/D5/V1 are
honest post-landing items with this document as the record.

### 2.7 Codex mirror lanes (post-migration, independent)

Both lanes reviewed the migrated tree at af3bae84 with this Part 2
stripped from their copy.

The design-critique mirror (verdict: needs revision) found Part 1's
PROSE trailing the implementation in places where the code was already
right - the 1.6 snippet showed metadata after the sweep stages that
read/write triRef (BuildImplFromTris sets it first; snippet fixed),
1.4 understated the finalize chain (SortGeometry/RemoveUnreferencedVerts
semantics; fixed), 1.5 overstated collider sharing (stored collider is
the INDEX, query boxes are fresh per call, post-sort validity; fixed),
and 1.1's "internal code never holds a Manifold" was false as a house
claim (Minkowski; softened to the real rationale - structural fallback
identity). Two findings drove CODE: the finalize sweep now passes ctx
into SortGeometry and converts a mid-sort cancel to
MakeEmpty(Error::Cancelled) (the Hull/LevelSet pattern), and the
WithContext doc no longer overclaims progress reporting for RSI
(cancellation-only in v1 - the progress-budget work is recorded
below).

The branch-vs-design mirror (verdict: divergent) independently
re-derived exactly the recorded remaining-debt set - D3 (diagnostic
stays RSI-local), V1 (push_back accumulation), D5 (file split) - and
confirmed FINE on the member seam, the three-outcome optional, the
emit sweep, AlphaBudgetEpsilon placement, the Collider fix, and the
zero-parameter bindings. Its one NEW item is the progress-accounting
gap above (now D6). It argues D3/V1/D5 belong BEFORE the upstream PR,
against this document's post-landing recommendation - recorded as
dissent for the user's PR-scoping call. It also flagged that only the
C binding has a smoke test (Python/WASM exposure unpinned) -
post-landing test debt.

- **D6: ctx integration is cancellation-only.** No phase/progress
  accounting (no ADVANCE_PHASE_OR_RETURN budget); a caller watching
  Progress() during RSI sees no contribution. v1 posture: the public
  WithContext doc states cancellation-only for RSI instead of
  overclaiming. Implementing a phase budget is M, post-landing,
  after the stage list stabilizes.
