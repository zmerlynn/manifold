# Overlap Removal (RemoveSelfIntersections)

`Manifold::RemoveSelfIntersections()` returns a new manifold with no more
geometric self-intersections than the input. Boolean operations produce
topology-manifold output (every edge shared by exactly two triangles), but the
result may still contain pairs of triangles whose interiors cross. This pass
detects and removes those pierces while preserving topology manifold-ness, or
returns the input unchanged when it cannot.

It is implemented in `src/overlap_removal.{h,cpp}` (+ `src/overlap_removal_
internal.h`) with the per-vertex two-sided winding classifier in
`src/self_mesh_analysis.{h,cpp}`. The geometric SoS kernels it relies on are
Boolean3's own (`Shadow01`, `Kernel02`, `Kernel11`, `Kernel12` in
`src/boolean3.cpp`), reused through `Manifold::Impl::RayCast`.

## Goals

- Stricter than Boolean: Boolean3 guarantees topology-manifold output; this pass
  additionally drives down geometric self-intersection pierces.
- Safe by construction: never return output worse than the input. If the
  pipeline cannot improve the mesh, it returns the input (see Guarantee).
- Opt-in and side-effect-free: a pure query method; it does not change Boolean.

## Algorithm Outline

The pipeline follows Smith's robust-arrangement framework (issue #289), adapting
its 2D arrangement to a 3D triangle mesh:

1. `MergeVertsEps` - epsilon-merge near-coincident verts.
2-5. Enumerate edges, build the on-edge and in-triangle vertex lists.
6. `FindEdgeTriIntersections` - edge-pierces-triangle events (Moller-Trumbore
   narrow phase over a BVH broad phase), snapped to existing verts where close.
7. `GenerateChordEdges` - one chord edge per intersecting triangle pair.
8-10. Propagate interior verts onto chords; subdivide each triangle's halfedges
   at the on-edge verts.
11. `WalkPolygons` - assemble the sub-polygons of each split triangle by an
   angle-sorted next-around-face walk; only closed cycles are emitted.
12-13. `TriangulateAndEmit` - classify each sub-polygon keep/drop via the
   two-sided winding classifier, triangulate the kept ones, and run a
   pierce-aware surface-cap walker plus pre/post-cap pierce reducers to close
   residual boundaries and drop overlapping output triangles.

A pair-symmetric chord-enforcement step (Phases 1-3.5) reconciles the keep
decisions across both triangles of every chord so the output stays manifold.

## The Classifier

`AnalyzeSelfMesh` computes, per vertex, the winding numbers of the mesh just
above and just below the surface (probed at `v +/- eps * n(v)` and ray-cast).
A sub-polygon is kept when it lies on the boundary of the `winding >= 1` region.
Components are found by a `DisjointSets` flood-fill that breaks at pierced
("broken") halfedges, so each component needs only one ray-cast. Probe scale and
direction are shared (`ProbeMeshScale`, `kProbeRayDir`) so the per-vertex and
per-polygon probes agree.

## Guarantee

- **Pierce-monotonicity.** The returned manifold's self-intersection count (a
  tolerance-thresholded interior-pierce count) never exceeds the input's. The
  count may stay equal while the geometry is re-triangulated; a pierce-free input
  is returned as an equivalent re-meshed manifold, not bit-identical.
- **Fallback.** Two paths return the input. An *in-band gate reject* - pipeline
  output non-manifold, output volume drift over 50%, sign flip, or more pierces
  than the input - returns the input either unchanged or in its epsilon-merged
  form, whichever has no more pierces; a sign-flipped output is reverse-wound
  and re-checked before this fallback. An *internal exception* (a debug
  assertion, allocation failure, or a throwing constructor) returns the original
  input unchanged: no epsilon-merge runs on the throw path, since MergeVertsEps
  can itself assert or produce a non-manifold.

The gate enforcing this lives at the end of `RunOverlapRemovalImpl`; the pierce
metric is `CheckSelfIntersection`.

## Relationship to #289

Smith's framework (UCAM-CL-TR-766) resolves an arrangement of primitives with
Simulation-of-Simplicity (SoS) so every orientation test is decided
unambiguously. The 3D kernels - `Shadow01` (0D-1D), `Kernel02` (0D-2D),
`Kernel11` (1D-1D), `Kernel12` (1D-2D) - are the same ones Boolean3 uses for its
own intersection cascade; overlap removal reuses them through `RayCast` rather
than re-implementing them. Steps 1-13 above map onto Emmett Lalish's 13-step
#289 sketch; the winding classification is step 13.

## Validation

`test/manifold_test.cpp` (filter `Manifold.RemoveSelfIntersections*`) covers:
the API smoke and clean-input passthrough; genuinely self-intersecting fixtures
(hull-body minus hull-mask ~31 pierces -> 0; the self_intersect ovoids at 661
pierces, the dense-sliver fallback class) asserted via the white-box
`CheckSelfIntersection`; empty input; idempotence; determinism (repeated runs
produce identical output); and far-from-origin scaling. Each monotonicity test
carries an `ASSERT_GT(InteriorPierces(input), 0)` premise so it fails loudly if
its input ever stops piercing.

## Known Limitations

Two input classes fall back by design; both are limitations of the underlying
Boolean engine, not of this pass:

1. **Inflated tolerance.** An input whose `tolerance_` was inherited from an
   extreme-scale Boolean operand makes every geometric predicate treat all
   positions as coincident. The pass resets an obviously-polluted tolerance to a
   per-bbox value, but the worst cases still fall back.
2. **Dense slivers.** Near-coincident face boundaries produce many sliver
   triangles in the chord region (Boolean3's SoS emits these by construction).
   The cap walker absorbs a few but not hundreds, so these fall back unchanged.

A third limitation is internal rather than an input class. The per-triangle
halfedge graph does not fully resolve every chord: a chord with both endpoints
interior to one parent triangle leaves an isolated halfedge pair that stalls the
polygon walk, so that triangle emits a single polygon and is auto-kept whole;
and two intersections that snap to one shared vertex record no chord at all.
These cases are common, not rare - on the hull-mask fixture 26 of the 86
chord-bearing triangles auto-keep and 56 carry stalled halfedges - but the
downstream cap walker and pierce reducers close the residual boundaries the
incomplete arrangement leaves, so the fixture still cleans to zero. The cost
surfaces under precision stress: translated to 1e4 more chords fall into the
both-endpoints-interior path (42 of 55 auto-keep) and recovery is only partial
(5 residual pierces, down from ~31), still strictly monotonic. Gating on stalled
halfedges and splitting interior-only chords (a complete arrangement) is the
principled fix, tracked under #289; the monotonicity gate bounds the worst case
until then.

The surface-cap walker that does this recovery triangulates non-planar boundary
cycles by combinatorial fan/ear-clip without a planarity check, which can
over-inflate or leave slivers on pathological inputs; the volume-drift gate
catches the gross case (a best-fit-plane fill is the principled fix, tracked
under #289).
