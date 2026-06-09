# seed63: pipeline-induces-pierce case (investigated, no source fix)

> **Conclusion (May 2026):** `OverlapRemoval`'s gate already handles
> this case correctly via the pierceWorse fallback. No source-level
> fix is cheaper than what we have. Documented so future-us doesn't
> re-investigate.

## What it is

Running `--advfuzz --advseeds 200` reproduces a single case
(class15 mixed-op-chain seed63 at the default master seed 31415)
where:

- Input is `(sphere + cube) - cylinder` at kPow=30 displacement
  (bbox.Scale ~= 1.6e9, tolerance ~= 3e-3).
- Input has 3 above-tolerance pierces forming a 3-tri chain
  (tris 160, 163, 337). Max pierce magnitude 4e-3 (1.34x tolerance).
- Pipeline runs all phases successfully, output is manifold.
- Output has 4 above-tolerance pierces (= net regression of 1).
- Gate detects `siOut > siInput`, falls back to `mr.manifold`.
- User-facing behavior: gets input with original 3 pierces, no
  regression vs starting state.

Saved input lives at `seed31415_class15_seed63_pre2_fallback.obj`
when `--advfuzz` is run with `OVERLAP3D_FUZZ_SAVE=1`.

Reproduce a single inspection via the `--advload` mode:

```
OVERLAP3D_FUZZ_SAVE=1 build/extras/overlap3d_proto --advfuzz --advseeds 200
build/extras/overlap3d_proto --advload \
  /tmp/fuzz/seed31415_class15_seed63_pre2_fallback.obj
```

## Why it happens

Step 7p2 (`EmitNewVertsAndEdges`) creates 1 new vert at the geometric
intersection of an input edge-tri pierce. The new vert position is
mathematically the line-plane intersection of an edge with the
pierced tri's plane.

Step 13's triangulator then fan-triangulates the pierced tri around
the new vert. The 3 sub-tris are NOMINALLY coplanar with the parent
tri, but FP error in the new vert's position gives the sub-tris
slightly different normals. At near-tolerance pierce magnitudes,
that FP error is comparable to the pierce magnitude itself, so the
sub-tris produce new pierces of the same order against neighbor tris.

This is fundamental to the local fix architecture: any insertion of
a computed intersection point carries FP error bounded by tolerance,
and at the pierce-mag ~ tolerance boundary that error becomes
geometrically significant.

## What was tried (and why each path doesn't help)

1. **Pierce-aware cap walker rejection of placement.** The cap
   walker uses only existing cycle verts (no placement). Its
   `wouldPierce` check is already on and was not the source of the
   regression. Confirmed by toggling `OVERLAP3D_NO_PIERCE_CAP=1` and
   `OVERLAP3D_NO_PIERCE_GUARD=1` (same surface-cap output either
   way).

2. **Pierce-graph clustering as a predictor.** Built a connected-
   components analysis on the (tris=nodes, pierces=edges) graph.
   All 8 piercing cases in the 200-seed sweep have similar shapes
   (mostly 3-tri chains); shape does NOT distinguish seed63 from
   the 7 cases that fix cleanly. Committed as analysis tooling
   (commit `a59edf97`) but is not predictive of pipeline outcome.

3. **Static magnitude/tolerance threshold suppression in step 6.**
   Any K large enough to catch seed63 (mag/tol = 1.34) also catches
   class4 seed51 (1.49) which currently fixes cleanly. Trades one
   gate-handled regression for one lost fix. Net worse.

4. **Forward-checking at step 6 (per-etIsect BVH query).** The
   problem isn't sub-edges piercing other tris (those stay on the
   original edge's line and don't introduce new pierces). The
   problem is sub-tri normals around the new vert being FP-imprecise.
   Simulating that requires running the triangulator's fan-around-P
   logic per candidate, then computing sub-tri normals and pierce-
   testing neighbors. Substantial work, and the cheapest realistic
   forward-check is "run the whole pipeline and look at the output"
   which is exactly what the gate does.

5. **Smarter gate (magnitude-weighted fallback).** For seed63, input
   max pierce mag and output max pierce mag are the same (4e-3); all
   pierces sit in the same FP-precision-limited regime. No
   reasonable mag-weighted comparison distinguishes input from
   output. Slack ("allow N extra pierces") produces strictly worse
   user output.

## Scale of the issue

At 2000 seeds per class (16 classes = 30031 cases with different
RNG state), zero pipeline failures. The 200-seed sweep's seed63 is
an RNG-state artifact, not representative of a class of failures.
At 200 seeds: 1/3095 = 0.03%. At 2000 seeds: 0/30031.

The pipeline is robust at marginal mag/tol ratios in general; the
2000-seed sweep contains multiple cases at mag/tol = 1.05x, 1.12x,
1.20x that all fix cleanly. seed63 was a one-off, not a regime.

## Status

Open. No source-level fix planned. The gate's fallback is the
correct user-facing behavior. If a class of similar failures shows
up (= many cases that need the gate fallback at scale), revisit
with forward-checking at step 13's triangulator output instead of
at step 6.
