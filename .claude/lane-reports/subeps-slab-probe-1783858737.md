# SUB-EPS SLAB PROBE - lab notebook

Branch explore/sweep-plane-3d-v3 @ b6443958. Solo, clean context. Sole editor.
Suite green at start: 602 = 601 pass + 1 documented Gate4c skip.

PROBE EXPERIMENT (measure, do not commit a design). Env-gated code change in
the tree; reverted at the end; only this notebook committed.

## The question (owner, challenging a settled assumption)

BuildSlabs skips slabs of width <= eps ("degenerate slabs"). But eps =
EpsilonFromScale(scale, 1000) is the Smith alpha-BUDGET tolerance (~1000x the FP
noise floor at scale), NOT machine precision. A 0.5-eps slab is ~500 ulps wide,
its midpoint exactly computable, its section exact geometry. What actually
breaks if we just BUILD those slabs? The dead-zone-c / chain-plane / width-guard
machinery exists to work AROUND skipped runs; if sub-eps slabs process fine,
that machinery may dissolve, and wall-A junction spreads (slope x
extension-distance) may shrink since extension distances shrink.

## THE CHANGE (experimental, env-gated: getenv OV3_BUILD_SUBEPS, "// TEMP PROBE")

Lower the slab-build gate (overlap3_sweep.cpp:125) from `width > eps` to strict
FP ordering: build iff `xLo < xMid < xHi` in double arithmetic (the true FP
floor - the midpoint must be a representable double strictly between the
bounds; ~2+ ulps wide).

## Predictions (recorded BEFORE running)

Knock-ons through the pipeline:

1. NO SKIPPED RUNS (mostly): with the strict-FP gate nearly every critical gap
   builds. For critical ci: li=ci-1 built, ri=ci built, so ci==li+1 ALWAYS ->
   every critical is PAIR-CANONICAL. Residual skipped runs only at ~1-ulp gaps
   (xMid not strictly between), which are machine-precision, not eps-scale.

2. CHAIN-PLANE RULE degenerates gracefully: EmitCaps sets chains[ci-1].hiX =
   crits[ci] = slabs[ci-1].xHi and chains[ci].loX = crits[ci] = slabs[ci].xLo.
   loX/hiX == slab bounds -> ZipperEmit spans nothing -> the rule is a no-op.

3. ComputeCap wideRun GUARD never fires: wideRun = crits[ri]-crits[li+1] > eps.
   With runs only ~1 ulp wide (if any), crits[ri]-crits[li+1] << eps -> never
   true. So the SubEpsFeature "macro cap over wide run" guard is dead code under
   the probe -> openscad + GenericTwin7863 (currently trip THIS guard) should
   stop tripping it. Whether they RESOLVE (oracle-true) or hit a different wall
   is the open question.

4. SUB-EPS STRIPS weld-collapse in assembly (balance-preserving drop, per
   P4-era analysis): a thin slab (say 0.5 eps) emits strips xLo->xHi width
   0.5 eps. If the track slope is modest, the two ends land within eps in yz ->
   corners weld -> degenerate triangle dropped in BuildImpl (v0==v1). Claim to
   verify: dropping preserves halfedge balance. If the slope is STEEP (~100),
   0.5 eps in x -> ~50 eps in yz -> corners do NOT weld -> a thin-but-real quad
   survives. So steep-track thin slabs produce real thin strips (the wall-A
   family shrinks because extension distance shrank from run-width to slab-width,
   but does not necessarily vanish).

5. CAPS at sub-eps-adjacent criticals may PARTIALLY WELD across planes: caps at
   c0, c1 (<eps apart in x) emit triangles in planes x=c0, x=c1. A cap vertex
   with the same yz at both -> within eps in 3D -> welds. If the transition is
   spurious (geometry unchanged across the thin slab), the two caps face +x/-x,
   coincident-opposite -> NOT deduped (dedup preserves orientation) -> could
   feed SplitTouchingSheets or the 2-manifold gate. OPEN: measure.

6. PERF: built-slab count grows. openscad 3492/8794 built -> ~8794 (2.5x). But
   self_A/B: only ~14% of 63k criticals are built-slab boundaries -> building
   all = ~7x built slabs + ~7x cap arrangements. Predict openscad ~2.5x
   (manageable), self_A/B possibly 120s+ TIMEOUT. Offset1 ~6x.

## Machinery that DISSOLVES if the probe is clean (to name in recommendation)

- BuildSlabs SubEpsFeature coverage guard (no unbuilt runs -> every straddling
  face is covered).
- ComputeCap wideRun / SubEpsFeature "macro cap over wide run" guard.
- CHAIN-PLANE RULE (StripChains loX/hiX, EmitCaps bind, ZipperEmit plane
  placement) - all no-ops when loX/hiX == slab bounds.
- EmitCaps non-canonical skip (ci != li+1) - never taken.
Load-bearing thing the eps gate protects (if probe is NOT clean): name it.

## Log

### Step 0: baseline + probe wired (DONE)

OFF baseline: Overlap3 filter 54 pass + 1 skip (Gate4c), 77.8s. Probe wired as
env OV3_BUILD_SUBEPS in BuildSlabs (strict-FP gate), a separate OV3_CONFLICT_DBG
dump at the EngineIdConflict site. Driver /tmp/subeps_driver.cpp: per-case
compose -> TestHooks (slab counts, respects env) -> RemoveOverlaps3D (timed) ->
outcome/vol/capArr. One subprocess per case, timeout-bounded. OFF driver numbers
reproduce the journals (Cray 1.576e116, Offset2/3/4 recorded vols, GT7863
SubEpsFeature, Offset1/self NonManifoldEmission, openscad SubEpsFeature,
hull/GT7081 EngineIdConflict, Havoc NonManifoldEmission).

### Step 1: 12-case evidence table, OFF vs ON (DONE)

built/tot = built slabs / total slabs (from TestHooks; a BuildSlabs fatal
returns 0 slabs, so built=0 means the conflict fired INSIDE BuildSlabs).

| case     | OFF class          | OFF vol      | OFF b/tot   | OFF t  | ON class          | ON vol       | ON b/tot     | ON t     | delta                    |
|----------|--------------------|--------------|-------------|--------|-------------------|--------------|--------------|----------|--------------------------|
| Cray     | RESOLVE            | 1.576079971e116 | 1/1      | 0.00s  | RESOLVE           | 1.576079971e116 | 1/1       | 0.00s    | none (1 slab, huge eps)  |
| Offset2  | RESOLVE            | 209026.719   | 583/2387    | 0.24s  | RESOLVE           | 209026.719   | 1283/2387    | 0.34s    | vol=0; 2.2x slabs 1.4x t |
| Offset3  | RESOLVE            | 10803.47269  | 15/73       | 0.00s  | RESOLVE           | 10803.47269  | 32/73        | 0.01s    | vol=0; 2.1x slabs        |
| Offset4  | RESOLVE            | 15240.58292  | 21/118      | 0.01s  | RESOLVE           | 15240.58292  | 53/118       | 0.01s    | vol=0; 2.5x slabs        |
| GT7863   | SubEpsFeature      | -            | 245/890     | 0.15s  | **EngineIdConflict** | -         | 0 (in Slabs) | 0.01s    | CLASS FLIP               |
| Offset1  | NonManifoldEmission| -            | 1598/9885   | 1.60s  | **EngineIdConflict** | -         | 0 (in Slabs) | 0.13s    | CLASS FLIP (fails fast)  |
| openscad | SubEpsFeature      | -            | 3492/8794   | 4.45s  | **EngineIdConflict** | -         | 0 (in Slabs) | 0.08s    | CLASS FLIP (fails fast)  |
| GT7081   | EngineIdConflict   | -            | 0           | 11.24s | EngineIdConflict  | -            | 0            | 16.48s   | same class; 1.5x t       |
| hull     | EngineIdConflict   | -            | 0           | 6.13s  | EngineIdConflict  | -            | 0            | 33.51s   | same class; 5.5x t       |
| Havoc    | NonManifoldEmission| -            | 153/532     | 0.07s  | NonManifoldEmission | -          | 409/532      | 0.08s    | same class; 2.7x slabs   |
| self_B   | NonManifoldEmission| -            | 8926/57394  | 32.49s | NonManifoldEmission | -          | 43652/57394  | 210.70s  | same class; 4.9x slabs 6.5x t |
| self_A   | NonManifoldEmission| -            | 8932/63089  | 36.14s | NonManifoldEmission | -          | 54067/63089  | 434.19s  | same class; 6.1x slabs 12x t |

Split: (i) 4 CLEAN cases RESOLVE both modes, bitwise-identical volumes, valid
manifold - Cray/Offset2/3/4; (ii) 3 CLASS FLIP to EngineIdConflict - GT7863,
Offset1, openscad; (iii) 2 wall-B EngineIdConflict both modes, much slower ON -
hull, GT7081; (iv) 3 wall-A NonManifoldEmission both modes, bigger/slower ON -
Havoc, self_B, self_A (self_A 434s, needed a 600s cap; NOT a class flip). ZERO
cases RESOLVED that fail-closed OFF; ZERO oracle-wrong resolves in either mode.

### Step 2: full-suite ON census (DONE)

Overlap3 filter minus self_intersect (self_A hangs, self_B 210s): 48 pass, 1
skip, 4 FAILED. By class:
- Gate4a_Wedges8_MustResolve: was GREEN -> ON fatal=2 EngineIdConflict
  "coincident sections from unrelated faces". A currently-green SYNTHETIC gate
  broken. (8 wedges converging on an axis 1e-3 offset; near the concurrence the
  faces are near-coplanar and a thin slab sections them coincidentally.)
- Pin_InPlaneSkeletonCriticals: was GREEN -> ON h.fatal EngineIdConflict (same
  string) in BuildSlabs. Two rotated cubes; near-coplanar z=0 group content.
- Pin_ChainPlaneRule_WideRunResolves: precondition VANISHES - widestInteriorRun
  drops from >eps to 1.11e-16. Not a breakage: this pin exercises the exact
  chain-plane machinery the probe DISSOLVES (no wide skipped run exists when
  sub-eps slabs build). It would retire with that machinery.
- Corpus_GenericTwin7863_Recorded: ON EngineIdConflict, not in the pair gate's
  accept-list {NonManifoldEmission, SubEpsFeature} -> wrong-guard FAIL. Same
  near-coplanar mechanism.
(Corpus_SelfIntersectA would hang the suite; self_B/Offset1/openscad pass since
CorpusSingleGate accepts ANY fatal.)

### Step 3: breakage mechanism, instrumented (DONE, 1 iteration)

OV3_CONFLICT_DBG dump at the BuildSlabs conflict site (coincident section edges
w/ different srcId + source-face dihedral + slab width). GT7863, the clean case:
```
CONFLICT slab[86] w/eps=0.000 srcId(15,16) face(15,16) dihedral=0.002051deg coincident=same
```
MECHANISM CONFIRMED: the conflicting slab is ~0 eps wide (a slab the eps gate
SKIPS). Faces 15 and 16 are near-coplanar (0.002 deg dihedral) but NOT in one
coplanar group - two shallow-angle planes separate by more than the grouping
eps test over any macroscopic span (L*sin(0.002deg) > eps for L > ~1e-3), so
they were never unioned. In the thin slab both straddle and their section
segments COINCIDE -> two srcIds on one edge -> the engine's id-conflict detector
fires. This is the KNOWN near-coplanar wall (wall B: hull/GT7081, the Gate4c
narrowing) made PERVASIVELY reachable, exactly as the eps gate previously
contained it. Gate4a/Pin_InPlaneSkeleton conflict the same way but via collinear
OVERLAP (not identical endpoints), so the crude same/rev dump missed them; the
engine's own detail string classifies them identically.

Note prediction 3 refined: the wideRun/SubEpsFeature guard does NOT "dissolve
gracefully" on GT7863/openscad - it is PREEMPTED. Those cases now fatal EARLIER
(EngineIdConflict in BuildSlabs) before reaching ComputeCap, so the guard is
simply never reached.

### Step 4: predictions scorecard (DONE)

1. No skipped runs / all pair-canonical: CONFIRMED on clean geom (Offset2 built
   583->1283; chain-plane loX/hiX == slab bounds). Residual ~1-ulp gaps only
   (Pin_ChainPlaneRule widestInteriorRun = 1.11e-16).
2. Chain-plane rule degenerates to a no-op: CONFIRMED (its pin's wide-run
   precondition can no longer be constructed).
3. wideRun guard never fires: CONFIRMED but by PREEMPTION (earlier
   EngineIdConflict), not graceful dissolution.
4. Sub-eps strips weld-collapse, balance-preserving: CONFIRMED EMPIRICALLY -
   the 4 clean cases resolve to VALID manifolds at bitwise-identical volumes
   despite 2-2.5x more thin slabs, so the extra thin strips collapse (dropped
   degenerate) or thin-weld without breaking topology or changing the result.
5. Caps partially weld across sub-eps planes: no separate defect observed on
   clean geometry (identical volumes); on near-coplanar geometry the pipeline
   dies earlier at EngineIdConflict, upstream of any cap weld.
6. Perf 2-3x or worse: CONFIRMED WORSE for dense self-intersect - self_B 4.9x
   slabs / 6.5x time, self_A 6.1x slabs / 12x time (434s). Clean cases ~2-2.5x
   slabs, ~1.4x time. openscad's full slab explosion is unmeasurable ON (it
   conflicts before building them; theoretical 8794/3492 = 2.5x). Time grows
   super-linearly in built slabs (weld + per-cap engine both scale up).

## RECOMMENDATION: (c) KEEP the eps gate (for now); its dissolution is COUPLED to the near-coplanar arc + a perf fix

The owner's underlying intuition is CORRECT where it was aimed: a sub-eps slab's
midpoint is exactly computable and its section is exact geometry. On CLEAN
(non-near-coplanar) inputs, building sub-eps slabs is CORRECTNESS-PRESERVING -
Cray and Offset2/3/4 resolve to valid manifolds at BITWISE-IDENTICAL volumes
with 2-2.5x more slabs. There is no sub-eps-slab correctness defect. Prediction
4 (balance-preserving weld-collapse) held empirically.

But the eps gate is LOAD-BEARING for a DIFFERENT, orthogonal reason, and this is
the evidence for keeping it: it contains NEAR-COPLANAR EngineIdConflict. Two
faces at a shallow dihedral (measured 0.002 deg on GT7863) are within eps of
each other only inside a thin x-band; over any macroscopic span they exceed the
coplanar-grouping eps test, so they are not grouped. The eps gate SKIPS exactly
the sub-eps slabs inside that band, so their coincident sections are never built
and never conflict. Remove the gate and the conflict surfaces the instant such a
band exists:
- 3 corpus cases FLIP from resolve/benign-guard to EngineIdConflict (GT7863,
  Offset1, openscad).
- 2 currently-GREEN synthetic gates break the same way (Gate4a_Wedges8,
  Pin_InPlaneSkeletonCriticals) - this is not just a real-mesh artifact.
This is the SAME wall B that hull/GT7081 and the Gate4c narrowing already
document; the gate keeps it contained to genuine built-slab boundaries instead
of surfacing it at every sub-eps slab.

Second load: PERF. Building all sub-eps slabs is 2-2.5x slabs/engine-runs on
clean cases (1.4x time) but 4.9x slabs/6.5x time on self_B and 6.1x slabs/12x
time (36s -> 434s) on self_A. The slab/cap/weld counts scale with built slabs,
and dense real geometry packs 84-86% of criticals into sub-eps runs.

WHAT WOULD DISSOLVE if adopted (named, for the record): the CHAIN-PLANE RULE
(StripChains loX/hiX, EmitCaps bind, ZipperEmit plane placement), the ComputeCap
wideRun/SubEpsFeature "macro cap over wide run" guard, the BuildSlabs
SubEpsFeature coverage guard, and the EmitCaps non-canonical-critical skip - all
become no-ops or dead code when there are no wide skipped runs. That is a genuine
simplification and the dead-zone-c class it fights goes away. So the prize is
real; the blocker is that adoption REQUIRES, as hard preconditions (option (b)
conditions, both currently out of scope / research-grade):
1. A LOCAL near-coplanar conflict tolerance (per-section / overlap-region-scoped
   grouping) so a thin near-coplanar section no longer fatals - the exact "next
   arc's opening problem" recorded at the COPLANAR implementation close and the
   3D-IDENTITY / near-coplanar arc. WITHOUT it, the strict-FP gate converts
   today's dead-zone-c and benign SubEpsFeature outcomes into EngineIdConflict
   on the same inputs - a lateral move, not progress.
2. A perf strategy for the sub-eps slab explosion (self_A 6x slabs / 12x time),
   e.g.
   building only sub-eps slabs that actually carry a geometry change, which
   reintroduces a width/critical judgment - i.e. a smarter gate, not no gate.

Net: keep the eps gate now. Re-run this exact probe AFTER the near-coplanar arc
lands a local conflict tolerance; at that point the strict-FP gate plus a perf
bound is the natural way to retire the dead-zone-c/chain-plane machinery, and the
4 clean-case bitwise-identical resolves say the assembly side is already ready
for it.

### Step 5: revert + verify + commit (DONE)

Reverted src/overlap3_sweep.cpp to HEAD (git restore): git diff HEAD -- src/ is
empty, zero OV3_BUILD_SUBEPS/OV3_CONFLICT_DBG/TEMP PROBE remnants. Rebuilt clean;
Overlap3 filter 54 pass + 1 skip (Gate4c), 77.3s = baseline (identical tree).
Committed only this notebook; not pushed. Probe driver /tmp/subeps_driver.cpp
left in /tmp (re-buildable per the journal command).
</content>
</invoke>
