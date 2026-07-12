# WALL-B EngineIdConflict demotion - adversarial verification lane

Artifact: landing f0dab0fb on explore/sweep-plane-3d-v3. Claim under review: the
2D-engine source-id conflict guarded ONLY an attribution label (srcId) that no
3D consumer reads, so demoting FatalReason::EngineIdConflict to a benign counter
is sound. Default skeptical. Canonical repo READ-ONLY; build in
/tmp/wallb-verify-adjud (fresh cmake, Release, -j4). GT7081 never run (OOM).

VERDICT: SURVIVE.

## Attack 1 - unconsumed-label argument, adversarially (traced every channel)

srcId/sourceId write sites (grep, whole tree minus boolean2.h):
- overlap3_sweep.cpp:170-171 - the ONLY 3D writer (group-id seeding into section
  EdgeM.srcId). No other 3D writer.
- boolean2_sweep.cpp - internal engine plumbing (PolyVal.srcId, MergeSrcId,
  EmitBoundary capture, materialize).

The -1 sentinel (MergeSrcId, boolean2_sweep.cpp:132-140) flows out of the
BuildSlabs SweepWinding call into exactly two carriers:
- SweepCapture.sourceId (slab.pieces), via EmitBoundary:392.
- The discarded OutEdge return of materialize:648 (BuildSlabs ignores the return
  value entirely - only &slab.pieces and &conflictCount are taken).

Every READ of slab.pieces / SweepCapture.sourceId in the 3D pipeline:
- overlap3.cpp:545-551 (SlabResolver ctor) - reads piece.from/piece.to ONLY
  (positions; clusters + resolves tracks geometrically).
- overlap3.cpp:777-780 (BuildCapEdgeSet) - reads piece.from/piece.to ONLY, then
  RE-SEEDS cap edges with a fresh EdgeM {i0,i1,+1|-1} => srcId defaults to 0. The
  conflicted -1 is DROPPED here; it never propagates into the cap arrangement.
- overlap3.cpp:681 (EmitStrips) - reads .size() only.
Conclusion (1a): srcId=-1 never reaches an array index or map lookup. It is
inert. No consumer assumes a valid face/group id.

(1b) CONFIRMED and it is a consistency argument FOR the demotion: the cap-side
RemoveOverlaps2D (overlap3.cpp:844) has NO conflictCount parameter at all, and
internally (boolean2.cpp:954) calls SweepWinding(..., /*conflictCount=*/nullptr,
...). Caps have ALWAYS swallowed source-id coincidence silently. Moreover every
cap input edge carries srcId=0 (overlap3.cpp:828,830), so MergeSrcId(0,0) can
never even flag a conflict on the cap side. The caps already tolerated exactly
what BuildSlabs used to fail closed on.

(1c) No telemetry/debug path branches on srcId or engineIdConflicts. grep of
src/ + test/ (minus boolean2): engineIdConflicts is written once
(overlap3_sweep.cpp:207), copied into result.counters, and read by NO code and
NO test (only counters.capArrangements is asserted, test:1276). RemoveOverlaps3D
is not wired into the main boolean path (overlap3.h:17), so no user sees it.

(1d) negEdges second measure: BuildSlabs passes negEdges=nullptr, so the -1
never reaches a negated measure. negEdges is exercised only by the cap
(boolean2.cpp:954) where all srcId=0; and negEdges' OutEdge.srcId is consumed by
OutEdgesToPolygons (boolean2.cpp:159), which reads v0/v1 + positions only, never
srcId. Channel closed.

Attack 1: no hole. The unconsumed-label argument holds under a full channel
trace; the emitted boundary is srcId-independent (from/to set by net
multiplicity sign + vertex positions in materialize/EmitBoundary).

## Attack 2 - flap backstop, adversarially (EMPIRICAL, fresh build)

Built fixtures in /tmp/wallb-verify-adjud, swept the dihedral a few decades.

(i) The exact wall-B shape on the defended point (Coplanar_SameOriented tilted
by theta about x, theta in 1e-5..0.5 deg): conflicts=0 at every theta. The
near-coplanar overlap resolves through ordinary TriTriSeam machinery to the
Boolean oracle EXACTLY (dVol ~ 0, is2mf=1, genus 0/0, comps=1). At sub-degenerate
theta (1e-5,1e-4) it fails closed on NonManifoldEmission "unresolvable sheet
contact" (SplitTouchingSheets, overlap3.cpp:1165 - a zero-volume sheet, honest).
Stacked-plate variant: conflicts=0 at every theta, oracle-exact resolve.
=> The imagined near-coplanar flap does NOT reach the demoted guard; it is
handled cleanly by seams. This matches the notebook's flap re-adjudication.

(ii) Real conflict carriers + a near-degenerate battery (ConflictCarrierHunt):
- hull (Gate4c body+mask): conflicts=47 -> FATAL SubEpsFeature "macro cap content
  over a skipped run wider than eps" (the CHAIN-PLANE wide-run guard = the
  documented residual). NOT garbage. Confirms notebook step 3/5.
- GT7863: conflicts=0 -> FATAL SubEpsFeature (same class).
- apex3box battery (three cubes meeting near one edge, d=1e-11..1e-6):
  conflicts=0, all resolve is2mf=1 vol correct.
- edgeTilt (two cubes abutting at x=1-1e-9, B rotated deg about z): THE
  conflict>0 EMITTERS. deg=5e-4,8e-4,1e-3 each fire conflicts=3 AND EMIT.

(iii) The money case - oracle-checked emit of a formerly-fatal input. edgeTilt
deg in {5e-4,8e-4,1e-3} (conflicts=3, i.e. pre-f0dab0fb these DIED on
EngineIdConflict) now emit: is2mf=1, comps=1, genus 0/0, vol matches Boolean
oracle a+b to dVol 5.8e-15 / 9.3e-15 / 1.2e-14, and OracleCompare's 17^3
winding-number grid AGREES everywhere (test PASSED, no "oracle disagree"). The
demotion converts a former fail-closed into a CORRECT resolve - exactly what the
soundness argument predicts (winding is multiplicity-based; the conflicted label
is a passenger, so the emitted boundary is what it always would have been).

No BREAK: every reachable conflict outcome is either (a) a correct oracle-exact
emit, or (b) a fail-closed on an accepted named guard (SubEpsFeature /
NonManifoldEmission). No conflict>0 case produced self-intersecting or
wrong-volume geometry. The 69dc83e0 backstop is genuinely superseded.

Full Overlap3 suite (minus my probes) on the landed code: 54 PASS + 1 SKIP
(Gate4c on chain-plane guard) - matches the landing's claim, no regression.

## Attack 3 - spec honesty

docs/SweepEmit3D.md [WALL-B]: states the unsoundness condition crisply ("If a
future change gives srcId a 3D CONSUMER ... this demotion must be
re-adjudicated"); names BOTH residuals (HULL->SubEpsFeature/chain-plane;
GT7081->emission bad_alloc, wall-A-shaped, explicitly NOT pinned). No overclaim:
it does NOT claim hull/GT7081 resolve - hull "lands on the CHAIN-PLANE wide-run
guard" (a fail-closed), GT7081 "EXHAUSTS memory". ASCII-clean (grep -P
'[^\x00-\x7F]' finds nothing in the section or in the diff additions). No raw
counts (magnitudes only; "8*eps" is a factor). Gate4c test comment (test:580-589)
honestly explains the SubEpsFeature skip (the former guard was a dead diagnostic;
hull now lands on the chain-plane guard, verified on real cap content).

## Attack 4 - counter semantics

cnt.engineIdConflicts is documented at its declaration (overlap3.h:58, "total
engine id-conflict events") inside a struct headed "Non-fatal counter
accumulator" (overlap3.h:55), so the non-fatal framing is set. No code/test
branches on it or surfaces it as an error. COSMETIC nit only: the field comment
does not itself say "benign / unconsumed (see [WALL-B])", which a future reader
skimming just that line could misread as an error tally; the struct header
mitigates. NEED-CHANGE at most, not load-bearing.

## Residual soundness note (not a hole in the demotion, but the honest edge)

The single load-bearing precondition is srcId staying UNCONSUMED in 3D. This is
the spec's own stated unsoundness condition and it is currently TRUE (grep: one
writer, zero readers). If per-source property transfer / provenance-attributed
emission ever reads srcId, the -1 label becomes load-bearing and this must be
re-adjudicated. Recorded, not a present defect.
