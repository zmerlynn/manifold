# VERIFY lane: phantom-crossing hardening (e55082e1), Regularize3D

Branch explore/sweep-plane-3d-v5. Gate the push of the exact seam-crossing
predicate. Range 4a4d1564 (pushed, verified) .. HEAD. READ-ONLY except this
notebook. Grounding: .claude/lane-reports/phantom-x-1784140994.md.

## VERDICT: PASS (push-safe)

Every claim reproduced from a from-scratch Release build at BOTH SHAs. Do NOT
push performed by this lane (verify-only).

## Range hygiene
- Only code commit in range: e55082e1 (overlap3.cpp +docs). `git log` on
  src/test/include/CMakeLists confirms it is the sole non-notebook commit.
- Range src+docs diff == exactly the e55082e1 diff (overlap3.cpp, Regularize3D.md).
- No test/ change in the range (same 32 Overlap3 tests).
- base 4a4d1564 src == 2b509faa src (the grounding's base): identical.
- Mid-lane HEAD advanced 37aef4f3 -> 52e42dcb (concurrent attack22 NOTEBOOK-only
  commit). src/docs/test byte-identical since 37aef4f3; 4a4d1564 and e55082e1
  both still ancestors of HEAD. The reviewed code is unchanged; verdict stands.

## Builds (from-scratch Release, Unix Makefiles, -O3 -ffp-contract=off
## -fexcess-precision=standard, gtest src reused from build-base/_deps)
- /tmp/phantom-verify-base (4a4d1564) -> /tmp/pv-base-off  (PAR=OFF)
- /tmp/phantom-verify     (37aef4f3) -> /tmp/pv-head-off   (PAR=OFF)
- FNV copies (add_fnv.py hook on COPIES only; canonical never instrumented):
  pv-base-fnv, pv-head-fnv (PAR=OFF), pv-head-par (PAR=ON). Hook presence
  confirmed in each .so (2 ZZZFNV markers) and ABSENT in pristine pv-head-off
  (0) - rules out the grounding's stale-rebuild trap; all binaries fresh-linked.

## 1. Suites (exact counts, no skips)
- HEAD Overlap3 PAR=OFF: 32/32 PASSED, 0 SKIPPED, 0 FAILED (139.5s).
- BASE Overlap3 PAR=OFF: 32/32 PASSED, 0 SKIPPED (222.4s).
- HEAD Overlap3 PAR=ON : 32/32 PASSED, 0 SKIPPED (86.8s).  [both PAR settings]
- GT7081 alone, ulimit -v 4000000, timeout 900: PASSED both trees.
- openscad, ulimit -v 6000000, timeout 900: PASSED both trees. The fail-closed
  PIN asserts refusal (r.fatal == NonManifoldEmission, "unresolvable sheet
  contact", no impl) - test OK == refusal held.
- openscad census (F4B_DUMP, large component):
    BASE F4B_TRIPLES distinct=168 incidences=962 ; openEdges=22 ; okfaces=574
    HEAD F4B_TRIPLES distinct=162 incidences=948 ; openEdges=22 ; okfaces=574
  => 6 phantom triples removed (168->162); openEdges=22 UNCHANGED (fail-closed
  census identical); fan census (okfaces=574) byte-identical.

## 2. FNV base vs HEAD (independent, full Overlap3.* with ZZZ_FNV, PAR=OFF)
- 52 FNV lines each, DIFF EMPTY, same SHA1 e4b3a364. Byte-identical on every
  resolving carrier including GT7081 (nv 15754 nhe 94380) and siA/siB-scale
  (nv 8474). GT7081's large line present in both sets -> not OOM-skipped.

## 3. Diff review e55082e1
(i)  sos:: namespace byte-identical base vs HEAD (same SHA1 1fb1eb87); diff
     touches no `sos` line.
(ii) axis-drop projection argument (read + derived):
     - ExactOrient2DDrop drops the DOMINANT normal axis (DominantAxis = argmax
       |faceN component|), so |n[axis]| != 0 and the axis-drop projection of the
       (coplanar, in-face-plane) seam endpoints is NON-DEGENERATE (bijection
       plane->R^2). Orientation parity: proj may be orientation-reversing and
       the z=+1 lift adds a further -1 factor, but ExactSegProperCross uses ONLY
       relative signs (o1!=o2, o3!=o4), invariant under any global sign flip ->
       crossing decision preserved.
     - All predicate args are the shared 3D seam endpoints A.faceSeams[f][k].p0/p1
       (RecordSeams once-only constructions, pushed once per face pair L1754-55),
       read directly - no fresh per-call construction. The projection is pure
       coordinate SELECTION + a unit lift (no arithmetic, no new rounding), so
       "no exact-on-constructed" holds: only level-0 doubles enter the ONE
       blessed FORM, FILTER-FIRST (Orient3DFilterSign -> Orient3DExactSign).
     - seg[]/pf.proj retained ONLY for SegLineIntersect2D (the split position x,
       a double keyed by the once-only Intersect3Planes pos, NOT a decision).
     Semantics vs base: removes PlaneFrame-projection rounding + la::cross
     roundoff; exact w.r.t. the stored endpoint doubles. Byte-identity on the
     resolving corpus is the safety net that the axis-drop choice flips no
     resolving decision (empirically confirmed, section 2).
(iii) inventories BOTH updated in e55082e1: src tripwire overlap3.cpp ~L782 and
     docs/Regularize3D.md ~L376; consistent phrasing, honest (byte-identical on
     resolving corpus, removes only openscad phantom triples, filter-first).
(iv) SegProperCross2D DEAD: only a past-tense comment mention remains anywhere
     in the HEAD tree; no live caller.
(v)  no leftover instrumentation (new helpers carry no getenv/fprintf; F4B_DUMP
     is pre-existing measurement); diff pure ASCII (src + docs).

## 4. Determinism spot
- FNV 1t vs 4t (PAR=ON, taskset -c 0 vs 0-3; nproc confirms 1 vs 4) on
  siA + GT7081: 4 lines each, DIFF EMPTY, same SHA1 10ffe028 - thread invariant.
- openscad census stable across 2 HEAD runs (distinct=162, openEdges=22 both).

## Note
Verify lane only. This notebook is the sole commit; no push.
