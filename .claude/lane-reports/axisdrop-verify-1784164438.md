# VERIFY lane: axis-drop landing (REBLESS) - bc350f70..4423266e

Adversarial gate for the owner-approved re-bless landing (f94d005e seamed
axis-drop, 3f098522 fold+clean migration + PlaneFrame deletion, 4423266e docs).
Canonical READ-ONLY except this notebook. Grounding: axisdrop-1784159614.md.
All builds FROM-SCRATCH Release (cmake configure + build) in /tmp worktrees at
the exact SHAs: /tmp/axisdrop-verify (HEAD PAR=ON), -base (base PAR=ON),
-nopar (HEAD PAR=OFF), -base-nopar (base PAR=OFF). 4-core box, -j4.

## VERDICT: PASS (push-safe)

## 1. Suites + budgets

- Overlap3 suite 32/32 PASSED on ALL FOUR from-scratch builds: base-PAR,
  head-PAR, head-NOPAR, base-NOPAR (runs under ulimit -v 6000000).
- GT7081 @ HEAD (ulimit -v 4000000; timeout 900): PASSED, wall 92.0s,
  RSS 85 MB. (Slower than the lane's 67s - a concurrent lane was burning CPU
  during my run; well inside budget either way.)
- openscad @ HEAD (ulimit -v 6000000; timeout 900): PASSED fail-closed,
  wall 0.53s, RSS 15 MB.

## 2+3. Independent equivalence oracle (NOT the lane's numbers)

Own harness (env-gated dump test appended to the throwaway worktree copies of
overlap3_test.cpp, both SHAs identical; canonical never touched): per carrier,
raw double bits of every output vert + tri triples + outcome class + counters.
Own Python comparator: order-independent canonical soup (min-vertex rotation,
orientation preserved) -> FNV-1a; exact fsum signed volume; closed-manifold
check (every directed edge paired); GWN membership (Van Oosterom-Strackee) on
MY OWN fixed grid - 2600 points, seed 12345, 80% concentrated near the base
surface at graded offsets (1e-4..8e-3 of scale), 20% uniform; a point counts
only if BOTH soups classify it unambiguously.

- IDENTITY carriers (bar = byte-identity): GT7863, BarsCrossZ, BridgedCaps,
  GT7081 - canonical soup BYTE-IDENTICAL base vs HEAD. PASS.
- EQUIVALENCE carriers (bar = equivalence): siA, siB, PokedCube,
  PushedCapSphere, EntangledBars, EntangledBarsRotated - all byte-DIFFERENT as
  disclosed, with: identical outcome class + dispatch counters, identical
  nVert/nTri (siA/siB 8474/16948, PokedCube 11/14, PushedCapSphere 80/160,
  both EntangledBars 32/60), volume agreement rel <= 5e-16 (exact-zero on
  siA/siB/PokedCube/EntangledBarsRotated), GWN membership 2600/2600 checked
  with ZERO disagreements on every carrier, closed-manifold both sides,
  IsSelfIntersecting identical. PASS.
- openscad: fail-closed BOTH SHAs, fatal class identical (NonManifoldEmission,
  detail "unresolvable sheet contact" = the SplitTouchingSheets refusal),
  counters identical (2 components, 2 dirty, 1 regularized, 1 failClosed).
  Census unchanged. PASS.
- Comparator false alarm, resolved: my first pass showed GT7081 vol rel
  3.2e-9 DESPITE byte-identical soups - that was my own summation noise
  (det[a,b,c] is cyclic-invariant only in exact arithmetic; base/head index
  rotations differ). Recomputed from the canonical rotation: exactly equal.
  Also my status check first flagged GT7863/GT7081 selfint=1 - that is the
  composed cross-component overlap BY DESIGN (per the tests' own comments),
  identical base vs head, and both carriers are byte-identical anyway.

## 4. Code review (all three commits, full diff read)

(i) AxisDropFrame soundness, re-derived: DominantAxis on the UNNORMALIZED
faceN (scale-invariant, so it equals nHat's dominant axis); |nHat[axis]| >=
1/sqrt(3) for any unit vector, so the lift division and the projection are
never degenerate for a valid (nonzero, finite) normal - the edge-on
DEBUG_ASSERT is unreachable except for non-finite input, the same failure
class the old frame had. Parity: the swap lives in ONE place (proj), shared by
every consumer; per-consumer audit: EmitSeamedFace has a runtime CCW guard
(parity error = fail closed, never inverted emission); FoldCoplanarClusters'
member sign s = sign(dot(faceN, nHat)) is 3D and parity-free, PointInTri2D is
orientation-agnostic, and the wjump==m cross-check plus the re-gate backstop
the emission; EmitCleanFaces re-triangulates the face's own CCW loop (a
parity flip would make it CW and Triangulate throws = fail closed);
ExactSegProperCross takes the raw axis and is chirality-invariant (strict
straddle). Lift is applied at EXACTLY two sites, both genuinely 2D-born: the
fold's getP new-crossing vertices and the fold's interior classify probe
(probe only, never emitted). Every other emitted vertex carries canon3/pts3
3D-source bits; the seamed classify point is now the DIRECT 3D centroid; the
triple-point production path is once-only Intersect3Planes (lift there only
under the F4B_PERFACE measurement env).
(ii) PlaneFrame/BuildPlaneFrame: ZERO references at HEAD (grep src/).
(iii) sos:: namespace extracted from both SHAs and cmp'd: BYTE-IDENTICAL.
Exact-predicate caller inventory: the same 10 call sites at both SHAs (line
shifts only) - no new exact callers; kernel untouched.
(iv) INCIDENT SWEEP: each commit touches ONLY src/overlap3.cpp,
docs/Regularize3D.md, .claude/lane-reports/. All 17 src hunks map to
axis-drop scope (read in full). No added env flags (getenv inventory identical
base vs HEAD; F4P_PROJ absent from canonical at both SHAs). NO FOREIGN HUNKS.
(v) Docs: L1 marked RESOLVED with the re-bless described honestly as
equivalence-not-regression; byte-different framing by kind, no counts-as-facts
in docs; landed diff is ASCII-clean. Net src lines +40 (107 ins / 67 del) -
matches the lane's claim.

## 5. Determinism at HEAD

Dump files (raw vert bits + tri order) BITWISE IDENTICAL across taskset -c 0
(1 core) vs -c 0-3 (4 cores) under PAR=ON, AND across PAR=ON vs PAR=OFF, for
ALL 11 carriers including GT7081 and both required equivalence carriers.
Stronger than the order-independent-FNV claim: even emission order is
invariant.

## 6. Mutation (the 486-class divergence, re-measured NOW)

The lane's F4P_PROJ instrumentation patch (/tmp/attack22-proj-instr.patch)
applies cleanly to base bc350f70; fresh worktree + from-scratch PAR=OFF build;
F4B_DUMP=1 F4P_PROJ=1 on the openscad carrier reproduces the divergence
EXACTLY: xTotal=33231 xAgree=32745 xPFonly=319 xADonly=167 (=486) across
268/574 faces; eps-merge projCulprit=4 normalDrop=0 mADonly=0; small comp
28/0. Matches attack22-proj and the landing notebook to the digit - the
number is current, not stale, and the rounded overlay it convicts is exactly
what this landing deletes.

## Perf spot check

siA whole-test wall, best-of-5 interleaved, 4 cores: base 7.82s -> head 7.67s
(~2% on the whole test, consistent with ~3% on the resolve fraction).
Direction confirmed, no regression.

## Artifacts

Throwaway worktrees /tmp/axisdrop-src-{base,head,mut}, builds
/tmp/axisdrop-verify*, /tmp/axisdrop-mut; harness + comparator + dumps in the
session scratchpad. Canonical repo never modified except this notebook.
