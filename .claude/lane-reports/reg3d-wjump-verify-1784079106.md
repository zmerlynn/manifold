# reg3d-wjump / shares-vertex-recovery ADVERSARIAL VERIFICATION notebook

Verifier lane. Branch explore/sweep-plane-3d-v5, HEAD 9b1fbf0d (2 commits atop
1e55b700, UNPUSHED). Canonical READ-ONLY; all work in an rsync copy under
scratchpad/verify-copy with a fresh cmake -B vbuild Release. This notebook is the
only artifact I write; no commits. ASCII, halfedge vocab.

Under audit:
- 9b1fbf0d reg3d-wjump: shares-vertex genuine-crossing recovery (src/overlap3.cpp
  RecordSeams ~56 lines; test flips; docs ledger).
- b7c97708 reg3d-wjump: notebook + Gap ledger bake-off corrections (docs+notebook).

The recovery (RecordSeams, lines ~1169-1214 + ~1378-1389): when two tris share a
vertex, do NOT skip if an edge with NEITHER endpoint a shared vertex genuinely
pierces the other tri's interior (EdgePiercesTri==1, or EdgePiercesTriSoS==1 on
the -1 residue). Recovered seam = [V, offVertexP], V the shared corner supplied
as pts[1] when nPts==1 && nShared==1 && |P-V|>eps.

## Build
- rsync copy + fresh cmake Ninja Release (MANIFOLD_TEST/PAR/DEBUG/ASSERT ON,
  matching canonical CMakeCache). Deps reused from canonical _deps/-src (read-only).
- manifold_test built -j4 (4 cores/7GB; -j8 would thrash), 68/68 linked OK.


## AUDIT 1 - THE RECOVERY'S FALSE-POSITIVE DIRECTION (sharpest): PASS
Method: 8 adversarial genuinely-adjacent, NON-crossing valid solids driven through
ResolveComponentDirect (which FORCES RecordSeams on any soup). stderr census
counter added at `recoveredSV = true` (arrangement-only; behaviorally inert -
verified FNV unchanged vs pristine). FNV compared HEAD vs parent (recovery off).
Fixtures: Sphere_32, Sphere_80 (high-valence poles + valence-6), Cylinder_64/128
(COPLANAR cap fans - stress the SoS edge-in-plane guard), ConeCyl_48 (apex fan),
Cube, Tetrahedron (hub corners), StarExtrude_6 (reflex "bellows" fold-backs).
RESULT: recovery fires ZERO on all 8 (attributed census: every fixture 0; the
PokedCube positive control fires exactly 4). Every one of the 8 round-trips
EXACTLY (outVol==inVol, outSI=0) and its output FNV is BITWISE-IDENTICAL parent
vs HEAD:
  Sphere_32 6184d5db66200e90  Sphere_80 b0b7042b085516c1  Cyl_64 5a5abdee29114301
  Cyl_128 4f303f4112b8cef3  ConeCyl_48 d8c1a64d5db1045d  Cube a17ffe5f57b07c4f
  Tet e505d27a238c1aab  StarExtrude_6 b57e8cfeca389fe7   (all identical both libs)
The SoS fallback's edge-in-plane guard (overlap3.cpp:777) makes coplanar cap-fan
edges return 0, so no legitimate adjacency manufactures a phantom pierce.
=> No false positive on legitimate adjacency; recovery is purely additive on
genuine off-vertex transversal crossings.

## AUDIT 2 - THE THREE FLIPS, independently: PASS
Independent MC-GWN oracle (my own Van Oosterom-Strackee + RANDOM uniform sampling,
distinct from the test's grid GwnVolumeOracle):
- PokedCube: RemoveOverlaps3D resolvedVol = 0.250000 EXACT; my MC-GWN(4M pts,
  2 seeds) = 0.249187 / 0.249984 (~0.25, MC noise); tol-invariant across
  eps*{0.25,0.5,2,4} all = 0.250000000; outSI=0, Manifold::Status NoError.
- GT7863 8-edge-hole component (comp#0, 216 tris, SELF-INTERSECTING): resolves,
  signedVol 855.014, resolved vol 855.036, tol-invariant (eps and eps/2 both
  855.0359), outSI=0.
- GT7081 both dirty shells: full test resolves (96.9s under ulimit -v 4G),
  vol-preserved + tol-invariant + per-component non-SI.
MUTATION (install parent 1e55b700 overlap3.cpp = recovery disabled, rebuild):
all three revert to fail-closed - PokedCube -> NonManifoldEmission (no output),
DirtySingleComponent / ExactZeroTie_Constructed / MultiComponent / GT7863-hole
all -> "unresolvable sheet contact"; GT7081 -> "unresolvable sheet contact"
(37s). Every resolve is load-bearing on the recovery.
s7b single-face no-op PIN: `git show 9b1fbf0d -- src/overlap3.cpp` = exactly TWO
hunks, BOTH inside RecordSeams; the only "fold" token in the diff is a comment.
EmitCleanFaces/EmitSeamedFace/fold retention rule is byte-untouched => the
+1-only rule is unchanged, the arrangement is what completes.

## AUDIT 3 - THE REFUTATION TRAIL: PASS
Notebook reg3d-wjump-1784074100.md records all five refuted prior diagnoses with
evidence: (sliver) "no admissible near-coplanar target; 887-2259 ULP strips are
CROSS-mesh pass-through, mislabeled"; (weld) "thinness-weld REFUTED - ZERO weld
collapses corpus-wide"; (plane-based-rep) + (kernel) "no new predicate FORM, no
kernel escalation - REFUTES ... plane-based-rep / Shewchuk-tripwire"; (double-
sheet-rule) "DOUBLE-SHEET SUBSUMED by (1) + reverted-as-unforced (s7b was right)".
CENSUS CROSS-CHECK (decisive): the recovered off-vertex crossings ARE exactly the
prior "unmatched seam pierce points." Instrumented WJUMP_SEAM dump on PokedCube:
P = {(0.0714,-0.0714,-0.5),(-0.0714,-0.5,0.0714),(0.0714,-0.5,-0.0714),
(-0.5,-0.5,-0.25)} - the SAME SET (order aside) as the notebook's four points;
all nShared=1.
Ledger deltas accurate: carriers table PokedCube + GT7081 -> RESOLVES wjump;
GT7863 -> 1 EMISSION (1a') DirtyComponentUnresolved (narrowed); closure plan
point 3 C-1a CLOSED + bake-off correction; point 4 C-1b CLOSED.

## AUDIT 4 - GT7863 comp#1 HONEST WALL: PASS
Whole compose: fatal=DirtyComponentUnresolved, hasImpl=0 (NO partial output),
comps=4 clean=2 dirty=2 regularized=1 failClosed=1.
comp#1 (144 tris, NOT self-intersecting, within-component coplanar): fails closed
DirtyComponentUnresolved, no partial; detail = "exact-coplanar fold declined
(coplanar/transversal entanglement...)". My independent measurement: signedVol
2829.920 and MC-GWN{w>=1} 2826.5 AGREE at ~2830 => a valid single-sheet {w>=1}
of vol ~2830 provably EXISTS, yet the doubled connectivity (fold self-check m=2
vs winding jump 1) is not a representable single-sheet manifold => honest
representability wall, refused not silently mis-resolved. Banked fold-escalation
prototype is NOT landed: grep of overlap3.cpp for
fold.?escalat/reliability-certified/wide-probe/trust.*winding = empty.

## AUDIT 5 - ACCOUNTING: PASS
Full Overlap3 suite: 32/32 PASS (one run at close, ulimit -v 6000000, openscad
fail-closed 203ms, GT7081 resolves 96.9s, siA/siB in the MC volume band). Zero
failures. Change is ISOLATED: overlap3.h/RemoveOverlaps3D/ResolveComponentDirect
referenced ONLY by overlap3.cpp + overlap3_test.cpp + src/CMakeLists.txt - no
boolean/manifold/impl path calls it (header: "public wiring is post-prototype"),
so no other suite can be affected; nothing weakened. FNV bitwise-identical on
untouched (zero-recovery) resolving fixtures (Audit 1). All artifacts pure ASCII
(overlap3.cpp/.h, ledger, notebook, test: 0 non-ASCII lines). Halfedge vocab
(recovery reads A.vid/A.tri from in.halfedge_.Start). openscad measured once.

## VERDICT: PASS (all 5 audits green)
The shares-vertex genuine-crossing recovery is oracle-correct and does NOT
false-fire on legitimate adjacency: zero fires + bitwise-identical output on 8
fan/bellows/high-valence/coplanar-fan solids; the 4 PokedCube fires are exactly
the previously-unmatched pierce points. The three carrier flips are independently
oracle-graded (my own MC-GWN) + tol-invariant + mutation-anchored; the retention
rule is byte-untouched (s7b holds). GT7863 comp#1 is an HONEST fail-closed wall
(GWN-proven ~2830 single sheet exists but the coplanar double-sheet is not
representable; escalation banked-unbuilt, confirmed absent). Full suite green,
change isolated, ASCII. Push gate: CLEAR from this verifier's seat.
Canonical READ-ONLY preserved; all experiments in an rsync copy; copy reverted.
