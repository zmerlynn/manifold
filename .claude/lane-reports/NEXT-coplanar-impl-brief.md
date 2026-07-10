# NEXT: coplanar-family implementation brief

DESIGN: docs/SweepEmit3D.md section "COPLANAR (extension design,
round 2 - post D1/D2 fold)" - crucible-closed (no unrefuted BREAK;
trail in coplanar-crucible-1783698000.md). Probe bypasses were TEMP
and are reverted - the probe driver survives at
/tmp/coplanar_probe.cpp (rebuild against build/src/libmanifold.so,
-DMANIFOLD_PAR=1).

IMPLEMENT (house discipline: gates + mechanism pins RED FIRST, then
mechanisms; fence rule absolute):

1. GATES first (all red until mechanisms land): oracle gates
   stacked-perpendicular (P1), stacked-parallel + shared-wall (P2),
   same-oriented partial overlap (P2c), three-face group,
   mixed-orientation (+ + -, - - +), inverted-shell stacking,
   closed-shell folded FLAP (single solid, expects regularized
   single cover), rotated skeleton pin (mechanism 4 red-first),
   eps-chain (3 faces), razor-band fixture (2-10 eps separation -
   records behavior, decides threshold-vs-guard), group corner
   within eps of a vertex plane. Gate4d evolves to
   resolve-with-oracle-or-named-guard; Gate4c skip narrows to the
   open-walk guard, MustResolve on empirical success.
2. Mechanism 1: grouping pre-pass in FindSeams (symmetric-OR
   membership, union-find, INCLUDES shared-edge pairs, before the
   exemption). CoplanarOverlap retires from the taxonomy.
3. Mechanism 2: group-id seeding in BuildSlabs (ids nFaces+k).
   EngineIdConflict narrows to the near-coplanar guard.
4. Mechanism 4 (before 3 - criticals feed slabs): in-plane
   SKELETON crossings -> criticalXs (member edges + in-plane seams
   of transversal-x-member; pairwise crossing x's).
5. Mechanism 3 (the big one): per-vertex extension resolution
   replacing per-piece ExtendPtWithSeams - resolve once per
   section vertex; candidates class-i face edges / class-ii seams /
   class-iii in-plane member edges / weld; 3D-identity preference
   at criticals (extension IS the arr.verts vertex when the
   governing track terminates there); deterministic ties (lowest
   class, lowest member id); unresolvable -> fail closed.
   BuildCapEdgeSet and strips consume the SAME resolved table.
6. EdgeInPlane fatal retires; the FindSeams flap comment updates to
   the resolve-as-+2 re-adjudication (the doc records it).
7. ACCEPTANCE: Gate4c hulls. D1's instrumented case (cap
   x=-107.465, verts 32/45) is the canary - its junction must
   resolve to ONE vertex. Expect more real-geometry surprises
   behind it; the open-walk guard is the honest boundary for
   whatever remains.

Suite at brief time: 35+1 Overlap3, 580+1 full, all green on the
reverted tree. Crucible notebooks: coplanar-crucible-1783698000.md
(orchestrator, incl. the D3 attacks executed main-agent after two
lane output-cap deaths), lane-D1/D2 notebooks in the
/tmp/codex-branch-review trees.
