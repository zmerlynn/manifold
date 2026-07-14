# reg3d C-2bx VERIFY - lean adversarial verification (evidence notebook)

Independent verifier. Canonical READ-ONLY; all builds/runs in an rsync copy under
scratchpad (vcopy, fresh `cmake -B vbuild Release`, -j8). No commits from this
lane. Target: branch explore/sweep-plane-3d-v5 HEAD e2e00600 (2 commits unpushed:
2c9390f8 escalation + e2e00600 notebook close). ASCII; halfedge vocab.

VERDICT: PASS. All four audits hold on rebuild-from-clean. 2c9390f8 is a filter-
first, zero-new-arithmetic PURE REFINEMENT that closes the GT7081 winding-probe
residue and honestly re-routes it to the pre-existing C-1a emission wall.

## Method
- rsync copy, configure+build clean (rc=0).
- Scoped filters during audits; ONE full suite run at close.
- GT7081 ulimit -v 4000000; openscad ulimit -v 6000000; timeout 900 throughout.
- Instrumented a throwaway build (WindingAt fire/zero counters + per-call
  RemoveOverlaps3D/ResolveComponentDirect dump of fires/zeros/fatal/detail) to
  reproduce fire counts and read fatal strings; then a mutation build (escalation
  no-op'd); then `git checkout -- src/overlap3.cpp` to restore pristine for the
  full run. Copy reverted at close.

## Audit 1 - THE ESCALATION (kernel-adjacent)
- FILTER-FIRST, verified in code (overlap3.cpp WindingAt): each of the 5
  escalation statements is gated `if (da==0)/(db==0)/(o1..o3==0)` - exact fires
  ONLY on a filter-0. Same holds for the other two production callers:
  Orient3DSoS (`if (s != 0) return s;` before exact) and EdgePiercesTriSoS
  (`if (fu==0 && fv==0 && ...)`).
- SAME PREDICATE, ZERO NEW ARITHMETIC: all 5 escalation calls invoke the
  identical `Orient3DExactSign`. The commit touches only that function's COMMENT,
  not its body (grep + diff). Exact args mirror the filter args verbatim
  (a,b,c,p / a,b,c,seed / p,seed,{ab,bc,ca}) - no reordered/altered predicate.
- FIRE COUNTS REPRODUCED (instrumented): GT7081 fires=5521 zeros=0 (every graze
  decidably NONZERO) - matches the notebook's 5521/0 exactly. siA=0, siB=0,
  openscad=0 fires.
- GENUINE EXACT-ZERO -> nullopt/re-seed, behavior identical pre/post: live tie
  fixture BridgedCaps fires=96 zeros=96 - all 96 escalations hit a genuine
  exact-zero, returned nullopt, RobustWinding re-seeded, and it STILL resolves.
  Identical to the old filter-0 -> nullopt path.

## Audit 2 - PURE-REFINEMENT CLAIM
- The escalation changes a WindingAt return value ONLY when fires > zeros (exact
  decides a filter-0 as NONZERO). Measured on every resolving fixture, fires ==
  zeros or fires == 0: siA/siB = 0, EntangledBars = 0, CoplanarFold Mult2/Mult3 =
  0, BridgedCaps = 96==96. So on resolvers the escalation returns the SAME nullopt
  the parent (90b4b247) returned => output bit-identical BY CONSTRUCTION. This is
  strictly stronger than an FNV spot-check: the executed return values are proven
  unchanged, not merely hashed-equal. Pass-through pins (Cray/Havoc/Offsets) that
  ALREADY assert bit-identical vertex multisets pass on the pristine build (full
  run).
- Only behavioral delta = GT7081 (5521 fires, 0 zeros). MUTATION-VERIFIED: no-op
  the 5 escalation lines, rebuild -> GT7081 reverts to fatal=2
  DirtyComponentUnresolved, detail "seam sub-face arrangement not exactly
  resolvable ..." and the pin REDs (NonManifoldEmission expected; "unresolvable
  sheet contact" absent). Restore -> green. The reason-string flip is exactly as
  claimed and load-bearing.

## Audit 3 - THE HONEST WALL
- GT7081 now dies at NonManifoldEmission, detail "unresolvable sheet contact"
  (instrumented dump + the pin's own substring assertion, both green).
- GT7863 dies at the SAME site: instrumented fires=0 fatal=1 detail="unresolvable
  sheet contact". Same string, same FatalReason - GT7081 reaches it with 5521
  fires, GT7863 with 0; identical terminal wall. Claim holds.
- openscad UNCHANGED: fires=0, fatal=2, detail "... non-2-endpoint seam
  (degenerate incidence) ..." - fails UPSTREAM at RecordSeams; the winding probe
  is never reached. Matches the notebook.
- LEDGER reclassification accurately recorded (docs/Regularize3D.md): table row
  GT7081 3 WINDING-PROBE(2b) -> 1 EMISSION (1a); Cluster 1a text now "GT7863 AND
  GT7081"; Cluster 3 winding-probe residue marked CLOSED (reg3d-c2bx); closure
  plan item 2 "CRUCIBLE C-2b (LANDED escalation; residue -> C-1a)".
- TRIPWIRE caller inventory matches code: 3 production callers - Orient3DSoS
  (738), EdgePiercesTriSoS (777-778), WindingAt (1003-1012) - plus the test probe
  (Orient3DExactSignProbe) and the BANKED (unbuilt) >2-sheet radial rule.

## Audit 4 - ACCOUNTING
- FULL SUITE: see close line below (single run, ulimit -v 8000000, timeout 900).
- NOTHING WEAKENED: the pin flip is reason-NARROWING - both outcomes are
  fail-closed with NO output; NonManifoldEmission is a more specific / later wall
  than the broad DirtyComponentUnresolved decline bucket. Zero oracle-wrong: the
  resolve pins carry independent volume-band + tol-invariance + GWN oracles and
  all pass; the only fixture whose behavior changed (GT7081) still fails closed,
  narrower.
- v3-era fixtures: covered by the full run, unaffected.
- ASCII: the committed code+test+docs diff and both notebooks are pure ASCII
  (grep -P '[^\x00-\x7F]' empty). halfedge vocab preserved (halfedge_.Start).

## CLOSE
- Full suite (ulimit -v 8000000, timeout 900): 579 tests / 17 suites, ALL PASSED,
  zero FAILED (241.7s). Matches the claimed 579-green.
- Copy reverted (git checkout -- src/overlap3.cpp; git diff --stat empty).

VERDICT: PASS - land/push cleared.
