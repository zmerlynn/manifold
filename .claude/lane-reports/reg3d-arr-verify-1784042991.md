# reg3d ARRANGEMENT-CRUCIBLE stage-1 landing - ADVERSARIAL VERIFICATION

Verifier notebook (no commits). Target: explore/sweep-plane-3d-v5 HEAD a14f0f00
(743df61a clean-face per-face + graze-sound; a14f0f00 Stage-2/3 notebook), atop
40fa795d. Isolated rsync copy, fresh cmake -B vbuild -G Ninja Release -j8. All
mutations in the copy only; canonical byte-untouched (git status clean, HEAD
unchanged). Copy reverted to canonical src at close.

VERDICT: PASS. All five audits hold empirically.

## Audit 1 - THE RE-SAMPLING FIX (decision-completion, not heuristic)
EmitCleanFaces (src/overlap3.cpp ~2859-2880): centroid-first over a STATIC
constexpr kBary[6] barycentric table, break on first non-grazing sample, else
`return false`.
- DETERMINISTIC: fixed constexpr table, no RNG, points are barycentric combos of
  the face's OWN vertices A.tri[t] -> input-derived, same input -> same points.
- CERTIFIED-OR-CLOSED: any successful sample yields a directly-measured winding
  (soundness rests on the constant-cell claim: no face crossed at height eps over
  an uncrossed post-stage-5 triangle). No guess: grazes-everywhere -> fail closed.
- ADVERSARIAL all-graze: the code path exists (`if(!g) return false`); constructed
  the degenerate empirically via MUTATION-A (force centroid-only): BarsCrossZ hits
  it -> fatal "clean-face winding probe was filter-uncertain (SoS) - fail-closed".
  So the all-graze branch fails closed, never emits unverified geometry. Confirmed.
- MUTATION (break re-sampling = centroid only, rebuilt): TJunction + BridgedCaps +
  BarsCrossZ REGRESS to fail-closed exactly as the lane found; RetainsEvertedBoundary
  stays GREEN (re-sampling is orthogonal to per-face). The landed multi-sample
  CURES all three (green on landed). Fix is load-bearing and real.

## Audit 2 - BarsCrossZ FailClosed -> Resolves (independent re-oracle)
Wrote an INDEPENDENT oracle (VerifyBarsCrossZ_IndependentOracle): ANALYTIC box-union
membership (bar1 x[0,6]y[2,4]z[0,4], bar2 x[2,4]y[0,6]z[1,3]) - independent of the
committed test's atan2 GWN - with MY OWN 40000 points (seed 0xA11CE5EED) and a
FRESH eps (EpsilonFromScale ÷777, neither the pin's 1000 nor 500).
- Result: resolves, single solid, not self-intersecting, Volume = 64.0 EXACT
  (48+24-8), 6694 inside, disagree = 0 against analytic truth.
- CURED not masked: the graze is dodged by sampling the constant cell above the
  uncrossed clean triangle (a bounded decision-completion for CLEAN faces only);
  the SEAMED per-cell centroid graze is still open+documented (doc lines 286-290).
- Old pin intent preserved: under mutation A/B the carrier fails closed again with
  the winding-probe residue -> the fail-closed contract is intact for the uncured
  path. New pin mutation-verified (both mutations RED it). The conversion is
  oracle-STRENGTHENING (was: assert fatal; now: assert manifold + vol + GWN +
  decompose==1 + tol-invariance), not acceptance-widening.

## Audit 3 - THE REFUTATION TRAIL (s7b bitwise claim)
- s7b (40fa795d, read-only) claimed the per-face keeper form "bitwise-identical on
  siA/siB/BridgedCaps/SlantPlug". The arr notebook HONESTLY records BOTH sides:
  the FIRST-TRY naive per-face was bitwise-safe on siA/siB but REGRESSED
  BridgedCaps+TJunction to fail-closed -> s7b's "on all resolving fixtures"
  overclaim refuted; Design C added to keep them green. Both sides recorded.
- FNV RE-CHECK (mine, over the final welded MeshGL64 = strictly downstream of
  EmitCleanFaces, so equality proves clean-face-set identity):
    landed siA fnv=2221783380158756670  siB fnv=16197787909251774317
    flood  siA fnv=2221783380158756670  siB fnv=16197787909251774317   -> IDENTICAL
  Landed Design-C == flood on siA/siB, bitwise. (My scheme differs from the
  notebook's cited emitted-coord numbers; I verified the EQUALITY, the load-bearing
  claim, end-to-end.) Centroid-only mutation gives the same siA/siB hashes too
  (no graze on mild self-intersectors).
- Fixed regressions BridgedCaps + TJunction GREEN on landed.

## Audit 4 - TRIPWIRE + RESIDUES
- >2-sheet radial rule NOT built: grep of the landing diff adds NO Orient3DExactSign
  / ExactSign / int256 / radial call site ("NONE added"); the 9 existing predicate
  refs are pre-existing (one predicate, filter-0 site + probe hook). Second call
  site absent.
- Tripwire surfaced with owner-decision framing: docs/Regularize3D.md lines 275-282
  "EXACT-KERNEL SURFACE, QUEUED FOR REVISIT (owner contract, reluctant acceptance)
  ... ONE predicate at ONE call site ... a SECOND ... VENDOR Shewchuk's
  predicates.c." Lane Stage 3 escalates the radial rule to it (STOP + surface).
- Fail-closed carriers, NARROWER names + NO partial output (all run on landed):
    Constructed(PokedCube): NonManifoldEmission, impl absent, regularized=0, failClosed=1
    GT7863: NonManifoldEmission, whole-compose suppressed (no partial)
    openscad (ulimit -v 6G, timeout 900): fail-closed 164ms, no OOM
    EntangledBars: fail-closed (>2-sheet triple point, A.ok false pre-emission)
    GT7081 (ulimit -v 4G, timeout 900): ArrangementBudget, "retained section content
      exceeds budget", 32s under 4G (v3 sweep path, untouched)

## Audit 5 - ACCOUNTING
- FULL FILTERED RUN (ulimit -v 6G, timeout 900): Overlap3.* = 84 tests, 83 PASSED,
  0 FAILED, 1 SKIPPED (Gate4c_HullMask: pre-existing conditional GTEST_SKIP on
  arrangement-budget in the untouched hull/sweep path; 3 of the 84 are my own
  green instruments). INNER_EXIT=0.
- Named delta: exactly ONE contract change - BarsCrossZ FailClosed->Resolves
  (removed assertions are all the old fail-closed pins, replaced by oracle checks);
  RetainsEvertedBoundary is purely additive. No other EXPECT/ASSERT loosened.
- v3 sweep byte-untouched: overlap3_sweep.cpp not in the diff; landing touches only
  {notebook, doc, overlap3.cpp, overlap3.h, overlap3_test.cpp}.
- ASCII: diff (code + doc + notebook) has zero non-ASCII bytes.
- halfedge vocab: notebook's open-boundary anatomy uses halfedge/openHE; per-face
  emit legitimately traverses no adjacency (no edge-vocab regression).

## RESIDUES / HONESTY
- Soundness of the constant-cell claim (winding constant across an uncrossed
  post-stage-5 triangle at height eps) is a DESIGN argument, not machine-checked;
  it is the basis for "any successful sample is the true winding". The fail-closed
  floor bounds the downside (grazes-everywhere never emits). Fine for a landing.
- My FNV scheme != the notebook's cited numbers (different hash target); I verified
  the invariant (landed==flood on siA/siB), which is what the audit asked.
