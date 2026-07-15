# VERIFY lane: junction-registry landing (gate for the push)

Range fb0fadc7 (pushed, verified) .. HEAD 96326c11. Src/test commits under gate:
b42e0adf (both-sides retention), aa18c77c (global junction registry),
7a728385 (ledger + openscad pin narrowing). Adversarial, from-scratch Release
builds at BOTH SHAs (MANIFOLD_PAR=ON, -j8): /tmp/junction-verify (HEAD),
/tmp/junction-verify-base (fb0fadc7). Canonical READ-ONLY; this notebook is the
only commit. NOT pushed.

VERDICT: PASS (push-safe).

## 1. Full Overlap3 suite, both trees (exact counts)

- Bulk (32 total minus the 2 heavy): BASE 30/30 PASSED, HEAD 30/30 PASSED. No
  FAILED, no SKIPPED (all models present).
- GT7081 (ulimit -v 4000000; timeout 900): BASE OK (63.9s), HEAD OK (76.0s).
- openscad (ulimit -v 6000000; timeout 900): BASE OK (306ms), HEAD OK (473ms) -
  the pin PASSES and still asserts the REFUSAL (NonManifoldEmission +
  "unresolvable sheet contact").
- Total: BASE 32/32, HEAD 32/32.

## 2. Independent FNV base vs HEAD, all resolving carriers (memcmp)

Standalone harness (reg_fnv), same source compiled against BOTH build trees,
reproduces the test carrier constructions exactly (free ReadOBJ + compose) and
hashes the resolved Impl on BitIdenticalMesh's bit-basis (raw vertPos_ +
halfedge Start + sizes). 1 thread. base == HEAD, every carrier:

    pokedcube e82e444cca63307e   entangled 38c8778de9a90df2
    siA       e9fd9fb18ac2ce71   siB       ad15857512fd41a2
    gt7863    d83bc1bde57feb85   gt7081    cf324fdeda3a4e71
    openscad  FATAL reason=1 (NonManifoldEmission) "unresolvable sheet contact"

All six resolving carriers byte-identical; openscad fail-closed identical. The
both-sides + registry changes are a byte-for-byte no-op off openscad.

## 3. Diff review of the 3 src/test commits

(i)   sos:: byte-identical base vs HEAD (219==219 lines, diff clean). Kernel
      exact-sign paths (Orient3DExactSign / ExactOrient3D / EdgePiercesTriSoS)
      untouched in the range. Registry math is la::dot/length (level-0 double);
      no new exact arithmetic on constructed points.
(ii)  both-sides rule (EmitSeamedFace): probes w_above AND w_below; drops when
      aboveIn==belowIn (both inside or both outside), emits when exactly one,
      oriented to the solid (-nHat side -> original, else reversed). REDUCES to
      the old w_above==0 rule at jump==1: for mult-1, w_below=w_above+1 so
      aboveIn!=belowIn holds EXACTLY at w_above==0, where belowIn=true ->
      original orientation == the former rule byte-for-byte. (Empirically the
      no-op above.)
(iii) registry keying deterministic/order-independent: raw endpoints+triples
      sorted by (x,y,z) then greedy within-eps dedup to the sorted-lowest rep -
      face iteration order cannot change the result. The split DECISION is a
      pure 3D on-segment test on the segment's SHARED 3D endpoints, so both
      incident faces of a seam/mesh edge decide identically (endpoint swap ->
      t<->1-t, symmetric interval; same accepted vertex set). Empirically:
      1t==4t determinism below.
(iv)  pos2in acceptance keyed-only: add/addAt key vidx by the 3D tuple
      {P.x,P.y,P.z}; registry split points are EXISTING once-only constructions
      (seam endpoints, Intersect3Planes triples), so no new preimage enters the
      acceptance map. TriangulateIdx adds no Steiner points.
(v)   fail-closed arms intact: 12 `ok=false` arms preserved (12==12); the seamed
      probe now fails closed on !g||!gb (stricter); +1 NEW arm - a try/catch
      returning false on a malformed clean-face split polygon. No arm removed.
(vi)  test-pin change (7a728385, overlap3_test.cpp): comment-only. Every changed
      line is a `//` comment; TEST body and assertions unchanged.

## 4. Load-bearing mutations (F4B_DUMP census, openscad large dirty component)

    both-sides OFF (BASE build)        : openEdges=91  ties=9  overlaps=5
    registry OFF (HEAD, F4J_NOREG=1)   : openEdges=94  ties=0  overlaps=0
    both ON (HEAD)                     : openEdges=30  ties=0  overlaps=0

Small dirty component openEdges=0 throughout. both-sides lever: ties 9->0,
overlaps 5->0. registry lever: 94->30 (registry off reopens to 94). Trajectory
91 -> 94 -> 30 exactly as claimed. Both levers load-bearing.

## 5. Determinism

HEAD, tbb global_control 1 vs 4 threads:
    siA    1t==4t  e9fd9fb18ac2ce71
    gt7081 1t==4t  cf324fdeda3a4e71
openscad refusal deterministic across 2 runs (FATAL reason=1, same detail).

## Notes

- Build config PAR=ON (matches the 32/32 rail); FNV/mutation/census under
  single-thread or tbb global_control for a clean isolation of the code no-op
  from thread nondeterminism. Determinism rail separately exercises 1t vs 4t.
- No push. Canonical tree unmodified except this notebook.
