# reg3d-ent ENTANGLEMENT landing - LEAN ADVERSARIAL VERIFICATION

Branch explore/sweep-plane-3d-v5, HEAD 213d5559 (6 unpushed commits). Read-only
canonical; all builds/mutations in an rsync copy under scratchpad/vsrc (fresh
`cmake -B vbuild` Release, -j8, deps reused). Verdict gates the push.

Artifacts read: the 6-commit diff (src/overlap3.cpp +46/-3, test +175, docs +30,
lane-report), .claude/lane-reports/reg3d-ent-1784044725.md, docs/Regularize3D.md.

Baseline (copy, HEAD as-is): EntangledBars_Resolves / EntangledBarsRotated_Resolves
/ BarsCrossZ_Resolves all GREEN. Initial full build 2:22.

## Audit 1 - THE DEDUP (positional vs construction-keyed) -> PASS (with noted hardening)

The dedup (overlap3.cpp:2349) IS genuinely positional, not construction-keyed:
`la::length(pts[a]-pts[b]) <= eps`. The two symmetric double-pierce points come
from DIFFERENT pierce-cache keys (owner=i cap-edge INTO tgt=j vs owner=j cap-edge
INTO tgt=i -> different SegPlanePoint constructions), so it is a distance merge of
two distinct derivations, NOT a shared-vertex reference. The audit's sharpest
concern is CONFIRMED as to mechanism.

Why it is nonetheless SOUND (bounded, not fragile-to-silent-wrong):
- Threshold == the weld radius eps. The dedup can only merge points the downstream
  BuildImpl emission weld would ALSO merge, so it never manufactures a topological
  decision the weld wouldn't. Merging within eps is sound by the standing
  eps-valid contract.
- Colinearity: for a single triangle pair every recorded point lies on the
  plane_i ^ plane_j line, and two convex triangles meet in ONE segment (2 true
  endpoints); a 3rd/4th recorded point is a near-duplicate of a true endpoint or a
  recovery-double, never a spurious-far point (the recovery only records ON the
  seam region, inside tgt's triangle).
- Failure DIRECTION is safe: over-merge (a genuine sub-eps seam collapses) only
  drops nPts below 2 -> the nPts!=2 fail-closed. Never nPts==2-with-wrong-endpoints.
- EMPIRICAL (ZZZDEDUP merge-separation logger, build E): on the two RESOLVING
  carriers every merge is ULP-scale - EntangledBars(axis) sep == 0.000e+00 (8x,
  bit-identical); Rotated sep in {0, 2.2e-16, 3.1e-16, 4.4e-16}, i.e. <= ~2e-5 * eps
  (eps=1.1e-11). So on the resolving path the dedup removes FP twins of the SAME
  reentrant corner, not genuinely-distinct points. openscad's max merge is 0.244*eps
  (still sub-eps, well below the weld radius) and that carrier fails closed anyway.
- The direct adversary (MUT-mine, audit 3: displace the recorded endpoint by 2*eps)
  fails CLOSED ("unresolvable sheet contact"), never a silent wrong resolve.

NOTED HARDENING (not blocking): the design's own words are "once-only keyed,
construction shared" - the reentrant corner is a vertex of the fold's in-plane
arrangement (an edgeSubdiv overlap-corner). The clean form is for FoldCoplanarClusters
to expose that overlay vertex by construction identity and have BOTH walls' recovered
endpoints REFERENCE it (bit-identical), retiring the positional dedup entirely. The
current impl realizes "give the endpoint the fold's identity" via a weld-radius
positional merge instead of a shared construction; sound today, but the gap between
stated design and impl is real.

## Audit 2 - ROTATED CARRIER + WELD RELIANCE -> PASS on correctness; NEW mesh-quality finding

Independent re-oracle (build F/G, my own seed 0x5EED42, 40k points, my own GWN BOTH
sides, fresh eps): vol = 37.55674168534177 == (A+B) library-union volume EXACTLY
(absdiff 0, reldiff 0); independent-GWN checked=40000 disagree=0; selfIntersecting=0;
Decompose=1. The rotated resolve is the exact boolean union - reproduced independently.

Near-coincident-sheet scan (the R1 self-fold direction): 0 pairs within 100*eps;
min sheet gap between nearby-centroid non-shared-vertex tris = 0.218 (>> eps). No
sub-eps parallel sheets -> not on the harmful side of R1 here.

NEW FINDING (the audit's target - a silent thin artifact IsSelfIntersecting is blind
to): the rotated output carries ONE DEGENERATE (zero-area) triangle. minArea=3.8e-16,
minQual=5.2e-16; worst tri = three near-colinear verts on a wall line x=-2.14450692,
y=-1, z in {4, 3.7148, 2.2851} (x agrees to 15 digits, ~1 ULP spread -> colinear).
ATTRIBUTION (build G survey) is clean and rotation/weld-specific:
  EntangledBars(axis)  minQual=0.433  0 degenerate  (bit-exact junctions)
  BarsCrossZ(no-cluster) minQual=0.170 0 degenerate  (recovery dead)
  EntangledBars(rot50) minQual=5.2e-16 1 DEGENERATE
  REF library union of the SAME rotated bars: minQual=0.083, 0 degenerate.
So it is NOT a general candidate-B trait and NOT in the library's own union - it is the
materialized R1 weld artifact: the sub-eps mismatch between the recovered cap endpoint
and the fold overlay vertex leaves a 1-ULP-distinct wall vertex that the emission weld
does not merge, producing a zero-area sliver. It is INERT (vol/GWN/topology all exact;
contributes 0 solid angle, 0 volume) but it is a genuine mesh-quality defect that the
re-gate (IsSelfIntersecting + validity) does NOT catch - an empirical demonstration
that R1/R2(i) is LIVE. The lane's rotated pin asserts vol/GWN/one-solid/not-SI but does
NOT assert min triangle quality, and the notebook/docs claim "RESOLVES oracle-true"
without naming the sliver. Recommend: name the artifact and add a min-quality guard to
the rotated pin (and/or the re-gate). Confined to the SYNTHETIC adversarial pin - no
corpus carrier (EntangledBars axis is clean).

## Audit 3 - MUTATIONS -> PASS (all fail-closed, none silently wrong)

- MUT-recovery (build B, neutralize the cap-plane recovery hit): EntangledBars AND
  Rotated RED with "degenerate incidence" fail-closed. Recovery load-bearing.
- MUT-dedup (build C, disable the dedup loop): both RED; ZZZTRUNC histogram flips to
  nPts>=3 = 16 (8 per carrier), nPts==1 = 0 - the un-deduped symmetric double-pierce
  makes nPts==3 -> fail-closed. Dedup load-bearing, exactly as the lane claims.
- MUT-mine (build D, displace the recorded recovered endpoint by 2*eps): both RED via
  FAIL-CLOSED ("unresolvable sheet contact"), r.fatal set, no output. Never a silent
  wrong resolve. Confirms the safety direction.

## Audit 4 - OPENSCAD NARROWING -> PASS

Apples-to-apples ZZZTRUNC counter (same build, recovery toggled):
  recovery OFF (MUT-recovery) = 496 nPts==1 truncations
  recovery ON            = 344 nPts==1 truncations  -> 152 seams completed. MATCHES.
  nPts>=3 = 0 in BOTH -> zero >2-sheet junctions (the reduction holds; no radial rule).
openscad stays FAIL-CLOSED with recovery on: FatalReason::DirtyComponentUnresolved,
detail names "degenerate incidence", r.impl empty (no partial output). Residue is the
strictly-narrower near-coplanar-sliver / cap-interior wall (separate axes). Also
EntangledBars axis-aligned reproduced oracle-true: vol=40 exact (24+24-8), GWN 0/>3000,
one solid, not self-intersecting, tol-invariant.

## Audit 5 - ACCOUNTING -> PASS

- v3 sweep byte-untouched: git diff 51ab74d3..HEAD -- src/overlap3_sweep.cpp = 0 lines.
  Whole src delta is overlap3.cpp (+46/-3).
- No test weakened: test deltas are EXACTLY the two named pin conversions/additions -
  EntangledBars_FailClosed -> _Resolves (rename, coverage upgraded not dropped) and
  new EntangledBarsRotated_Resolves; openscad test kept (name unchanged), still
  fail-closed, only comment + narrower-residue prose. Zero DISABLED_/GTEST_SKIP added,
  zero coverage deletions.
- Bitwise on the uncrossed carrier (FNV-1a over raw vertPos_ bits + halfedge starts,
  BarsCrossZ, resolves in both): OLD 51ab74d3 = 0x25f465667ff7fdae, NEW HEAD =
  0x25f465667ff7fdae -> BIT-IDENTICAL. The recovery is structurally dead there (no
  clusters) and the dedup never fires (endpoints >> eps). (siA not separately hashed;
  BarsCrossZ is the named uncrossed control and asserts coplanarClusterFaces==0.)
- ASCII-clean: no non-ASCII bytes in the 6-commit diff (src/test/docs). halfedge vocab
  correct (halfedge_.Pair/Start; no DCEL/twin misuse in the added code).
- GT7081 (ulimit -v 4000000): GREEN, 34.6s, no OOM (ArrangementBudget fail-closed path).
- ONE full filtered run Overlap3.*:Boolean2.*:CrossSection.* (ulimit -v 6000000):
  193 ran, 192 PASSED, 1 SKIPPED (Gate4c_HullMask_MustResolve - the pre-existing
  conditional skip in the untouched hull/sweep path), 0 FAILED. Matches the fence.

## VERDICT: PASS on correctness (push defensible); ONE disclosure/guard item

Every audit passes on correctness and the zero-oracle-wrong ABSOLUTE holds: no silent
wrong resolve found anywhere; all three mutations (recovery off, dedup off, +2eps
endpoint) fail CLOSED; EntangledBars(axis) resolves oracle-true with a CLEAN mesh;
the rotated resolve is the exact library boolean union (independently re-oracled);
openscad narrows 496->344 and stays fail-closed; v3 sweep byte-untouched; uncrossed
control bit-identical; full filtered fence green.

The one non-correctness finding (audit 2), which the lane and docs do NOT disclose:
the ROTATED synthetic pin's output carries a single DEGENERATE zero-area triangle -
the materialized R1 weld artifact (a 1-ULP-distinct wall vertex the emission weld
leaves unmerged), inert (vol/GWN/topology exact) but invisible to IsSelfIntersecting
and to the re-gate. It is confined to the synthetic 50-degree pin (no corpus carrier
hits it; the library's own union of the same bars is clean). "RESOLVES oracle-true"
overstates it: the resolve is geometrically correct but not mesh-clean. Before/with the
push, name the R1 artifact on the rotated pin and add a min-triangle-quality assertion
(and consider the once-only fold-vertex identity from audit 1, which would remove both
the positional dedup and this weld-reliance). None of this is a correctness regression.

