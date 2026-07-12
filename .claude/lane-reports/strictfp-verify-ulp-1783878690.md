# STRICT-FP ULP-RUN ADVERSARIAL VERIFICATION - lab notebook

Lane: attack the generation lane's central adjudication (strictfp-arc-1783872322.md):
"the run machinery does NOT dissolve; unbuilt runs persist at ulp scale; chain-plane
placement is still load-bearing; wide-run guard is unreachable-on-realistic-input but
kept as a fidelity backstop." Attack from both sides: is retention justified, and are
the ulp-regime claims TRUE?

Artifact: 94eb88df (landing) on explore/sweep-plane-3d-v3. Canonical READ-ONLY.
Workspace: /tmp/strictfp-verify-ulp (rsync copy), fresh Release build (vbuild/),
`-ffp-contract=off -fexcess-precision=standard -DMANIFOLD_PAR=-1`. All probes env-gated,
reverted (git checkout src/ at end). Machine: 7GB RAM total (limits big fixtures).

Drivers (compiled against the probed libmanifold.so):
- runhist_driver <case>: corpus/pair cases, RUNHIST + outcome.
- perp_driver: PerpFaces repro + slab dump.
- build_wide2 <N>: adversarial >eps run (staircase prism + offset box B), valid manifold.
Probes (all reverted): OV3_RUNHIST (run histogram from crits alone, pre-loop),
OV3_DUMPCRITS, OV3_SLABBOUND (revert ZipperEmit to slab bounds = the chain-plane
deletion), OV3_SPLITDBG (SplitTouchingSheets + BuildImpl weld/drop diagnostics),
OV3_NOCOVER (skip the coverage guard), OV3_NOSKIP (remove the ci!=li+1 cap skip),
OV3_STAGEDBG (SweepEmit stage boundaries + wideRun-true print).

====================================================================
## Q1 - ULP-RUN CHARACTERIZATION
====================================================================

### Structural fact CONFIRMED: unbuilt <=> exactly adjacent doubles (1 ulp).
RUNHIST checks `nextafter(crits[k],crits[k+1]) == crits[k+1]` for every unbuilt slab.
`nNotAdjacent1ulp == 0` on EVERY case (suite + full corpus + synthetics). So an unbuilt
slab is ALWAYS xHi == nextafter(xLo). A run of N unbuilt slabs = N+1 criticals at
consecutive representable doubles; run width ~= N ulp. This is exact, not empirical luck:
xMid=(lo+hi)/2 rounds onto a bound iff lo,hi are 1 ulp apart; >=2 ulp always gives a
strictly-interior representable midpoint -> builds.

### eps/ulp ratio (the "1000 ulp" estimate is scale-dependent and LARGER).
eps = ldexp((k+1)*kAlphaCoeff*kU, expBits), k=1000, kAlphaCoeff=12.37, kU=2^-53.
=> eps/ulp(at scale) = 6191 * 2^expBits: ~12382 for scale in [1,2), ~24764 for [2,4).
So a run > eps needs ~12000-25000 consecutive-double criticals, NOT ~1000.

### Run histogram (built from crits alone; independent of the section/budget path):
| case      | eps      | nSlabs | nUnbuilt | nNotAdj1ulp | nNonCanonCrit | maxRunSlabs | widestInteriorRun/eps |
|-----------|----------|--------|----------|-------------|---------------|-------------|-----------------------|
| Cray      | 4.7e+26  | 1      | 0        | 0           | 0             | 0           | 0                     |
| Offset2   | 4.5e-8   | 2387   | 1104     | 0           | 1101          | 8           | 1.6e-4                |
| Offset3   | 8.8e-11  | 73     | 41       | 0           | 38            | 5           | 2.0e-4                |
| Offset4   | 8.8e-11  | 118    | 65       | 0           | 57            | 7           | 2.8e-4                |
| GT7863    | 4.5e-8   | 890    | 149      | 0           | 149           | 9           | 7.3e-4                |
| Havoc     | 2.3e-8   | 532    | 123      | 0           | 121           | 6           | 4.8e-4                |
| Offset1   | 4.5e-8   | 9885   | 3839     | 0           | 3834          | 12          | 2.4e-4                |
| openscad  | 3.5e-10  | 8794   | 1028     | 0           | 1028          | 19          | 1.9e-4                |
| hull      | 3.5e-10  | 36316  | 11265    | 0           | 11264         | 16          | 5.7e-4                |
| self_A    | 2.7e-12  | 63089  | 9022     | 0           | 9018          | 6           | 2.4e-4                |
| self_B    | 2.7e-12  | 57394  | 13742    | 0           | 13739         | 7           | 5.7e-4                |
| **GT7081**| 4.5e-8   |231616  | 74199    | 0           | 74195         | **6862**    | **0.554**             |
| PerpFaces | 2.7e-12  | 4      | 1        | 0           | 1             | 1           | 4.0e-5                |

### THE HEADLINE: runs DO chain, and GT7081 (real corpus) reaches 0.554 eps.
maxRunSlabs=6862 consecutive 1-ulp slabs, widest interior run = 0.554 eps - more than
HALF of eps, on a REAL reported bug mesh. So the generation lane's "unreachable on
realistic input" is TOO STRONG in the strong sense: GT7081 gets halfway there; a
modestly denser corpus case (or a different eps) crosses.

### CONSTRUCTED a valid > eps run.
build_wide2 N (staircase prism whose right edge is a nextafter chain, macro-separated in
y so canonicalize keeps them distinct; offset box B). At scale 1.4 (eps/ulp=12382),
N=15000 gives an INTERIOR unbuilt run of 15000 1-ulp slabs, widthInteriorRun/eps=1.2114,
WIDERUN_REACHABLE=1, INPUT status=NoError (valid manifold input). So the guard's
precondition (crits[ri]-crits[li+1] > eps) is definitively constructible - NOT
structurally impossible. There is no bound making runWidth > eps impossible; the only
requirement is ~eps/ulp consecutive-double criticals, which real degenerate clusters
approach (GT7081) and adversarial input meets.

### BUT: the guard's FIRE POINT is structurally SHADOWED by the arrangement budget/OOM.
Firing the wide-run guard needs BOTH (run > eps) AND (non-empty cap = macro geometry
change across the run). The chain criticals are eps-clustered in x, so any NON-DEGENERATE
face that touches the chain verts must extend MACRO in x - into a built flank - where it
is sectioned. The section/cap arrangement then carries ~one piece per chain critical
(~12000-25000), and RemoveOverlaps2D on that dense arrangement OOMs (or with the retained
budget would ArrangementBudget) BEFORE the wide-run critical's cap is evaluated. Verified:
- build_wide2 N=12500 (1.0095 eps): STAGEDBG shows BuildSlabs OK + "calling EmitCaps"
  reached, then bad_alloc in the FIRST (exterior) cap arrangement of the dense flank -
  before the wide-run ci. So even reaching EmitCaps, the dense flank cap dies first.
- Three attempts to build a cheap "in-run donor" (thin wafer / bipyramid / (1,0)-fan)
  all produced status!=NoError: a face confined to the sub-eps x-run is DEGENERATE; a
  non-degenerate face necessarily spans into a flank. This Catch-22 is the structural
  reason.
- This exactly mirrors GT7081: its 0.554-eps run coexists with density that fatals
  ArrangementBudget (RUNHIST heavy-case run: GT7081 -> ArrangementBudget) - it never
  reaches caps.

Q1 VERDICT: The wide-run guard is NOT dead code. Its precondition (run > eps) is
reachable in principle (constructed as valid manifold) and half-reached on real input
(GT7081 0.554 eps). Its FIRE point is not reached by any fixture because the density that
makes a wide run drives the section/cap arrangement to OOM/ArrangementBudget FIRST. The
spec wording should change from "unreachable on realistic input / a pathological pile of
~1000 criticals" to: "the run precondition is constructible and GT7081 reaches 0.55 eps;
the guard's fire point is shadowed by the section/cap ArrangementBudget that the requisite
density triggers first (~12000-25000 consecutive-double criticals per eps, not ~1000)."

====================================================================
## Q2 - THE LATENT FRAGILITY: root-cause of the 1-ulp weld failure
====================================================================

PerpFaces slab structure (OV3_DUMPCRITS), eps=2.749438e-12:
  slab0 [0, 1-eps]                       width 1.0     BUILT   (a interior)
  slab1 [1-eps, 1-eps+1ulp]              width 1 ulp   UNBUILT (the run; a seam crit
                                                                sits 1 ulp above b's wall)
  slab2 [1-eps+1ulp, 1]                  width ~eps    BUILT
  slab3 [1, 2-eps]                       width 1.0     BUILT   (b interior)
The chain-plane rule binds slab2's LO strip at the run's canonical critical
crits[1]=1-eps (the cap plane), not at slab2.xLo=crits[2]=1-eps+1ulp. Displacement = 1 ulp.

Experiment (OV3_SPLITDBG; identical vertices, only ZipperEmit's x differs):
- CORRECT (chain-plane): 44 rawTris, 24 weldedVerts, 44 kept, 0 degen, **0 boundaryEdges,
  0 fanEdges** -> clean closed surface -> RESOLVE vol=2 (oracle).
- MUTATED (OV3_SLABBOUND, slab bounds): SAME 44 rawTris, **SAME 24 weldedVerts**, but
  **2 tris DEGENERATE -> dropped**, 42 kept, **4 boundaryEdges** -> SplitTouchingSheets
  returns false at `FAIL=count edge(6,11) fwd=0 bwd=1` (open boundary edge) ->
  NonManifoldEmission "unresolvable sheet contact".

The two degenerate triangles (pre-weld corners) are slab2's strip SLIVER at b's front-
bottom edge (y=0.3, z 0..1), where b's front wall gives a constant-(y,z) edge:
  tri1: (1-eps+1ulp,0.3,0) & (1,0.3,0) both weld to vert 11  (dx=eps-... <= eps)
  tri2: same collapse.
The strip sliver's x-width there = (slab2.xHi=1) - (lo placement):
  chain-plane lo=crits[1]=1-eps  -> width 2.74947e-12 > eps -> sliver SURVIVES the weld.
  slab-bound lo=crits[2]=1-eps+1ulp (1 ulp higher) -> width 2.74936e-12 <= eps -> the two
  corners weld together -> sliver COLLAPSES -> dropped -> holes.

DISCRIMINATION (mission's hypotheses):
- (a) different eps-cluster / min-index picking different reps: REFUTED. The welded vertex
  set is IDENTICAL (24 verts) in both modes; min-index weld is bit-identical. The 1 ulp
  does not change which verts exist.
- (b) adjacency/pairing change: it is a CONSEQUENCE (the dropped tris open the fan) but not
  the cause.
- (c) exact-duplicate/degenerate-tri drop eats a sliver and leaves an open edge: CONFIRMED.
  Precise mechanism: the slab-bound placement truncates a constant-(y,z) strip sliver's
  x-width from just-ABOVE eps to just-BELOW eps, so the assembly eps-weld collapses it,
  and the degenerate-tri drop opens 4 boundary holes the sheet splitter rejects.

Is it a shallow bug whose fix would let chain-plane DELETE? NO. The chain-plane placement
(emit the strip at the cap plane = where the pieces were paired) is the CORRECT placement;
the slab-bound placement is simply wrong (it disconnects the strip from the cap by the run
width and truncates the sliver under the weld radius). There is no shallow bug being
papered over - deletion is just incorrect emission. A more robust assembly that merged
sub-eps slivers into neighbours instead of dropping them would be a non-trivial change and
would still leave chain-plane as the right placement. So chain-plane retention is the
honest, DURABLE answer (recorded finding, not a fix for the next arc).

====================================================================
## Q3 - RETENTION ADJUDICATION
====================================================================

### CHAIN-PLANE RULE (loX/hiX placement) -> KEEP (confirmed load-bearing).
Re-ran the deletion experiment myself: OV3_SLABBOUND (= revert ZipperEmit to slab bounds).
- PerpFaces test (MUST-RESOLVE) reds: fatal=2 "unresolvable sheet contact" (root-caused in
  Q2). Clean run PASSES.
- Fast Overlap3 suite: clean 47 PASS; OV3_SLABBOUND 46 PASS / 1 FAIL (only PerpFaces).
So chain-plane is load-bearing and PerpFaces is its SOLE pin in the synthetic suite. KEEP.

### WIDE-RUN GUARD (ComputeCap wideRun, overlap3.cpp:871/952) -> KEEP + REWORD SPEC.
- Precondition reachable (Q1): constructed a valid > eps run (1.2114 eps); GT7081 real
  corpus reaches 0.554 eps. NOT provably dead - there is no structural bound on run width.
- Fire point shadowed (Q1): the density needed for a > eps run drives the section/cap
  arrangement to OOM/ArrangementBudget before the cap over the wide run is evaluated. So no
  fixture reaches the guard's FIRE, but not because runs can't be wide.
- It guards a genuine fidelity hole: with the interpolating chain-plane rule kept, a wide
  run with a non-empty cap that DID reach the caps stage would otherwise resolve silently
  oracle-wrong. Fidelity gate is absolute.
=> KEEP. NOT delete-with-record (not provably dead) and NOT downgrade-to-assert (it is a
   fail-closed on a reachable precondition, not an invariant). But the spec comments
   (ComputeCap :816, EmitCaps :947) must be reworded: replace "unreachable on realistic
   input / ~1000 consecutive adjacent-double criticals" with the accurate "run precondition
   is constructible; GT7081 reaches 0.55 eps; ~12000-25000 consecutive-double criticals per
   eps; the guard's fire point is shadowed by the section/cap ArrangementBudget the requisite
   density triggers first." Retention becomes PRINCIPLED (reachable precondition, budget-
   shadowed fire), not sentimental.

### COVERAGE GUARD (single-face SubEpsFeature, overlap3_sweep.cpp:243/251-281) -> KEEP.
REACHABLE (not dead): build_wide2 N=500 fires it directly - "face in unbuilt x-range with
area > perimeter*eps" - because the staircase side-quads live wholly in the run and are
macro in z. On suite+corpus it is DORMANT: no corpus case reaches it (GT7863/Havoc/Offset1/
openscad -> NonManifoldEmission; hull/self_A/self_B/GT7081 -> ArrangementBudget; nothing
reports "face in unbuilt x-range"). Reachable-but-dormant. KEEP.

### NON-CANONICAL CAP SKIP (ci != li+1, overlap3.cpp:945) -> KEEP.
- EXERCISED everywhere: nNonCanonCrit > 0 on every non-trivial case (Offset2 1101, GT7081
  74195). ci != li+1 genuinely occurs (multi-slab ulp runs).
- Removing it (OV3_NOSKIP): geometric correctness on Offset2/Offset3 and the fast suite is
  UNCHANGED (assembly de-dup swallows the doubled caps at ulp scale) - a nuance vs the
  generation lane's "doubled caps break it" framing. BUT it breaks the white-box pin
  Pin_OneArrangementPerCritical (doubled cap arrangements), and it is a major PERF/RESOURCE
  guard: on GT7081 it would compute 74195 EXTRA cap arrangements (each a 2D boolean). KEEP
  (perf/resource load-bearing + white-box pinned; correctness happens to be de-dup-covered
  under strict-FP).

====================================================================
## SUMMARY
====================================================================
No run machinery is structurally dead. Chain-plane: KEEP (load-bearing, root-caused as
correct placement, PerpFaces-pinned). Wide-run guard: KEEP (precondition reachable, fire
shadowed by budget) - REWORD the spec's "unreachable" claim. Coverage guard: KEEP
(reachable synthetically, dormant on corpus). Non-canonical skip: KEEP (exercised, perf/
resource load-bearing, white-box pinned).

The generation lane's adjudication SURVIVES on the retention decisions. The one correction:
its ulp-regime CLAIM that the wide-run guard is "unreachable on realistic input" is
overstated - the run precondition is constructible and GT7081 reaches 0.55 eps; the honest
statement is "budget-shadowed fire," which strengthens (not weakens) the retention case.
