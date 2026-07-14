# reg3d PLANE-REP PROBE lane (owner-commissioned: "probe it first, cost estimate")

Branch explore/sweep-plane-3d-v5 (HEAD 77812a76 confirmed). PROBE ONLY: offline
Python exact rationals in scratchpad, NO production code changes (src byte-clean).
ASCII; halfedge vocab; zero-oracle-wrong absolute.

## Mission (two deliverables)
1. THE PROBE. Represent the failing dirty region symbolically: input face planes
   exact (rationals from double verts), constructed verts = PLANE TRIPLES
   (references, never coords), edges = plane pairs. Build local arrangement
   symbolically; classify cells by exact w_S. DECISIVE: does the symbolic
   arrangement PAIR THE FANS (8-edge hole closes with correct winding, nothing
   rounded mid-pipeline)? Grade vs exact oracle (w_S ground truth; a+b volume).
   Then OUTPUT-EXTRACTION: round symbolic verts to double at the end (naive
   first) - does the rounded mesh stay manifold/oracle-true, or does the sliver's
   thinness re-emerge (which verts collide/degenerate; would a snap-consistency
   rule fix it - CHARACTERIZE, do not build)?
   Carriers: GT7863 (8-edge near-tangent sliver hole), PokedCube (w=-1/w=1
   double-sheet at degenerate near-triple-point).
2. THE COST ESTIMATE (measured, not theory): (a) predicate FORMS needed
   symbolically (which determinants over plane coeffs, degree in input coords);
   (b) what of the landed resolver survives vs is rewritten (stage map);
   (c) LOC estimate vs campaign record (B core ~1.9k over stages 1-3, exact
   accumulator 218); (d) output-extraction difficulty class; (e) carrier coverage;
   (f) perf order-of-magnitude (symbolic vs current, measured in probe).

## Step 0: ORIENT - grounding read DONE
Sources read at HEAD: reg3d-s7 (plane-rep probe SPEC, the pre-agreed escalation),
reg3d-c1a (proven-open anatomies), reg3d-census + Gap ledger, docs/Regularize3D.md,
src/overlap3.cpp (kernel + winding + seam + emit).

Key grounding facts (measured by prior lanes, not re-derived):
- ALL emission carriers fail SplitTouchingSheets at B1 UNBALANCED = OPEN BOUNDARY
  (a hole / dropped faces), NOT a touching/doubled sheet. Zero surviving weld
  twins eps..16eps; zero non-manifold edges. Fans need MORE faces, not merged.
- GT7863/GT7081 (1a): near-coplanar sliver. Two near-parallel seamed faces, gap
  ABOVE the exact-coplanar fold threshold but the sub-face is thin; per-face
  emission in DOUBLES drops a boundary sub-face -> 8-edge hole. x span ~0.76 at
  x~-31165.
- PokedCube (1b): +++ corner warped to (-1,-1,-1) everts, genuine w_S=-1 region
  (GWN oracle = -1). w=-1|w=1 junction is a double sheet the mult-1 per-face
  emission opens; 8 emitted tris -> 3 closed boundary loops incl 2 isolated
  islands whose connecting neighbours were dropped at the degenerate spike.
- Current kernel = ONE predicate FORM: orient3d (4-pt determinant, DEGREE-3 in
  input coords), filter (double) -> ExactOrient3D (integer mantissa*2^E, ~100 LOC
  + accumulator 218 lines total incl SoS) -> SoSOrient3D (symbolic perturbation).
  EdgePiercesTri, WindingAt all reduce to orient3d calls. Constructed positions
  (SegPlanePoint) are DOUBLE, built once.

PokedCube exact reconstruction (from src): Cube({1,1,1},true) = [-.5,.5]^3, 8
verts, 12 tris (impl.cpp:104); vertex (.5,.5,.5) warped to (-1,-1,-1). All
coords exact rationals (halves + integers). GT7863: OBJ pair (left 110v/212f,
right 126v/244f) at x~-31165; dirty region = near-coplanar overlap of the two.

## Plan
- Build probe_core.py: exact rational plane kernel (plane-from-3-pts, vertex =
  plane-triple via Cramer, exact w_S ray oracle, exact volume). Sanity-test.
- PokedCube probe (clean rationals, decisive, the harder shape): build local
  arrangement symbolically, classify cells, test fan pairing + output rounding.
- GT7863 probe: extract near-coplanar sliver face pair from OBJ, local
  arrangement, test closure + output collapse.
- Measure predicate forms/degrees actually evaluated; write cost table.

## Step 1: exact kernel + PokedCube winding field - DONE
probe_core.py: exact rational plane kernel (plane-from-3-pts, vertex=plane-triple
via Cramer, exact w_S = signed segment-crossing oracle matching WindingAt, exact
mesh volume). Sanity: unit cube vol/center-winding correct; plane-triple exact.

PokedCube exact reconstruction (poked.py): 12 tris, 9 distinct planes (the two
everted cone faces 10,11 have irrational-free integer planes (1,1,-2) and
(-2,1,1); apex (-1,-1,-1)).
- signed mesh volume = 0 (everted cone exactly cancels the body's +1 in the
  SIGNED integral; NOT the {w>=1} volume).
- exact winding histogram (15^3 generic-rational grid): w=1: 182, w=-1: 142,
  w=0: 3051, none: 0. So {w_S>=1} is NONEMPTY (a w=1 pocket, centroid
  ~(0.03,-0.28,0.06)); w=-1 spike pocket centroid ~(-0.51,-0.23,-0.47). Body
  CENTER (0,0,0) is w=0 (the everted cone passes through it, subtracting the
  body's +1). Reconciles with prior lanes' comparable w=1/w=-1 regions.
- KEY: winding is a consistent integer field (jumps +1 across each oriented face,
  -n side; verified no contradictions on 3375 exact samples). So d{w_S>=1} is a
  closed 2-manifold BY THEOREM (a0-verify-witness: topology is a pure function of
  the soup). "Does exact close the fans" = YES by winding-consistency, provided
  the arrangement is BUILT exactly at the concurrences.

## Step 2: PokedCube radial closure - PARTIAL (degenerate spike not fully built)
radial.py: exact radial-assembly closure check (incident half-faces around an
edge, exact angular sort, wedge windings, retained-face pairing) + arrangement
vertex enumeration as plane triples. poked_close.py driver.
- Pierce points at the spike = EXACT rationals -1/14, 1/14 (macroscopic, ~1/7
  apart, NOT sub-eps) - matches c1a's verts 6,8,10. So PokedCube's pierce points
  do NOT round-collapse; its wall is NOT sliver precision.
- Pairwise seam extraction found only 3 clean seams: the everted cone (10,11)
  piercing the y=-.5 / z=-.5 body faces. The cone-tri x body-tri pairs at the
  APEX region return no-seg (my clip drops the degenerate incidences where the
  intersection line runs through a shared apex vertex). These degenerate
  concurrences ARE the crux and need SoS to build - reimplementing that in Python
  = rebuilding the exact resolver, out of probe budget. Closure OK=26/OPEN=0 on
  the NON-degenerate edges (consistent with the theorem) but the spike edges were
  not fully exercised. HONEST: PokedCube's exact closure rests on the
  winding-consistency theorem + confirmed field, not a full Python arrangement.
- REFRAME (measured): PokedCube's failure is NOT pierce-point rounding (pierce pts
  are macroscopic rationals). It is the NEGATIVE-WINDING / double-sheet emission
  (s7 step1: mult-1 rule "keep iff w_above==0" correctly drops w_above=-1 cells;
  the {w>=1} boundary near the w=-1|w=1 region needs the coordinated double-sheet
  the per-face walk opens). This is a CLASSIFICATION/EMISSION completion,
  ORTHOGONAL to vertex precision - CONFIRMS s7's "plane-rep does NOT buy negative
  winding." Pivoting the constructive closure probe to GT7863 (2-face sliver, the
  mission's LEAD carrier, buildable exactly).

## Step 3: GT7863 sliver + OUTPUT-EXTRACTION - DONE (decisive)
Parsed OBJ pair to exact rationals via the doubles the strings denote (as the
resolver sees them). left 110v/212f, right 126v/244f, NO shared exact verts. The
two twins share the x=-31165.2676 wall EXACTLY coplanar (393+ gap=0 dih=0
cross-pairs) and INTERPENETRATE (266 transversal LxR crossings). Cross-component
(non-fusion pass-through) dominates; the within-component dirty region is a
near-coplanar sliver.
- eps = kPrecision*Scale = 1e-12*31165 = 3.12e-8. ULP at magnitude 31165 = 2^-38
  = 3.64e-12. eps/ULP = 8567. THE KEY RATIO.
- OUTPUT EXTRACTION (scale-faithful rounding, part A): an exact sub-eps strip
  (two long edges d apart) with d in [1e-11, 1e-8] rounds to DISTINCT doubles
  (rounds_identical=False for all d > ~1 ULP) but the rounded separation is < eps
  -> the ASSEMBLY WELD (uniform eps-grid merge, R1) collapses it, NOT the
  rounding. Only d < 1 ULP (~3.6e-12) round-collapses. d > eps survives.
- REAL GT7863 (part B): seam-endpoint-to-mesh-edge distances < 50 eps: 522 total.
  Bucketed: 515 EXACT-INCIDENT (< 1 ULP, essentially 0 = degenerate triple/
  edge-on-edge incidences, dominated by the shared coplanar wall between the
  interpenetrating twins); 0 round-collapse; 7 genuine thin strips in (1 ULP,
  eps): widths 3.2e-9..8.2e-9 = 0.10..0.26 eps = 887..2259 ULPs, all on hub face
  225 (198x225,153x225,135x225,133x225) = THE SLIVER. 0 survive.
- VERDICT (output extraction): the 8-edge hole's sliver strips round to DISTINCT
  doubles (rounding lossless, ~900-2260 ULPs apart) but fall WITHIN the weld
  radius -> the uniform eps-weld merges them -> collapse. "Thinness re-emerges at
  output" = TRUE, culprit = the WELD (R1), not naive rounding. Difficulty class:
  NEEDS SNAP-CONSISTENCY RULE (symbolic union-find keyed by the exact arrangement
  replacing the positional eps-weld; two exact-distinct verts within eps merged
  consistently across all incident faces). This is exactly stage-5's near-coplanar
  merge decision RELOCATED to output. Since a sub-eps strip is below
  double-manifold representability, the correct output MERGES the two sheets =
  what stage-5 already does on input. => for the SLIVER class, plane-rep +
  output-snap == stage-5 relocated; the wall is RELOCATED, not eliminated
  (confirms s7). Plane-rep's genuine GT7863 payoff is narrower: the 515 exact
  incidences get bit-identical canonical triple identity (no per-face-pair
  double-construction disagreement -> no dropped faces), i.e. the R1/once-only /
  triple-point payoff, NOT the sliver.
