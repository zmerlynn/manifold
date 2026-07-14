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
