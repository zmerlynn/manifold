# V5 exact-arithmetic ADVERSARIAL VERIFICATION lane notebook

Branch: explore/sweep-plane-3d-v5 @ 8a802aff (canonical READ-ONLY; probe copy at
/tmp/v5-verify-probe, untouched except reads).
Started: 2026-07-13. Target: docs/ExactArrangement3D.md + generation notebook
.claude/lane-reports/v5-design-1783909556.md. Attack the two foundational
empirical claims: the exact-winding probe and the "emission-manufactured, not
intrinsic" pivot.

Independent tooling (NOT the generation lane's single segment-cast script):
  - SA  = solid-angle summation (Van Oosterom-Strackee), mpmath dps=100. NO ray
          casting; a genuinely different algorithm. Reports fractional residual
          (distance to nearest int) as robustness margin.
  - RAY = exact rational semi-infinite ray casting, 9 independent generic
          directions, forward interior crossings; unanimity across directions is
          an internal ray-invariance proof.
  eps per model = EpsilonFromScale(bBox.Scale(), 1000) (overlap3.cpp:1276), the
  same eps the pipeline uses. Scripts in scratchpad/: indep_probe.py,
  havoc_detail.py, havoc_detail2.py, thin_shell.py, bisect_sheets.py,
  carrier_sep.py, carrier_exact.py, gt_surface.py, gt7081_planes.py,
  attack4_cost.py.

================================================================================
## ATTACK 1 - INDEPENDENT PROBE REPRODUCTION : VERDICT CONFIRMED
================================================================================
Faithfulness gates independently reproduced: left/right both edge-manifold
bad=0, 6*vol sign=+1 (closed oriented 2-manifold, positive volume).

Two independent methods agree bit-for-bit with EACH OTHER and with the
generation lane at every query point (SA residuals ~1e-100 => unambiguous
integers; RAY unanimous across all 9 directions):

  mission_lump (13544.5,-382.156,2252.83): w_left=1 w_right=1 w_S=2  MATERIAL
  right_centroid                          : w_left=1 w_right=1 w_S=2  MATERIAL
  left_centroid (retained region)         : w_left=1 w_right=0 w_S=1  MATERIAL
  far_outside                             : w_left=0 w_right=0 w_S=0

INTERIORITY (exact): these are NOT boundary-adjacent points - they sit millions
of eps deep inside their regions (mission_lump 10.7M eps from both surfaces;
centroids 133M+ eps). "Strictly inside" is exact and robust.

RIGHT-ONLY WEDGE confirmed: a w_left=0 w_right=1 (w_S=1) sliver exists (found in
the L75/R7 column below). ADVERSARIAL near-boundary points: winding is a clean
integer down to ~1 eps from the surface (margins ~1e-100); only a point evaluated
EXACTLY on a face is ambiguous (margin 0.23, correctly flagged). A wrong w_S near
a boundary would require evaluating literally on the surface.

The probe reproduces. No disagreement to root-cause.

================================================================================
## ATTACK 2 - THE PIVOT, PER CARRIER : VERDICT HAVOC-SPECIFIC
================================================================================
Direct exact measurement of the minimum INPUT feature separation per carrier.

  carrier   eps        closest input feature                     verdict
  -------   --------   ---------------------------------------   ------------
  Havoc     2.25e-8    50 eps (edge-edge); verts 33413 eps       MANUFACTURED
  GT7863    4.51e-8    0 eps  (R vert #2 exactly on L tri         INTRINSIC
                        (63,60,61), confirmed exact coplanar+in)
  Offset1   4.51e-8    0 eps  (2 exactly-coincident verts,        INTRINSIC
                        dx=dy=dz=0)
  openscad  3.52e-10   0 eps  (71 exactly-coincident vert pairs)  INTRINSIC
  GT7081    4.51e-8    verts 61239 eps; 1.37M seamLen<=eps        INTRINSIC
                        near-tangent contacts (density lane) +    (surface
                        ~200 near-coplanar overlapping L-R faces  near-tangency)
                        (my gt7081_planes probe)

HAVOC pivot HOLDS (no sub-eps input feature: true surface min 50 eps > 1 eps, so
the 1.7-eps emitted twin is genuinely manufactured), BUT the doc OVERSTATES the
input cleanliness:
  - doc: "parallel sheets ~7e-4 apart (tens of thousands of eps)" = the
    VERTEX-vertex gap (33413 eps). The true SURFACE-surface min is 50.2 eps
    (edge-edge) - ~650x closer. Emission amplifies 50 eps -> 1.7 eps (~30x),
    plausible; the "tens of thousands of eps" headline is a vertex artifact.
  - doc/notebook: "right verts >= 23.45 eps from nearest left face PLANE" - plane
    distance is a weak proxy (a vert can be 23 eps from a far triangle's infinite
    plane); the true R-vert-to-L-FACE min is 29894 eps. The 23-eps number is not
    evidence of near-tangency.

The doc GENERALIZES ("the hard carrier's near-degeneracy is a property of the eps
sweep-emission MACHINERY, not of the input" - the strongest structural argument
for exact arrangement). REFUTED as a universal claim: 4 of 5 carriers have
INTRINSIC sub-eps / exactly-coincident input degeneracy. Only Havoc is
emission-manufactured. The framing must be scoped PER CARRIER.

Architectural consequence: "emission-manufactured => exact PLACEMENT might
suffice" is a HAVOC-ONLY hope. The four intrinsic carriers (exact coincident
verts / faces / near-tangent surface contacts) require the local arrangement -
which STRENGTHENS variant C over the placement-only path, but INVALIDATES the
doc's stated "sidesteps the amplification" rationale as the general argument.

================================================================================
## ATTACK 3 - INTEGER SHELLS : VERDICT CLEAN (substance) / CAVEATED (evidence)
================================================================================
The shells ARE clean integers with a well-defined {w_S>=1} boundary and the lump
sliver survives - confirmed by EXACT evaluation everywhere off-surface (margins
~1e-80..1e-100). The architectural conclusion holds.

BUT the doc's specific "2 -> 1 -> 0" pattern is LOCATION-DEPENDENT. Measured:
  - closest vertex pair L75/R7 column: 0 -> 1 -> 0  (right-only sliver, ~13000 eps)
  - closest surface approach L66 (62 eps): 2 -> 1 -> 0 (left-only shell ~75 eps)
  - mission-lump xy column: 0 -> 1 -> 2  (w_S=1 shell only ~905 eps thick)
Shell thicknesses span ~50 eps (edge-edge surface min) to ~13000 eps.

CRITICAL sampling caveat: at the mission-lump xy the w_S=1 shell is ~905 eps
thick and was SKIPPED by a 400-sample AND a 700-sample scan of the transition
(both reported a clean 0 -> 2 jump with NO w_S=1). Only bisecting the left vs
right sheet crossings separately (2251.83257 vs 2251.83255, |dz|=905 eps)
resolved it. A thin shell CAN hide between samples. Therefore "manifold-by-
construction locally" is justified ONLY by the exact boundary evaluation (winding
is a clean integer off-surface), NOT by sampling. The doc/notebook Step 4 present
SAMPLING (resolution unrecorded) as the evidence - that specific evidence is
unsound; the conclusion is right for the right reason (exactness), not the stated
reason (a stratification sample). The doc should state the sampling caveat and
lean the manifold claim on exact evaluation only.

================================================================================
## ATTACK 4 - ENGINEERING REALITIES
================================================================================
COST (exact w_S per query), Python+mpmath, Havoc soup 176 tris:
  SA  : 38.4 ms/query (26 q/s), ~0.22 ms/tri
  RAY : 15.7 ms/query (64 q/s), ~0.09 ms/tri
  cost is O(ntri) PER QUERY. Projected linearly to GT7081 (31360 tris):
  ~6.8 s (SA) / ~2.8 s (RAY) per single winding query.
  Absolute numbers are Python-specific; a C++ adaptive-orient3d kernel (FP filter,
  exact only when the sign is uncertain) would be ~100-1000x faster in the common
  case. The O(ntri) full-soup scaling is REAL and makes the doc's B-scales-with-
  input / C-bounds-to-submesh adjudication concrete. The tree has NO adaptive/
  exact kernel (doc inventory confirmed by my grep) so this is new surface area.

RAY-DEGENERACY FREQUENCY (fixed dir returns None = SoS/re-choose needed):
  generic dir (2,3,5), 3000 RANDOM pts : 0/3000 = 0.000%
  axis dir (1,0,0),   3000 RANDOM pts : 0/3000 = 0.000%
  generic dir, 180 ON-EDGE points     : 180/180 = 100.0%
  => For the intended usage (classify a cell by its INTERIOR point) degeneracy is
     effectively 0% - the SoS machinery is almost never invoked. It is the RULE
     (100%) only if you evaluate w_S exactly on a boundary, which the design
     correctly avoids by using cell interiors. The SoS/re-choose load is light for
     cell classification; R4's concern is well-placed but not a runtime hotspot.

================================================================================
## SETTLED / NOT RE-LITIGATED
================================================================================
The probe's w_S values are CONFIRMED by two independent algorithms. Faithfulness
gates reproduced. No modification to canonical repo (reads only); probe copy
untouched; scripts confined to scratchpad. This notebook is the only write.
