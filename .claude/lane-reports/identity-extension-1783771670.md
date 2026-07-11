# 3D-IDENTITY EXTENSION arc

Started: 1783771670 (2026-07-11)
Branch: explore/sweep-plane-3d-v3 @ dc8bfb06
Suite baseline: 597 tests = 596 pass + 1 documented Gate4c skip

## Mission

At a critical, the L and R extensions of ONE 3D junction diverge 60-200 eps -
a forced-through WELD endpoint frozen at its section (y,z) versus its
track-extended twin on the other side - leaving a cap micro-edge no
subdivision scheme repairs. Identity must come from a SHARED CANONICAL SOURCE,
not independent interpolation/freezing.

Acceptance gate: flip Corpus_Havocglass8_Recorded and
Corpus_GenericTwin7863_Recorded to MustResolve, confirm RED, then make green
WITHOUT weakening any existing test (fence rule absolute).

Candidate shapes: (A) weld-endpoint identity, (B) cross-side unification at
the cap, (C) extend-through-the-junction (constructible crossing).

## Plan

1. Read grounding: SweepEmit3D.md, provenance-chains journal, coplanar-crucible
   Step 13, the code (overlap3.cpp SlabResolver/ComputeCap/EmitCaps,
   overlap3_sweep.cpp, boolean2_sweep.cpp).
2. Flip both gates to MustResolve, confirm RED, build.
3. Probe the exact failing Havocglass8 junction (instrument, measure the
   divergence) BEFORE designing.
4. Adjudicate A/B/C on the house simplicity bar; pick ONE.
5. Implement with a mutation-verified pin.
6. Full suite green; doc section; commit (no push).

## Log

### Step 1: grounding (DONE)

Read: SweepEmit3D.md (full), provenance-chains journal (full),
coplanar-crucible Step 4/13 (full), overlap3.h, overlap3.cpp
(SlabResolver, ComputeCap, BuildCapEdgeSet, EmitCaps, EmitStrips,
ZipperEmit, driver), overlap3_sweep.cpp (BuildSlabs, ComputeSectionSegment,
seamTracks), boolean2.h (SweepCapture, EdgeM), boolean2_sweep.cpp grep.

THE WALL, precisely (Havocglass8, cap x=13494.10429, eps~2.25e-8):
- The winding pass constructs T-junction section verts EXACTLY and never
  merges them, so ONE 3D junction J appears in a section as TWO twins ~62
  eps apart.
- SlabResolver clusters piece endpoints per 1x eps (line 546). 62-eps
  twins are NOT clustered -> resolved independently.
- One twin matches a track (class-i/ii/iii, InterpolateSafe lands ~on J at
  the critical); the OTHER twin matches NO track -> weld (constant (y,z),
  frozen at its section position, ~208 eps off J at the critical).
- The cap arrangement runs at capEps=8eps; the 208-eps micro-edge is kept
  -> sliver cap triangles (Havocglass8) or a strip-less cap edge
  (GenericTwin7863) -> SplitTouchingSheets fails closed.

KEY CODE LOCI:
- SlabResolver::Resolve (overlap3.cpp:572): face/seam/planar/weld cascade,
  eps_ ball match. Extend (:558): InterpolateSafe along track, or return pt
  (weld = frozen section coord).
- Resolver cluster (:533-554): per-eps-cluster, one track per ball.
- Tracks carry the true 3D endpoints (FaceTrack va/vb, SeamTrackEntry
  vA/vB, GroupedEdge a/b) - so a resolved track KNOWS its 3D geometry.
- Cap input built in BuildCapEdgeSet (:745): each piece's from/to extended
  independently by its slab's resolver; L verts then R verts into ONE
  arrangement at capEps=8eps.
- The 3D junctions ARE constructible: arr.verts (canonical 3D verts),
  arr.seams (3D segments), and track-crossings within a slab.

WHY THE THREE CANDIDATE SHAPES map to this code:
(A) weld-endpoint identity: a welded section vert that is the image of an
    arr.verts vertex / seam endpoint at the critical -> Extend TO that 3D
    point. Needs: associate the welded (y,z) at xMid with a 3D object whose
    (y,z) at the critical is the true junction. The track a NEIGHBOR twin
    found already names the 3D object.
(B) cross-side cap unification: unify the two cap-input verts that are one
    junction before MergeVerts (needs an identity signal spanning L/R).
(C) extend-through-the-junction: a track x track/seam crossing inside the
    slab is a constructible 3D point both sides compute identically.

Next: flip gates red-first (mandate 4), then PROBE the exact junction.

### Step 2: gates flipped RED (mandate 4, DONE)

Flipped the shared CorpusPairGate helper to MustResolve (dropped the
fatal-accepting arm). Build clean. Both RED:
- Corpus_Havocglass8_Recorded: fatal=3 "unresolvable sheet contact" (258ms)
- Corpus_GenericTwin7863_Recorded: fatal=3 "unresolvable sheet contact" (8.6s)
Iterate on Havocglass8. Red-first satisfied.

### Step 3: PROBE the exact junction (about to)

Reasoned geometry FIRST (to sharpen the probe):
- The section is taken at xMid; a 3D junction J lives at the critical c.
- Two 3D structures (seams/edges) that MEET at J are, at xMid (x != c), at
  two DIFFERENT (y,z) separated by (slope-difference)*(c - xMid) - the real
  62-eps section gap. NOT construction noise; slope-driven.
- Each such section vert SHOULD extend along its OWN track and BOTH tracks
  pass through J at c -> both land ~on J. Tracked-to-tracked across L/R
  agrees bitwise IFF they match the SAME seam (same vA,vB,c arithmetic).
- The bug is WELD-vs-TRACK: one twin finds its track (-> J), the other finds
  NO track -> welds (frozen at xMid, ends up ~= track-motion off J = 208eps).
- A T-junction (crossing of F,G section segments) lies on the F-G SEAM, so
  its (y,z) at xMid == seam's point at xMid within ~1eps. So a seam track
  SHOULD match it within eps_. If the weld happens, EITHER no seam F-G
  exists (coplanar/degenerate/seam-ended-at-c) OR it's a >2-segment /
  forced-through point.

PROBE QUESTION (decides A vs C): for the WELD twins, is there a seam/face
track NEAR (but beyond eps_) that its tracked twin matched, or genuinely NO
track? Instrumenting Resolve to print, per decision, nearest seam/face/
planar candidate distances (in eps units) + slab xMid.

### Step 3 RESULTS - the failing junction fully characterized

RESDUMP (Havocglass8): welds CONSISTENTLY have a seam/face track at ~8-50
eps (near-miss beyond eps_=1eps). Not trackless - near-miss.

CAPDUMP (failing cap xCap=13494.104291, R slab, 4 welds): the junction
carries section verts A/B/C within ~200 eps; B tracks, A/C weld.

SLABDUMP (definitive, failing R slab xMid=13494.10429, eps=2.25e-8):
The junction J is an arr.verts VERTEX V at the critical:
- seam [vA,vB]: vA=(13494.1043, -400.02159, 2214.03692), vB=(13497.04,...).
  vA.x == the critical (13494.1043). So V = vA is a 3D vertex AT the cap
  plane; its (y,z) = (-400.02159, 2214.03692) is the TRUE junction.
- face edge shared by f114,f116: [13552.74.., 13474.48..], crosses the
  section at B=(-400.021590465912, 2214.03692407955). Passes THROUGH V at c
  (both tracks meet at V).
- section verts at xMid:
    A=(-400.021590513751, 2214.03692389031)  WELD  (8.7 eps from B/face edge)
    B=(-400.021590465912, 2214.03692407955)  TRACK (ON the face edge, exact)
    C=(-400.021589364739, 2214.03692843553)  WELD  (200 eps from B; 53 eps
                                                     from the seam yzMid)

KEY STRUCTURE:
- ONE junction = arr.verts vertex V, V.x == critical c. Its yz is the exact
  shared 3D identity, computed identically by L and R (same arr.verts).
- The section (at xMid, ~7e-7 before c) spreads V into A,B,C over ~200 eps
  (spread/xdist ~6 = local slope). NO single track is near all three: A,B on
  the face edge; C on the seam; the two tracks MEET at V.
- So unification MUST be at V (the meeting point), not at any one track.
- The seam TERMINATES at V (vA.x==c): InterpolateSafe(vA,vB,c)=vA.yz EXACTLY
  (t=0). This is the exact 3D-identity target, radius-free.
- The welds are near-misses (8.7, 53 eps) beyond eps_=1; that's why they
  weld. Bumping eps_ uniformly to catch both needs ~60 eps (constant =
  suspect). The clean canonical is V.

ADJUDICATION (leaning A, realized as arr.verts@critical snap):
- (C) extend-through-junction: the junction IS constructible (V=arr.verts@c),
  but welds have no track/faces to compute a crossing from -> can't self-locate.
- (B) cap unification: correct locus (both sides present) but needs the same
  V identity signal; equivalent to A done at the cap.
- (A) weld-endpoint identity / 3D-identity: snap section verts to V=arr.verts
  vertex at c. V is the SHARED CANONICAL (bitwise L==R). This is exactly the
  spec's unimplemented mechanism-3 "3D-IDENTITY PREFERENCE".
  Association is the open question (welds are near-miss, not trackless).

### Step 3b - the junction is a NEAR-DEGENERATE arr.verts PAIR

arr.verts near J (yz within 5e-4):
  V1=(13494.1042908758, -400.021590424656, 2214.03692439446) dx=-36.4 eps
  V2=(13494.104292516,  -400.021590491914, 2214.03692375977) dx=+36.4 eps
V1.x, V2.x are the slab's two criticals (xLo, xHi), 72.8 eps apart in x,
SUB-EPS in yz. Canonicalize did NOT merge them (>eps in x). The failing cap
c = V1.x; V1 is THE junction there. The thin slab (73 eps wide) between the
degenerate pair spreads the junction into A/B/C (200-eps spread > slab
width). This is emission-time (slab built, no conflict) => IN scope.

SHARED CANONICAL = arr.verts vertex V with V.x == c (bitwise). L and R
compute V.yz identically (same arr.verts). Snapping the junction cluster to
V collapses the micro-edge.

SAFETY of a snap radius (fence rule): distinct arr.verts at EXACTLY x==c are
>eps apart in yz (Canonicalize 3D-merged anything closer). But the spread to
catch (C at 200 eps) forces R >> eps, and two distinct junctions at c could
be 10s-100s eps apart in dense geometry -> a constant R over-merges. Must
DERIVE R per-instance or guard by track-provenance.

### Step 4 - EMPIRICAL TEST: snap cap-input verts to arr.verts@c

Threading arr's per-critical vertex yz into BuildCapEdgeSet; snapping each
extended endpoint to the nearest arr.verts@c within R (env OVERLAP3_SNAPEPS,
eps units, 0=off). Sweep R to find gate-closure threshold + suite safety.

RESULT: x==c snap FIXES junction 1 (13494.10429) - confirmed via SPLITDBG,
the OVERLAP fail moved to a DIFFERENT junction 2 (12940.5137). Mechanism
validated for the AT-critical case.

Junction 2 (SPLITDBG OVERLAP k=4, cap edge dx=0, 55 eps): NOT fixed by
x==c snap even at R=5000. Root cause (CAPDUMP): vertsAtCap=0 - the junction
vertex V=(12940.5137,-470.842773,2319.81982) is at V.x=12940.5137 but the
cap critical xCap=12940.51369999 is ~4 eps LESS in x (a DISTINCT critical).
The tracks TERMINATE at V (V is their endpoint) but are evaluated at
xCap < V.x -> EXTRAPOLATE past V. L-side track lands 0.1 eps from V; R-side
track (steeper) extrapolates 54.8 eps past V. The x==c map has no vertex at
xCap -> no snap. Two flavors of ONE disease:
  J1: vertex AT the critical; section spreads it (incl. forced-through
      welds) -> proximity snap needed.
  J2: vertex NEAR the critical (4 eps in x); terminating tracks extrapolate
      past it and diverge -> snap to the terminus/nearby vertex.
UNIFIED FIX: snap extended cap verts to nearest arr.verts in 3D (compare
(c,y,z) to (V.x,V.y,V.z)), NOT requiring V.x==c. Handles welds (proximity),
terminating tracks (land ~on V), and extrapolating tracks (J2). Switching
to 3D-near snap; sweep R for closure + suite safety.

RESULT: naive 3D-near snap OVER-MERGES. R=30 -> tie; R>=100 -> UNBAL on a
MACRO edge (len=1.42, a real edge lost a triangle). Confirms the brief's
warning: a constant 3D radius merges distinct verts. Need PROVENANCE, not
bare proximity.

### Step 5 - provenance-grounded rules (two flavors, both -> shared arr.verts)

RULE 1 TERMINUS CLAMP (exact, radius-free): a resolved track [va,vb] has
endpoints that ARE arr.verts (face-edge / seam endpoints, shared). When
extending to c OUTSIDE [min,max](va.x,vb.x) (extrapolation past a terminus),
clamp to the nearer endpoint's yz. J2: R track [12942.04, V] with c<V.x
clamps to V exactly; L track [V, 12940.002] interpolates to ~V (0.1eps);
they merge. Fixes J2 by provenance (track terminus), no radius.
Concern: InterpolateSafe's intended unbuilt-run EXTRAPOLATION - does
clamping break it? Test the suite.

RULE 2 (J1 welds): forced-through welds AT a vertex-carrying critical need a
proximity snap to arr.verts@c. Safe BECAUSE restricted to x==c: distinct
arr.verts at exactly one critical are >eps apart in yz (Canonicalize
3D-merged closer). Add a clear-winner guard (2nd-nearest >> nearest).

Testing RULE 1 (clamp) first - env OVERLAP3_CLAMP.

### Step 6 - clamp too aggressive; safe snap partial; the residual classes

- CLAMP (universal terminus clamp) + snap: breaks a MACRO edge at 13878.9707
  (len 1.42) - clamping breaks legitimate wide-run extrapolation. Rejected.
- The junction-2 miss was NOT a terminus case: V.x is only 0.3 eps from
  xCap (a sub-eps critical pair), and my snap required v.x==xCap BITWISE.
  Relaxing to |v.x - c| <= eps fixes the association WITHOUT the clamp.

SNAP MECHANISM (|V.x-c|<=eps, nearest arr.vert, clear-winner guard
second>4*best, radius R):
- SAFETY (pivotal): full Overlap3 suite at R=30 = 47 pass, ONLY the 2 corpus
  gates fail. The snap does NOT regress any other fixture. Safe at R~30.
- Havocglass8: fixes junction 1 (13494) and junction 2 (12940.5137/-470.84);
  residual = TIE at nearby vertex (12940.5137/-464.58 at R=30..60,
  12941.28/-503.73 at R=100). R>=250 -> MACRO over-merge at 13878.
- GenericTwin7863: UNBAL 1F/0B at (-31165.1953, ...) persists (a strip-less
  cap edge, its own junction class).

DIAGNOSIS of the residual: arr.verts-only snap is INCOMPLETE. At R=250 the
13878 over-merge is a LEGIT point snapping to the nearest arr.vert when its
true canonical is a SEAM/EDGE crossing (not a vertex). The complete
canonical set at a cap plane = arr.verts@~c PLUS seam/edge crossings at c.
The dense Havocglass8 cluster has junctions of BOTH kinds within a few eps.
A full-canonical Voronoi-safe snap is the complete mechanism but larger.

Next: (a) check if the Havocglass8 TIEs are snap-created near-degenerate cap
triangles (droppable) or genuine; (b) decide complete-mechanism vs honest
partial per budget.

### Step 7 - the tie IS junction-2's pattern; the disease is uniform

CAPDUMP at the -464.58 tie: ONE arr.vert V=(12940.5137,-464.575409,...) at
dx=0.3eps. ALL pieces track to V (land 0.1-3.3 eps) EXCEPT R p5/p10 whose
track [a=(12942.04,-542.88,...), b=V] EXTRAPOLATES past V (c<V.x) landing
54.8 eps off - identical to junction 2. So the whole Havocglass8 dense
region is ONE disease: near-degenerate vertex clusters where the cap sits
0.3 eps from the true vertex V and the incident tracks either land on V or
extrapolate 54.8 eps past it (or weld). At R=30 the 54.8-eps ones miss; at
R=60 they snap but a near-tangency (gap 2.7e-14) remains -> snap creates a
degenerate cap triangle when the whole cluster collapses onto V.

The unified fix wants: (1) a per-junction (not constant) radius that catches
the spread when V is isolated but refuses in dense areas (Voronoi), and (2)
degenerate cap-triangle handling when a cluster collapses onto one V.
Trying: Voronoi-safe snap (snap P to nearest V@~c iff d1 < 0.5*d2, no global
R) + drop degenerate cap triangles.

### Step 8 - Voronoi snap is SAFE but insufficient; the wall is deeper

Voronoi-safe snap (snap P to nearest arr.vert V with |V.x-c|<=eps iff
best<0.5*second - a per-junction radius = half the local vertex spacing, no
global constant):
- SAFETY (strong): full Overlap3 suite at snapR=100000 = only the 2 corpus
  gates fail, zero other regressions. The ratio guard makes even an
  unbounded radius safe - it refuses to snap in dense areas.
- DROPDEGEN (drop near-degenerate cap slivers): BREAKS closure - the slivers
  are LOAD-BEARING (dropping them opens macro holes, len 70). Rejected.

The residuals after the safe snap are DEEPER than twin-divergence:
1. GenericTwin7863: the snap FIXES the primary micro-edge (the steep-track
   0.5-eps-critical-offset amplified to 62 eps - CAPDUMP confirmed:
   xCap 0.5 eps from V.x, track slope ~125), but EXPOSES a MACRO UNBAL
   (len 5.2e-4, a hole). Collapsing a junction cluster onto V creates
   VANISHING/degenerate pieces and can FLIP arrangement topology (moving a
   vert 62 eps crosses neighboring edges), unbalancing the emitted fan.
2. Havocglass8 13527: a DENSE cluster - FOUR arr.verts within ~7000 eps, TWO
   distinct vertices V1,V3 ~300 eps apart BOTH on the (sub-eps-adjacent) cap
   plane, plus a near-x-degenerate track. The tracks already land on their
   correct vertices (no divergence to fix); the defect is the emission
   itself around two near-coincident-but-distinct vertices. Snapping to
   "the" vertex is ill-defined (which of V1/V3?).

WHY THIS IS THE RESEARCH-GRADE WALL (adjudication A/B/C, final):
- (A) 3D-identity snap: the RIGHT shape for the PRIMARY divergence (safe,
  fixes it), but input-snapping alone cannot close the gates - it needs
  COORDINATED cap+strip RE-EMISSION around a collapsed junction (a vanishing
  piece must not leave a hole; a moved vert must not flip topology), which is
  a re-mesh, not a snap.
- (B) cap unification: same input-snapping limits; also can't distinguish
  "spread of one junction" (unify) from "distinct near-degenerate vertices"
  (V1,V3 - do NOT unify) at the 100-300 eps noise scale = the core hard
  problem, unchanged.
- (C) extend-through-junction: the dense clusters have no single constructible
  crossing; welds/degenerate-tracks self-locate nothing.

### Step 9 - FINAL STATE (committed honest partial)

Reverted src/overlap3.cpp to clean HEAD (all TEMP DEBUG + exploratory snap
threading removed - verified no remnants). Restored the corpus gates to
recorded-contract with the comment refined to the near-degenerate cluster
residual. Added the "3D-IDENTITY EXTENSION" section to docs/SweepEmit3D.md
(instrumented finding, Voronoi-safe mechanism, A/B/C adjudication, honest
residual). Full suite green: 597 = 596 pass + 1 skip (Gate4c) = baseline.
git diff: docs/SweepEmit3D.md + test/overlap3_test.cpp only.

WALL (precise, for the successor): every corpus failing junction is a
near-degenerate arr.verts cluster where the cap critical sits sub-eps from
the true 3D vertex and a steep track amplifies that offset to 60-200 eps. A
Voronoi-safe snap to the shared vertex (per-junction radius = half local
spacing, no constant) resolves the PRIMARY divergence and is provably SAFE
(zero non-corpus regressions at unbounded radius), but cannot close the gates
alone: (1) collapsing a cluster onto one vertex creates vanishing pieces and
FLIPS cap-arrangement topology globally (GenericTwin: fixes the micro-edge,
opens a macro hole); (2) dense clusters carry multiple distinct vertices
within a few hundred eps, indistinguishable from one spread junction at the
noise scale (Havocglass8 13527). Closing needs coordinated cap+strip
RE-EMISSION around collapsed junctions - a local re-mesh, not an input snap -
plus eps-scale spread-vs-distinct disambiguation (the near-coplanar arc's
problem). Not landed: the snap trades sliver-overlaps for holes and its pin
can't be met via the gates. Stopped at the iteration budget.

---

DECISION (honest partial, per fence rule + budget): the snap is SAFE and
correctly targets the primary divergence, but does NOT close either gate and
its secondary effects (holes/ties from cluster collapse) are not a clean net
improvement, so I do NOT land the exploratory snap (simplicity + the mechanism
pin can't be satisfied via the gates). I restore the gates to recorded-
contract (the pre-arc honest state, NOT a weakening), remove all TEMP DEBUG,
and document the refined wall. The predecessor's "twin-position divergence" is
REFINED: closing needs coordinated re-emission around collapsed/near-
degenerate junction clusters, and disambiguation of spread-vs-distinct at the
noise scale - beyond a bounded snap. Capped at the iteration budget.
