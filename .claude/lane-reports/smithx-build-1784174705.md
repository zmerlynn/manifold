# SMITH-EXTENSION lane 1 (construction): can P1/P2/P3 avoid the predicate AND close openscad?

Branch explore/sweep-plane-3d-v5, base CLEAN f53c9e3f (working tree has ~287 lines
of in-flight plumbing-lane drift in src/overlap3.cpp; I reset the copy's overlap3.cpp
to committed f53c9e3f and worked there consistently - stated per the task).  Copy
/tmp/smithx-a, from-scratch Release PAR=OFF.  Canonical READ-ONLY except this
notebook.  ASCII; file:line for code; [SPEC] marks speculation.

OWNER IDEAL: "extend Smith AND avoid the predicate AND close openscad."  Refuse
banned.  This lane hunted the real thing.  Deliverable is a POSITIVE (the predicate
IS avoidable) + a SHARPENED IMPOSSIBILITY (openscad closure is orthogonal to the
predicate; it is a representation wall no decision protocol reaches).

================================================================================
0. HEADLINE (the decoupling the "AND" hides)
================================================================================

The owner's three goals are NOT one problem - they are TWO ORTHOGONAL walls, and
the predicate stands between us and NEITHER:

  (A) "extend Smith + avoid the degree-9 predicate": ACHIEVABLE.  P1 (commit to
      once-only shared double positions) reduces every branch-node decision to
      level-0 input predicates + committed-double compares + SoS-by-input-index,
      which is JOINTLY REALIZABLE by construction (not merely consistent).  This
      SHARPENS act3-constructA: the degree-9 constructed-point orient2d it named as
      the E1 residue is AVOIDABLE, not just insufficient.

  (B) "close openscad": NOT a decision problem at all.  MEASURED this lane: the
      exact degree-9 predicate moves the triple SET (162 <-> 157) but NEVER the 21
      open edges (E1_OFF=21, default=21, E1_PURE=21).  The residue is the
      double-precision OUTPUT weld (SplitTouchingSheets, overlap3.cpp:77-216)
      collapsing sub-eps-thin true features born of near-parallel INPUT planes.
      That is a REPRESENTATION wall, downstream of every decision.

Because (B) is representation and P1/P2/P3 are all DECISION/perturbation protocols,
none of them closes openscad.  The predicate was never the thing blocking openscad
(act3-constructA/e1-close already showed it necessary-but-insufficient); this lane
adds that it is also AVOIDABLE, so openscad's wall is PURELY representation, fully
decoupled from the predicate question.

================================================================================
1. THE KEY STRUCTURAL FACT - verified, with one correction
================================================================================

Task claim: "a crossing {f,g} x {f,h} is a 2D event in EXACTLY ONE face's plane (f);
faces on g and h consume it only as a 1D position along their own seams."

CODE (EnumerateTriplePoints, overlap3.cpp:2511-2622): the loop is `for f: for k1:
for k2:` over seam pairs of face f.  A triple {f,g,h} is found in face f as
seam_g x seam_h, keyed by the SORTED plane triple {planeId[f],pg,ph} and constructed
ONCE via tripleTab (overlap3.cpp:2578-2596), position from Intersect3Planes on the
canonical reps (:2591-2593).  seamTriples[f][k1] and [f][k2] both receive the split
(:2618-2619).

CORRECTION (the fact is right only in the near-tangent special case): a GENERIC
triple point where all three planes meet carries THREE pairwise seams s_fg, s_fh,
s_gh.  In face f the pair (s_fg,s_fh) cross -> 2D event; in face g the pair
(s_fg,s_gh) cross -> ALSO a 2D event; in face h likewise.  So generically it is a
2D crossing in ALL THREE incident faces, and the enumeration finds it once in each
(deduped by the plane-triple key).  It is a 2D event in EXACTLY ONE face ONLY when
the third seam s_gh does not reach the point within its finite SEGMENT extent - i.e.
precisely the near-tangent regime.  This does not break the protocols' premise; it
sharpens WHERE the "1D consumption" lives: the once-only key + the junction registry
(BuildJunctionRegistry :2673, JunctionSplitsOnSegment :2718, threaded at :2857) is
what injects the shared point as a 1D split into any seam passing through it on ANY
face - the built mechanism for "consume as a 1D position."  VERIFIED present and
load-bearing on the resolving corpus.

So the structural substrate the protocols want (per-face 2D decisions Smith-safe;
per-seam 1D orders realizable; identity provenance) is real and already wired for
the enumerable triples.

================================================================================
2. THE JOINT-REALIZABILITY AUDIT (the owner's true question)
================================================================================

Owner's framing of the obstruction: "three faces' independently-perturbed views may
admit no single 3D perturbation."  This is EXACTLY Salesin-Stolfi-Guibas epsilon-
geometry's realizability gap (act3-lit-1784162607.md:232-239): a sign set can be
combinatorially derivable yet realized by NO nearby geometry.  Two facts settle it:

FACT 1 (SoS-by-input-index IS a single global 3D perturbation).  Edelsbrunner-Mucke
SoS keyed by input VERTEX index perturbs the input vertices by distinct
infinitesimals; every constructed point (a plane triple) inherits its perturbation
as a function of the input vertices it is built from.  So any decision resolved by
SoS-on-input-coords is jointly realizable BY CONSTRUCTION - there is a genuine
infinitesimal displacement of the input realizing all of them at once.  This is
boolean3's mechanism (sd-self-1784162422.md:154-159) and is already landed
(docs/Regularize3D.md stage-6 single global SoS).

FACT 2 (P1's commitment closes the near-degenerate gap SoS cannot).  SoS only fires
on EXACT ties; for near-degenerate (predicate within error of 0, sign genuinely
uncertain in double) SoS does nothing (act3-lit:226-231).  P1's move: COMMIT to one
set of once-only double positions as ground truth.  Then the per-face arrangement is
the EXACT arrangement of a set of genuine 3D double POINTS - trivially realizable by
those points - computed with CLASSIC orient2d on the exact axis-drop projection
(overlap3.cpp:2185-2187 notes the axis-drop is exact, no rounding), degree 2 on
plain doubles, not the degree-9 constructed-point form.  Ties among the committed
doubles are EXACT ties -> SoS-by-index -> Fact 1 -> jointly realizable.

CONSEQUENCE - the degree-9 predicate is AVOIDABLE.  The branch-node arrangement
needs only: (i) once-only shared positions (built), (ii) level-0 radial order from
input normals (attack22-c3-1784146879.md:65-70; f4-r4-1784128452.md:204-210 - 0
genuine ties), (iii) 1D seam order = arclength, monotone (act3-toy-1784163841.md:
109-121 - 0 flips corpus-wide), (iv) level-0 identity (nomerge-1784160552.md:101-110
- aliasGap=0), (v) SoS-by-index for exact ties.  NONE is the degree-9 form; (ii)-(iv)
are level-0 INPUT predicates, so P1 avoids even classic constructed-point orient2d.
This is act3-constructA Formulation-2's layered design; the NEW result here is that
it is not just CONSISTENT across faces but JOINTLY REALIZABLE (one input-index SoS
realizes every decision), which is the strong form the owner asked to exhibit.

So: P1 wins the "extend Smith + avoid the predicate" half OUTRIGHT.  P2 (global
lexicographic order) collapses into P1 - the induced per-carrier order already EQUALS
the geometric arclength order (act3-toy: 0 flips), so P2 adds no realizability the
committed-position order lacks and hits the same downstream wall.  P3 is deferred to
section 4.

================================================================================
3. THE RESIDUE IS REPRESENTATION, NOT DECISION (measured this lane)
================================================================================

openscad large dirty component, clean f53c9e3f (F4B_DUMP, ulimit -v 6000000):

  baseline (default, degree-9 exact extent):   distinct=157  openEdges=21
     (fan1=14 fan3=6 fan4=1; atTriple=4 oneTriple=10 noTriple=7)
  E1_OFF (pre-degree-9 rounded straddle+X-in-f): distinct=162  openEdges=21
  E1_PURE (ungated exact segment extent):        distinct=157  openEdges=21
  F4B_PERFACE (break once-only, per-face image):  perFace=1     openEdges=21

READ:
- The exact degree-9 predicate refines 162->157 (drops 5 rounded phantoms whose exact
  X falls outside seam g/h) but does NOT touch the 21 opens.  THREE predicate regimes,
  same 21.  The predicate is ORTHOGONAL to the residue.
- F4B_PERFACE (per-face back-projection instead of the once-only shared point) also
  stays 21: breaking cross-face CONSISTENCY does not move the count.  Consistency is
  NOT the residual wall either.
- Resolving carriers pass on this base (SelfIntersectA, GT7863, EntangledBars OK), so
  the measurement base is valid; on the resolving corpus the degree-9 refinement is a
  byte-clean no-op (e1-close-1784165580.md:85-88).

The wall is SplitTouchingSheets (overlap3.cpp:104-131): after BuildImpl welds output
vertices on a uniform eps grid (getVertIdx, :285-301), an output mesh edge with
fwd!=bwd half-edges = an open fan = fail closed (:108-109).  Prior exact-rational
reconstruction pinned WHY (f4-r4:114-131, f4-r5-1784304000.md:118-155): the 21 are
clusters of genuinely-distinct triple points spread ~1e-4 down to ~2*eps between
NEAR-PARALLEL walls, producing sub-eps-thin sub-faces (area to 1e-15, width ~1e-8
while eps ~3.5e-10).  The double eps-weld collapses them into balanced-but-non-
bounding fans.  The DECISIVE cross-check (e1-close:80-83, reconfirmed here as
E1_PURE=21): even splitting every exact crossing does not help because the
DOUBLE-PRECISION downstream cannot REPRESENT the sub-eps sub-cell - dropping the
resulting slivers makes it WORSE (f4-r5:132-143, 22->37).

This is sd-self's generic-triple-point obstruction realized at the OUTPUT: at a
generic triple point (3 sheets of ONE surface, un-perturbable) the collapsed 1-cycle
balances but does not bound (sd-self:202-233).

================================================================================
4. WHY EACH PROTOCOL MISSES IT (the sharpened impossibility)
================================================================================

The coupling that defeats all three, stated once: the OUTPUT is double-precision
vertices welded at eps (SplitTouchingSheets/getVertIdx), and at a GENERIC triple
point a sub-eps-thin sub-face - created by the near-parallelism of two INPUT planes,
a property of the input geometry, not of any decision - collapses at that weld into
a non-bounding fan.  This obstruction is DOWNSTREAM of every arrangement decision and
is a property of the INPUT PLANES' angles, so:

P1 (single-position) and P2 (global order) are DECISION protocols.  They make the
  arrangement consistent AND jointly realizable (section 2) - genuinely, and without
  the predicate - but that is necessary, not sufficient: the residue is not a wrong
  or inconsistent decision (F4B_PERFACE=21 proves consistency is not the wall; the
  three predicate regimes prove the sign is not the wall).  They bottom out unchanged
  at the eps-weld.  Measured: 21 -> 21 under every decision lever available.

P3 (local snap at branch nodes) is the only representation-flavored protocol, but it
  is at the WRONG LOCUS and forbidden from the right one:
  - WRONG LOCUS: the thin sub-faces' WIDTH is the distance between two near-parallel
    input-plane seams - set by the INPUT planes' angle, NOT by constructed-point
    spacing.  Snapping constructed triple points collapses sub-eps-CORNER clusters
    (the ~2*eps minority, e.g. f4-r4 edge 6 len 7.15e-10) but leaves the near-
    parallel-wall thin WEDGES (the dominant residue, spread ~1e-4 - corners far
    apart, only the width sub-eps).  It cannot remove a thin feature whose thinness
    lives in the input planes.
  - FORBIDDEN LOCUS: the fix that WOULD work is snapping the near-parallel INPUT
    planes exactly parallel/coincident so the thin cell never forms (f4-r5:180-186,
    the stage-5 wall-snap).  P3 explicitly excludes input, and act3-lit:326-333
    (friction #2) is why input snapping is contentious: for a self-overlap
    regularizer the coincident/near-coincident faces ARE the w_S signal; snapping
    them destroys it.  This is epsilon-geometry's OWN resolution (input perturbation,
    act3-lit:232-239) - the one known way through the realizability gap - and it is
    the one the non-input protocols cannot invoke.
  - r6's position-only quotient was already measured insufficient (sd-self:334-339);
    P3 differs (moves positions to unify DECISIONS, not identity) but the decisions
    are already unifiable without it (section 2), so it buys nothing the committed-
    position arrangement lacks, and the width residue survives it.

NO NEW PROTOCOL of this class escapes, because the class is "make the DECISIONS
right/consistent/realizable" and the wall is "the double OUTPUT cannot represent a
sub-eps INPUT-geometric feature."  The only routes that reach it (all outside the
class, all violating one contract term):
  (a) exact-rational OUTPUT arrangement + representability-aware rounding = 3D snap-
      rounding: provably manifold-preserving version is O(n^19), impractical
      (act3-lit:178-188); the deployable version (Valque-Lazard/CGAL 2025) gives NO
      manifold guarantee (act3-lit:199-208).  Violates the double-output-with-
      manifold-guarantee contract.
  (b) input near-parallel-wall snap: violates input-untouched + destroys w_S signal.
  (c) accept non-manifold sub-eps output: violates the gate.

================================================================================
5. VERDICT FOR THE OWNER
================================================================================

- EXTEND SMITH + AVOID THE PREDICATE: YES, achievable (P1).  The branch-node decision
  HAS a Smith-style form after all - not a "1D shadow of input coords" (sd-self was
  right there is none of THOSE), but a "1D/2D exact compare of COMMITTED once-only
  double positions + SoS-by-input-index," which is boolean3's exact-compare mechanism
  lifted to constructed points via once-only construction, and jointly realizable by
  construction.  The degree-9 constructed-point predicate is avoidable, and MEASURED
  orthogonal to openscad (162 vs 157 triples, always 21 opens).

- CLOSE OPENSCAD: NO decision/perturbation protocol closes it, because its residue is
  representation - sub-eps-thin sub-faces from near-parallel INPUT planes that the
  double-precision output weld collapses into non-bounding fans at generic triple
  points.  Confirmed orthogonal to the predicate (three regimes -> 21) and to
  consistency (per-face -> 21).  Closing it requires exact-rational output (no
  deployable manifold-preserving form exists) or input-plane snapping (violates
  input-untouched, destroys the w_S signal - epsilon-geometry's own and only
  resolution).

- THE "AND" IS THE TRAP.  The owner bundled predicate-avoidance and openscad-closure;
  they are orthogonal walls.  P1 gives the first for free and does nothing for the
  second; the second is not a predicate problem at all.  The one-predicate floor was
  never buying openscad closure either - it too stops at the same representation wall
  (e1-close:90-101).  So the honest terminal is: adopt P1 to make the pipeline
  predicate-free-and-realizable on the whole resolving corpus (a real Smith
  extension), and treat openscad's near-parallel-wall residue as a separate
  REPRESENTATION decision the owner must make (exact output, or input snap, or keep
  fail-closed) - it is not reachable by extending Smith.

[SPEC] The single un-ruled-out route to openscad-in-doubles is a boolean3-style
GRACEFUL collapse: when two sheets are within eps, cancel their emitted faces by
half-edge antisymmetry (as the EXACTLY-coplanar fold already does, f4-r4 A2).  The
fold handles exact coincidence; extending it to NEAR-parallel is snapping-coincident
= route (b).  Whether a cancellation that stays in doubles yet bounds at a GENERIC
triple point exists is, to my knowledge, open - but attack22-c3's winding-validated
completion probe (every candidate face distortion-banned or winding-rejected, 22->22)
is strong evidence against it.

================================================================================
RAILS
================================================================================
- Base: copy overlap3.cpp reset to committed f53c9e3f (working-tree drift excluded);
  from-scratch Release PAR=OFF.  Resolving carriers PASS (SelfIntersectA, GT7863,
  EntangledBars).  openscad fail-closed at 21, matching e1-close/f4 baselines.
- All measurements env-gated levers already in tree (F4B_DUMP, E1_OFF, E1_PURE,
  F4B_PERFACE); NO code added, NO predicate FORM touched, zero canonical edits.
- Zero-oracle-wrong preserved: no geometry emitted, no protocol shipped; the positive
  (P1 realizability) is an analysis + a citation to the built once-only/level-0
  substrate, the negative (openscad) is a measured orthogonality.
