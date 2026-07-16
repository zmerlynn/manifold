# CLOSURE BUILD: INPUT-EXACT basis unification (Regularize3D)

Branch explore/sweep-plane-3d-v5, base HEAD 7f8aec4e (src byte-identical to pushed
06784274). SOLE canonical editor. Owner law: refuse == counterexample; close or
per-edge proof. ASCII; FNV/memcmp never printed in committed code. Do NOT push
(adversarial verify follows).

NAME: INPUT-exact (planes composed exactly from input vertex doubles inside the
predicate), NOT in-exact. The lever is THE CAVEAT from maxdepth: the landed kernel
is exact-relative-to-ROUNDED-planes (faceN = cross of double edge vectors rounds at
construction). This lane tests the input-exact evaluation of the SAME homogeneous
form.

## THE THEOREM (dual theory lanes, both banked, triaged)
- Grammar CLOSED AT DEPTH 1: every arrangement vertex is {input vertex | plane
  triple f/\g/\h | plane-pair /\ point-pair-line}. No construction consumes a
  constructed point.
- Max predicate degree = 9 (plane-coeff basis) / 20 (input-coordinate basis,
  tight/homogeneous). Cramer triple point: X,Y,Z deg 7, W deg 6 in inputs.
- THE CAVEAT (the untested lever): landed CramerHPoint takes ROUNDED faceN doubles
  (PConst of the rounded plane coeff). Input-exact = compose n=cross, d=dot exactly
  from input vertex doubles INSIDE the predicate. Same one form, degree 20.
- depthk: compositional exact-integer evaluator (bi::Big), no monomial expansion,
  running-EBD filter is depth-general + sound; naive static bound unsound. But
  depthk deliberately used the ROUNDED plane-coeff basis; input-exact (input-coord
  basis) is the net-new instantiation this lane builds.

## THE WALLS TO MEASURE (each hypothesized a basis-mismatch artifact)
- f229/f257 grazes: foreign seam endpoints eps-injected onto faces they are only
  eps-near (graze lane: exact carrier membership drops them, build crosses, but
  emission opens 21->32).
- f467 chord: skew-near-crossing, two full-length seam edges whose carrier lines
  are 3D-skew; phantom crossing only after eps-snap.
- crowded cluster pairing: genuinely-distinct triple points spread ~1e-6..1e-5 on
  DIFFERENT plane triples within eps of a shared seam; per-face double arrangement
  not two-sided-consistent (stitch: 32->31, mechanism b survives).
- the 21 baseline opens: default RemoveOverlaps2D near-tangent residue (8
  default-default 1-ULP bit-diffs + the cluster).

## PLAN
1. INPUT-EXACT PREDICATES: compositional Big evaluator + input-exact CramerHPoint
   (n,d exact from input doubles) + orient2d deg-20 + running-EBD filter; existing
   callers untouched (bit-identity rail); validate vs independent oracle.
2. Re-land banked graze+stitch stack, route every exact decision input-exact.
3. Weld on flagged faces keys on exact identity.
4. MEASURE openscad opens trajectory per wall; resolve+oracle-grade OR honest
   per-edge terminal with theorem attached + feature INERT.

## LOG

### Baseline (build-ixbase, from-scratch Release PAR=OFF)
- openscad opens=21, distinct triples=156, b3faces=0, okfaces=574, fail-closed PASS.
- GT7081 RESOLVES (116s). GT7863/EntangledBars(+Rotated)/BarsCrossZ/SelfIntersectA/B/
  BridgedCaps/CoplanarFold family: 13 PASS.
- openscad eps (weld/merge radius) = 3.519e-10.

### DECISIVE OFFLINE MEASUREMENT (Python Fractions oracle over the 156 real triples)
Extracted every triple's sorted plane-triple + the 3 rep triangles' 9 input-vertex
doubles (IX_DUMP), computed the INPUT-EXACT position (exact rational Cramer over
exact planes) vs the pipeline's ROUNDED-plane position.
- rounded-vs-input-exact position gap: median 7.8e-11, MAX 4.999e-7 (8 triples > 1e-7).
  The plane-rounding gap reaches ~5e-7 = ~1400x eps in the wedge (condition-amplified).
- The dominant open cluster (11 triples around (-17.797,-2.698,-206.942), the named
  228/229 region incl {228,229,239/260/269/835} AND the off-229 {228,239,835},
  {228,257,260},{228,260,268/269}): EXACT min pairwise separation = 7.755e-9 =
  22x eps, 272840x the output ULP (2.84e-14). GENUINELY DISTINCT, ABOVE the merge
  radius, REPRESENTABLE in output doubles.
- f467 region (cluster {408,447,459/467},{408,459,467},{408,464,466/467}): exact
  min sep 3.849e-6 (well-separated, 10000x eps) - not sub-eps at all; the "skew
  crossing" is a straddle-decision question, not a representability one.
- ALL 11 near-coincident clusters: exact min sep from 3.8e-6 down to 7.8e-9, EVERY
  ONE > eps and >> ULP. NO cluster collapses to a genuine sub-eps/coincident point.

FINDING: the near-tangent clusters are NOT genuine sub-eps coincidences (the E1
representability reading). They are genuinely-distinct, representable points whose
ORDER the rounded-plane basis SCRAMBLES (plane-rounding gap 5e-7 >> the 7.8e-9
inter-point separation). This is the BASIS-MISMATCH signature the hypothesis
predicts: input-exact positions+decisions place and order them correctly, above the
merge radius. The stitch report's "spread ~1e-6..1e-5" was measured in the
contaminated rounded basis and overstated the true (input-exact) 7.8e-9 spread.
Proceeding to wire input-exact and measure openscad opens directly.

### WIRED: input-exact overlay (Big evaluator ported to overlap3.cpp; IX lever, default ON)
Ported the depthk compositional Big evaluator (validated vs Python Fractions: 159
openscad positions, 0.0 deviation). Wired input-exact into the wedge overlay:
- POSITIONS: CanonTriplePos routes every committed triple position through the
  input-exact Cramer (rep triangles -> exact n,d -> Cramer -> rounded to double).
  All keys (provOf / seamTriples / canon3) share ONE consistent input-exact pos.
- DECISIONS: buildH builds a parallel BigHPoint; orientH / properCrossH / the
  completeness check / ExtractCellsSymbolic cell-walk / signMul all use the
  input-exact BigOrient2D (degree-20-in-inputs; exact, no rounded-basis filter -
  which would certify a rounded sign that disagrees with input-exact truth).
- Completion crossings composed input-exact from planeTri rep faces.
Gated behind IX_OFF (default IX ON); byte-clean off the wedge population (provOf/
planeTri empty -> never reached; the resolving corpus has ZERO triples).

### OPENSCAD OPENS TRAJECTORY (per wall)
- IX_OFF (rounded overlay, reproduces stitch baseline): opens=31.
- IX ON (input-exact overlay):                          opens=21, b3faces=0.
  -> The input-exact overlay REMOVES the stitch/graze net-negative (31->21): the
     wedge faces no longer worsen openscad. The f229/f257 graze phantoms and f467
     skew-crossing DISSOLVE (b3faces=0, all 8 wedge faces build+complete); the
     crowded-cluster pairing divergence at the OVERLAY level is gone. The
     phantom/pairing walls are confirmed BASIS-MISMATCH artifacts.
- But 21 is the pure-default (EX2_OFF) baseline - openscad does NOT close.
  F4B_OPENDUMP: the 21 opens sit at the SAME crowded cluster (-17.797,-2.698,
  -206.942) in BOTH EX2_OFF and IX-on. They are the DEFAULT RemoveOverlaps2D
  residue, not overlay-induced.

### WHY 21 REMAINS (per-edge truth): overlay COVERAGE, not representability
- Only 8 faces are flagged for the overlay (the narrow gap-free-wedge detector);
  the crowded cluster spans dozens of faces (planes 228/229/239/257/260/268/269/
  835/...). The un-flagged cluster faces use double-precision RemoveOverlaps2D,
  which cannot subdivide the dense cluster two-sided-consistently.
- Forcing the overlay on ALL seamed faces (EX2_ALL): input-exact builds MORE faces
  than rounded (b3 11->6) and all 574 completing faces complete (completed=1); the
  6 residual b3 are the coplanar-cap seam-endpoint PROVENANCE gap (buildH, case (c)
  needs a transversal edge-neighbor) - a DISTINCT known issue, NOT cluster failure.
  But EX2_ALL opens 119 (IX) / 142 (rounded): overlay-on-all creates overlay-vs-
  default + skipped-b3 boundary holes. Not the path.

TERMINAL (refined): the residue is NOT the representation contract. The offline
oracle proved every cluster triple point is genuinely distinct, ABOVE the merge
radius (min 7.8e-9 = 22x eps), and representable in output doubles (272840x ULP);
the prior "spread down to 2*eps / eps-merge collapse" reading was the rounded-plane
error (1400x eps) masquerading as a spread. Input-exact predicates DISSOLVE the
overlay-level phantom/pairing walls but do NOT close openscad because the exact 2D
arrangement must COVER the whole cluster (replacing RemoveOverlaps2D on every
cluster face + two-sided-consistent overlay/default boundary), which is the
docs-E1 exact-rational arrangement ENGINE (priced ~600-1300 LOC net-new surface) -
NOT a deeper predicate (depth-1/degree-20 suffices; tripwire un-hit) and NOT a
representability limit. The input-exact kernel is the validated predicate basis
that engine needs.

### CANARIES (input-exact feature ON, build-ixbase)
- 15 resolving canaries PASS (GT7863, EntangledBars(+Rotated), BarsCrossZ,
  SelfIntersectA/B, BridgedCaps, CoplanarFold family, NegativeWinding, openscad
  fail-closed).
- GT7081 RESOLVES with IX on; F4B_TRIPLES=0 -> the input-exact path is NEVER
  reached off openscad (zero triples), so the feature is BYTE-CLEAN by
  construction on the whole resolving corpus. Byte-identity rail holds.
- EX2_ALL sanity: input-exact builds MORE than rounded (b3 11->6); the 6 residual
  are the coplanar-cap provenance gap, distinct from the cluster.

### DISPOSITION - RESIDUE (a FINDING); src REVERTED to byte-identical HEAD
openscad does not close (21 residue) -> per the fail-closed residue protocol the
feature ships INERT. src/overlap3.cpp REVERTED to byte-identical HEAD (git diff
empty). openscad stays fail-closed at its HEAD emission wall (21); resolving
carriers + GT7081 re-verified on the reverted build. Zero-oracle-wrong preserved;
NOTHING SHIPS (net LOC = 0). DID NOT push. The full input-exact stack (Big
evaluator + banked stitch overlay + input-exact positions/decisions) banked at
.claude/lane-reports/scratchpad/inexact-basis-1784205082.patch.

TRIPWIRE un-hit: the input-exact instantiation is the SAME homogeneous form (the
sign of an orientation determinant corrected by sign(prod W), on the compositional
exact-integer evaluator), degree 20 in inputs / degree 9 in plane coeffs - depth 1,
no nested construction, no third degree. The planes are composed exactly from input
doubles instead of the rounded faceN; existing rounded callers untouched.

### TERMINAL, folded to the ledger
Input-exact is the untested lever from THE CAVEAT (landed kernel = exact-relative-
to-rounded-planes). Measured: it DISSOLVES the graze/stitch overlay walls (the
f229/f257 grazes, f467 skew crossing, and the crowded-cluster pairing at the
overlay level: opens 31->21, b3faces=0, all 8 wedge faces build+complete) - the
phantom/pairing walls WERE basis-mismatch artifacts, as hypothesized. It also
REFUTES the "representation contract" reading of the E1 terminal: the offline
Fractions oracle proves every openscad cluster triple point is genuinely distinct,
above the merge radius (min 7.8e-9 = 22x eps), and representable in output doubles
(272840x ULP); the prior "spread down to 2*eps, eps-merge collapses" was the
rounded-plane error (1400x eps) masquerading as a spread. But input-exact does NOT
close openscad: the residue is the docs-E1 exact-rational 2D ARRANGEMENT ENGINE
COVERAGE - the exact arrangement must replace RemoveOverlaps2D on EVERY cluster
face (not just the 8 gap-free wedge faces) with a two-sided-consistent overlay/
default boundary; forcing it on all faces (EX2_ALL) drags in the coplanar-cap
provenance gap and creates boundary holes. That engine is the priced ~600-1300 LOC
net-new surface; the input-exact predicate basis built here (depth-1, degree-20,
validated) is the kernel it needs. Depth 1 remains sufficient; the tripwire is not
hit; the residue is an ENGINE gap, not a predicate or representability wall.

### CLOSE RAILS (reverted src, from-scratch Release build-ixfinal)
- src/overlap3.cpp byte-identical HEAD (git diff empty). Net LOC = 0 (INERT).
- Full Overlap3 suite: 32/32 PASSED (GT7081 resolves 116s; openscad fail-closed
  at 21; all resolving carriers + fold/near-coplanar/corpus pass-throughs).
- Feature-on rails (build-ixbase): 15 resolving canaries PASS + GT7081 resolves;
  byte-clean off openscad by construction (F4B_TRIPLES=0 -> input-exact path never
  entered off the wedge population). Fail-closed arms + contract-eps untouched.
- Banked patch: .claude/lane-reports/scratchpad/inexact-basis-1784205082.patch
  (Big evaluator + stitch overlay + input-exact positions/decisions, 2006 lines).
- DID NOT push. Only src/overlap3.cpp was ever edited (reverted); overlap3_sweep.cpp
  and overlap3_test.cpp were not touched by this lane.
