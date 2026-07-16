# HOMOG-DESIGN lane: homogeneous unification of the exact predicate (Regularize3D)

Branch explore/sweep-plane-3d-v5, task HEAD bc350f70. Canonical READ-ONLY except
this notebook; ALL work in rsync copy /tmp/homog (from-scratch cmake Release
PAR=OFF -j8). Standalone differential harnesses in /tmp/homog/homog_harness/.
ASCII; FNV/memcmp never printed.

OWNER COMMISSION: "any way we can pivot the existing predicate to use this new
form so it's the same code?" GENERALIZE the one blessed exact predicate FORM
(Orient3DExactSign) into a HOMOGENEOUS orientation, so the openscad closure's one
new decision (orient2d of {f,g,h}-constructed crossing points in a face) is the
SAME code, with today's input-point orient3d as the w==1 special case.

Ground truth read (in order): the kernel (src/overlap3.cpp sos:: namespace,
Orient3DFilterSign, Orient3DExactSign, Intersect3Planes, ExactOrient2DDrop,
ExactSegProperCross, the tripwire comment); f4-design-b (the priced
plane-based-symbolic-rational route, ~150-250 LOC); f4-r6 + attack22-c3 (the
terminal: orient2d-on-constructed is the ONE tripwire crossing; radial order and
identity are proven level-0 FREE); nomerge (identity/aliasing level-0, exact-zero
= genuine coincidence).

## VERDICT (one line each)
- w==1 orient3d instantiation is BIT-IDENTICAL to today's Orient3DExactSign
  (harness: 0 diff / 1.5e7 incl. 4.19M genuine ties; + real suite 13/13 identical).
- general orient2d (degree 9) matches exact rationals over 1.3M cases (1e6 random
  + 3e5 wedge/collinear adversarial): ZERO disagreements, ZERO wrong filter certs.
- accumulator: 12 limbs used on real data (112 suffices); degree-9 worst-case needs
  ~320 - widen the constant to keep totality-by-construction, zero runtime cost.
- filter (construction-aware, C~20 derived / 16.2 measured / 32 certified): 0%
  escalation on 1e6 generic, 22.68% on the near-parallel wedge family.
- SoS stays input-point-scoped; the new site's exact-zero = genuine coincidence
  (level-0 incidence, nomerge), no perturbation.
- LOC ~230-350 for the predicate surface; ZERO vendored code; ONE accumulator.

## 0. THE MATH (the unification)

Affine orient of d+1 points p_i in R^d = sign of the (d+1)x(d+1) determinant of
rows (coord_i..., 1). A HOMOGENEOUS point (X,Y,Z,W) denotes affine (X/W,Y/W,Z/W).
Multiply row i of the affine orient determinant by W_i:

  affine_orient = det[ homogeneous rows (X_i,Y_i,Z_i,W_i) ] / (prod_i W_i)

so, since sign(a/b) = sign(a)*sign(b),

  SIGN(orient) = sign( det[homogeneous points] ) * sign( prod_i W_i ).        (*)

That is the ONE FORM. The determinant is the same Leibniz expansion the kernel
already sums on its adaptive integer accumulator; the only new piece is the
W-product sign correction.

- ORIENT3D (d=3), the trivial instantiation: an INPUT double point (a.x,a.y,a.z)
  is the intersection of the three trivial axis planes {x=a.x, y=a.y, z=a.z};
  its homogeneous weight is W=1 (denominator exactly 1). Then prod W_i = 1 (sign
  +1) and the 4x4 determinant has an all-ones W-column -> it is EXACTLY today's
  orient3d determinant. So (*) with (d=3, all W==1) == Orient3DExactSign, and if
  the W-column is skipped as the literal ones-column (it is the constant 1), the
  arithmetic is bit-for-bit today's 24-permutation, <=3-factor accumulation.

- ORIENT2D (d=2), the general instantiation: a crossing point in face F is the
  triple point F /\ g /\ h. Its homogeneous coords come from Cramer over the
  three plane coefficient vectors; X,Y,Z,W are each a 3x3 determinant of plane
  coefficients (degree 3). Drop F's dominant normal axis -> 2D homogeneous point
  (say (X,Y,W)); orient2d of three of them = the 3x3 determinant of degree-3
  entries = degree 9, corrected by sign(W_0 W_1 W_2). At W==1 (the trivial 2D
  point) it reduces to the exact 2D orientation ExactOrient2DDrop computes today.

## THE DESIGN (HomogOrientSign, one form, two instantiations)

The ONE form: `sign( det[homogeneous rows] ) * sign( prod_i W_i )`, the determinant
summed on the existing adaptive integer accumulator (sos::SumSign), filter-first
with a PER-DEGREE compile-time error constant.

Structure for bit-identity (harness/h1.cpp): a template flag TrivialW.
  template<bool TrivialW> HomogOrient3DSign(pts[4][4]):
    for each of the 24 permutations, accumulate the product of the row factors;
    `if constexpr (TrivialW) { if (col==3) continue; }`  -> the W column is the
    LITERAL ones column, skipped exactly as ExactOrient3D does; and
    `if constexpr (TrivialW) return detSign;`            -> prod W == 1, sign +1.
  So TrivialW=true compiles to the SAME 24-permutation loop with the SAME
  col==3 skip and the SAME <=3-factor products as today's ExactOrient3D.
  TrivialW=false includes the W column (4-factor terms) and multiplies in the
  weight-product sign.
  HONESTY NUANCE: the prototype's MulMag4/SumSign4 form the <=3-factor product via
  a schoolbook 4-limb multiply (same INTEGER, mag[3]==0), not the historical
  specialized cnt==1/2/3 __int128 branches - so the trivial path is proven
  OUTPUT-bit-identical (harness: 0 diff / 1.5e7) but not instruction-identical. For
  the real landing, TEMPLATE SumSign/MulMag on the mag width so degree-3 reuses the
  EXACT existing width-3 MulMag/SumSign (literally the same instructions) and
  degree-9 uses the width-8 path - ONE routine, the trivial call unchanged down to
  the instruction. That is the bit-for-bit form the owner asked for; the prototype
  proves the arithmetic agrees.

## PROTOTYPE 1 - w==1 BIT-IDENTITY (DONE, PASS)

Standalone TU (harness/h1.cpp, h1b.cpp): verbatim-copied sos::ExactOrient3D +
Orient3DFilterSign (bc350f70) vs HomogOrient3DSign<true> + HomogFilterSign<true>.
  - h1: N=1e7, mix of wide-exponent random + adversarial near-coplanar (ULP nudges
    along the plane normal). filterDisagree=0, exactDisagree=0. ~6.1M certified by
    the filter, ~3.9M filter-uncertain (all exact-nonzero here).
  - h1b: N=5e6 FORCED ties/degeneracies (repeated points, shared-coord planes,
    exact in-plane integer combos). exactZero(genuine coplanar tie)=4.19M, and
    filterDisagree=0, exactDisagree=0 -> the exact-zero path also agrees bit-for-bit.
VERDICT: the w==1 instantiation is BIT-IDENTICAL to today's predicate on both the
filter verdict and the exact sign (incl. the genuine-tie path). The owner's
reduction holds by construction and is confirmed over 1.5e7 inputs.

### ECOLOGICAL CONFIRMATION (the real pipeline, not just random inputs)
Wired the COPY's production `Orient3DExactSign` to route through
`sos::HomogOrient3DSign<true>` (built pts[4][4] with W=1), rebuilt, ran the
overlap3 suite. Baseline and homogeneous-form: IDENTICAL 13/13 PASS -
Orient3DExactSign_PropertyPin, SelfIntersectA/B resolves, all ExactZeroTie
resolves (Constructed/GT7863/EntangledBars/BarsCrossZ), the CoplanarFold family,
and openscad still fail-closed unchanged. So the w==1 instantiation is a true
drop-in for the blessed predicate on real geometry, byte-unchanged. (Copy-only
edit; canonical untouched.)

## PROTOTYPE 2 - GENERAL orient2d on {F,g,h}-constructed points (DONE, PASS)

harness/h2_core.h + h2.cpp + gen_cases.py. Each crossing point P_i = F /\ g_i /\ h_i
built by Cramer over the plane coefficients (matches Intersect3Planes bit-for-bit:
c12=cross(g,h), c20=cross(h,F), c01=cross(F,g); W=F.c12; X=dF*c12+dg*c20+dh*c01).
Drop F's dominant normal axis -> 2D homogeneous (A,B,W); orient2d = sign of the 3x3
determinant * sign(W0 W1 W2).

(a) DEGREE ANALYSIS. Generators = the plane coefficient doubles (12 per point:
    nF,dF shared + ng,dg,nh,dh). Each Cramer coord X,Y,Z,W is a sum of 6 monomials,
    each a product of 3 generators -> DEGREE 3. The orient2d 3x3 determinant is a
    sum of 6 products of 3 such degree-3 entries -> DEGREE 9 (fully expanded: up to
    6*6^3 = 1296 monomials, each a product of 9 generator mantissas). The prototype
    builds exactly this monomial set (Poly arithmetic) and feeds it to the SAME
    accumulator (GeneralSumSign == sos::SumSign widened to the monomial length).
(b) DENOMINATOR-SIGN CORRECTION. Handled: the exact sign = sign(det) *
    sign(W0)*sign(W1)*sign(W2), each W-sign an independent GeneralSumSign of that
    point's degree-3 W polynomial. (In the padded-4x4 alternative the weight
    product is W0^2 W1 W2 -> a global sign flip only; handedness is moot, exactly as
    ExactOrient2DDrop notes today.)
(c) DIFFERENTIAL ORACLE. gen_cases.py computes the exact sign in Python Fractions
    over the bit-exact (.hex()) plane coefficients; h2 recomputes via the integer
    accumulator and diffs. FULL RUN: 1,300,000 cases = 1e6 random plane triples +
    3e5 adversarial (near-parallel wedge pairs at eps 1e-13..1e-4, near-collinear
    crossing families, both combined). RESULT: EXACT vs ORACLE disagreements = 0;
    FILTER wrong certifications = 0. (13k and 160k interim runs identical.)

## ACCUMULATOR-WIDTH VERDICT

Measured (h2 full 1.3M, same-scale plane coefficients incl. the wedge family):
max active accumulator width = 12 limbs -> the current 112-limb ceiling SUFFICES
with large headroom on real data (the per-call adaptive nLimbs already keeps it
tiny). BUT the PROVABLE totality bound differs by degree:
  - degree 3 (today): 192 mag bits + esum spread 6291 -> ~104 limbs -> 112 chosen.
  - degree 9 (general orient2d): 9*53 = 477 mag bits + worst-case esum spread
    9*(971-(-1126)) = 18873 bits -> (477+18873)/64 + 3 ~= 305 limbs.
So under ADVERSARIAL full-range exponents the degree-9 form's totality-by-
construction needs kAccumLimbs ~ 320 (up from 112). RECOMMENDATION: widen the
compile-time ceiling to ~320 to keep the "never window-fail, total by construction"
guarantee for the general form (cost: a fixed ~2.5KB stack array; ZERO runtime
change - the per-call active width stays ~12 on mesh data). Widen ONLY the constant;
the arithmetic is unchanged. (Alternatively document a per-degree exponent-spread
precondition; widening is cleaner and preserves the kernel's totality pride.)

## FILTER for the general form (derivation + escalation)

KEY FINDING: a naive "C * u * permanent-of-the-final-determinant" static filter is
UNSOUND for constructed points. In the wedge regime the constructed coordinates
X,Y,W are all ~eps (near-parallel planes: cross(g,h)~eps), formed by O(1)
cancellations, so their ABSOLUTE error is ~u while their magnitude is ~eps; the
final-determinant permanent collapses to ~eps^3 while the error floor stays ~u ->
ratio blows up (measured 1.4e15 with that permanent). The certified filter must be
CONSTRUCTION-AWARE.
  DERIVATION. Propagate the forward error through the construction: use the all-abs
  companion Pdet (every input |.|, every subtraction -> addition), a monotone upper
  bound on every intermediate AND on the sum of the 1296 monomials. Then
  |fl(det) - det| <= gamma_K * Pdet. Op-count along the deepest monomial path:
  cross (mult=1, sub=2) -> *d (3) -> sum-of-3 (5) [Cramer coord, ~gamma_5];
  det 2x2 of coords (each ~gamma_5): product ~gamma_11, sub ~gamma_12, *A ~gamma_18,
  sum-of-3 ~gamma_20. So K ~ 20 -> conservative static constant C = 32
  (errbound = 32*u*Pdet). MEASURED effective C = max(det.err_running / (u*Pdet))
  = 16.18 over the full 1.3M (< the derived 20 < the chosen 32) -> the static C=32
  is CERTIFIED-conservative on the tested family. The prototype filters with the
  exact running bound (EBD) - tighter and provably certified - and reports both.
  ESCALATION (full 1.3M): RANDOM 0 / 1,000,000 = 0.0000% (generic planes always
  certify - exact never fires off the near-tangent residue, matching the campaign),
  WEDGE 68,045 / 300,000 = 22.68% (the honest ill-conditioning cost; of those,
  W-sign-uncertain in only 261 = 0.38% - the escalation is dominated by the
  near-collinear det, not the near-parallel W).

## SoS: stays INPUT-POINT-SCOPED (why the new site needs none)

SoS perturbs by unique INPUT VERTEX INDICES to hand a consistent virtual sign to a
degenerate input configuration (vertex exactly on a face). The constructed crossing
points have NO input vertex index to perturb by - they are symbolic (a plane
triple). And an exact-zero orient2d of three crossing points is NOT a resolvable
degeneracy to perturb away: it means the three points are exactly collinear in the
face = the three seam lines are concurrent (a genuine >3-plane coincidence) or two
crossing points coincide (the same triple point). Per nomerge, that is exactly the
LEVEL-0 aliasing/incidence event: exact-zero = genuine coincidence, merged/kept by
GENERATOR IDENTITY (concurrency / same-triple), never by perturbation. So the new
site needs no SoS; SoS remains scoped to the input-point orient3d instantiation
(where the K==0 group of Orient3DSoS already IS this exact sign).

## TRIPWIRE COMMENT (draft for the real landing)

  // TRIPWIRE (owner contract, docs/Regularize3D.md open list): this is the ONE
  // blessed exact predicate FORM - the sign of a homogeneous orientation
  // determinant corrected by sign(prod W_i), summed on the ONE adaptive-width
  // integer accumulator (sos::SumSign). It has exactly TWO instantiations:
  //  (1) INPUT-POINT ORIENT3D (degree 3): four points with W==1 (an input double
  //      point == the intersection of its three trivial axis planes, denominator
  //      1). The W column is the literal ones column and the weight-product sign
  //      is +1, so this instantiation is BIT-IDENTICAL to the historical
  //      Orient3DExactSign (same 24-permutation, <=3-factor accumulation) - all
  //      existing callers unchanged.
  //  (2) CONSTRUCTED-POINT ORIENT2D (degree 9): three in-face crossing points,
  //      each the Cramer intersection of a plane triple {F,g,h}. Filter-first via
  //      the degree-9 construction-aware static bound (C*u*Pdet, all-abs Pdet); on
  //      a filter-0 the exact homogeneous sign fires; an exact zero is a GENUINE
  //      coincidence routed to the level-0 incidence path (nomerge), never a
  //      perturbation.
  // A THIRD instantiation of a NEW DEGREE (or any new constructed-point form) is an
  // OWNER DECISION - it widens the accumulator's proven totality bound and needs a
  // new per-degree filter constant; never add one silently. Vendoring rule intact:
  // if an exact primitive OUTSIDE this one form is ever needed, VENDOR Shewchuk's
  // public-domain predicates.c - do NOT rebuild expansion arithmetic piecemeal.
  // The certified fast path never touches the exact kernel.

## LOC (honest estimate for the real landing at post-axis-drop HEAD)

The PREDICATE surface only (the arrangement CONSUMER = f4-design-b's B1, ~120-180
LOC, is separate and already priced):
  - generalize sos:: kernel (DecomposeH; MulMag to <=9 factors -> 8-limb mag;
    SumSign to wider terms; widen kAccumLimbs to ~320; the TrivialW template)  ~90-130
  - HomogOrient2DSign exact (Cramer-coord Poly build + det expansion -> SumSign,
    or a leaner factored small-bignum eval)                                     ~70-120
  - the degree-9 construction-aware filter (all-abs Pdet + static C, or EBD)     ~60-90
  - SoS scoping (no code; the incidence route already exists) + tripwire text    ~5
  TOTAL predicate ~ 230-350 LOC. Matches f4-design-b's "B2 ~150-250 one new FORM"
  plus the accumulator widening + the certified constructed-point filter it did not
  cost. ZERO vendored code; ONE accumulator; the trivial path byte-unchanged.
