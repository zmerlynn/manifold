# Simplify A+B distillation - ADVERSARIAL VERIFICATION lane

Branch explore/sweep-plane-3d-v5 (HEAD 5c9ceae4). Canonical repo READ-ONLY;
this notebook is the only write, no commit. Target: docs/Regularize3D.md (the
distilled current design) + the DCEL->halfedge sweep of ExactArrangement3D.md +
generation notebook .claude/lane-reports/simplify-ab-1783947110.md. Default
skeptical: the artifact is a distillation, so the job is catching what
distillation LOST or dissolved wrongly, not re-litigating the design.

VERDICT: NEED-CHANGE. The design SURVIVES (no BREAK - both dissolutions are
defensible, the fail-closed posture holds, every fidelity spot-check traces). Two
doc-level distillation gaps need fixing before owner-ready: (1) the localizer
dissolution dropped a MEASURED cost residue the cited notebook flagged as new
surface area; (2) R1's "re-gate catches the failure" overclaims - the re-gate is
structurally blind to the exact weld-merge failure it is the backstop for, and
R2's harmless-direction mechanism is inverted.


================================================================================
## ATTACK 1a - SELECTIVE WELD DISSOLUTION : SURVIVE (sharpen R1)
================================================================================
Charge: the doc dissolves selective weld because B replaces whole components; did
it dissolve a real mechanism into a phrase?

NOT wrongly dissolved. The doc SPLITS it honestly: the EXTERIOR patch-to-exterior
stitch is dissolved (sound - B resolves the WHOLE component, so there is no
un-resolved exterior of the same shell to sew a local patch back into; the sweep
needed the stitch because it resolved only a local patch), and the B-INTERNAL
weld is explicitly KEPT as the load-bearing risk R1. The mechanism is not a
phrase; it is a named risk.

But R1 UNDERSELLS two things the cited wall-A / strict-FP arcs already proved:

(1) THE RE-GATE CANNOT CATCH THE FAILURE IT BACKSTOPS. R1 says "the re-gate
    catches the failure but does not repair it." I walked this against
    IsSelfIntersecting (properties.cpp:138). Its skip: `for i,j: if
    distance2(triVerts0[i], triVerts1[j]) <= (2*eps)^2 return;` - skip any
    tri-pair with ANY two verts within 2*eps. When the assembly weld MERGES two
    exact-distinct B verts (the R1 failure), the two triangles straddling the
    merge now share that vertex position (distance 0 <= (2eps)^2), so the pair is
    SKIPPED. A self-FOLD created by the merge, between triangles touching the
    merged vert, is INVISIBLE to IsSelfIntersecting - the merge manufactures
    exactly the shared-vertex config the gate skips. Only the validity half
    (IsManifold/Is2Manifold) backstops, and only for the NON-manifold weld
    outcome (a pinch/tear), not a fold. So a weld-merge that FOLDS rather than
    tears passes both halves = SILENT wrong, not caught. R1's "catches" is an
    overclaim; the honest statement is "validity catches non-manifold weld
    outcomes; IsSelfIntersecting is blind to self-folds at the merge."

(2) THE BOUNDED-WELD ESCAPE WAS ALREADY CLOSED. R1 offers "a bounded/exempt
    internal weld may still be needed." wall-A C3 (walla-arc, PHASE B) probed
    exactly this - BuildImpl weld radius x1.8/2.5/4.0 - and it KILLED at every
    multiplier: merging the ~1.7-eps twins "trades the micro-edge fan for
    degenerate-collapse defects (vanishing pieces / topology flips)". So a
    CONSTANT-radius bounded weld is a two-sided trap (too small re-manufactures
    the twin; too big collapses slivers to holes) - the doc's escape hatch is
    narrower than "may still be needed" implies.

FAIR TO B: B's genuine defense is real and should be R1's rebuttal - (a) the
ONCE-ONLY construction rule (each intersection built once, referenced everywhere)
directly attacks the sweep's twin-manufacturing, which came from TWO independent
projected derivations (a cap image via L, a cap image via R) of one junction; B
never makes two rounded images of one point. (b) radial-not-projection means B
has no O(seam^2) projected twins. But once-only only prevents twins of the SAME
point; it does NOT prevent two GENUINELY DISTINCT arrangement points from rounding
within eps and getting merged by the uniform assembly weld. That residual is open
and, per (1)+(2), neither cheaply weldable nor reliably re-gated.


================================================================================
## ATTACK 1b - LOCALIZER DISSOLUTION : NEED-CHANGE (cost hidden)
================================================================================
Charge: localizer dissolved into "the component is the scope"; GT7081 is ONE
large component, so B's arrangement+winding runs whole-component. Does the doc
state the O(ntri) cost honestly or hide it in the dissolution?

HIDDEN. The doc does NOT state B's whole-component winding cost anywhere - not in
the B-mechanism WINDING bullet, not in the open list, not in the risks. The only
scale mention (R2) is about the GATE (IsSelfIntersecting), not B's winding.

Yet the CITED evidence measured it. v5-verify-probe ATTACK 4: "cost is O(ntri)
PER QUERY. Projected to GT7081 (31360 tris): ~6.8s (SA) / ~2.8s (RAY) per single
winding query ... The O(ntri) full-soup scaling is REAL ... The tree has NO
adaptive/exact kernel so this is NEW SURFACE AREA." The generation notebook
(simplify-ab, LOCALIZER row) kept the residue explicitly: "the O(ntri)/query cost
survives, scoped to component." The distilled doc dropped BOTH.

Fair nuance (B is better than naive O(ntri)/cell): the Winding03 discipline the
doc cites amortizes - one seed ray per connected cell-component (O(ntri) each) +
integer-delta flood (cheap), so real cost is O(#cell-components * ntri) + flood,
not O(#cells * ntri). But on a large single component (GT7081) #cell-components
can be large and each seed is still O(ntri); GT7081 being ONE component is
precisely the localizer's lost-scoping case. The doc should carry this as the
localizer dissolution's honest cost residue + name the new-surface-area kernel
(no adaptive orient3d in the tree today - grep confirmed: zero orient2d/orient3d
in src/, include/).


================================================================================
## ATTACK 2 - R2 GATE-NOTION MISMATCH : NEED-CHANGE (sharpen; ii inverted)
================================================================================
Read IsSelfIntersecting fully (properties.cpp:138-191). ep=2*epsilon_; skips any
tri-pair with two verts within 2*eps; on a distance-0 contact it tries four
+-eps normal nudges and returns NON-intersecting if ANY nudge separates. So the
gate is SYSTEMATICALLY CLEAN-BIASED (both the (eps,2eps) vert-proximity annulus
AND the eps-nudge leniency push toward "clean"). That asymmetry is the key the
doc's R2 does not state.

(i) GATE CLEAN but B DIRTY - THE SHARP ONE, correctly the dangerous direction.
    Two non-adjacent triangles that genuinely cross but have two verts in the
    (eps, 2eps) band (distinct enough not to be merged by the eps weld, close
    enough to trip the 2*eps skip), OR a crossing shallower than one eps nudge:
    gate skips/clears it -> component EARLY-EXITS as clean -> unregularized
    self-overlap passes through SILENTLY (worse than fail-closed). B's arrangement
    (which skips only TOPOLOGICALLY-adjacent pairs) would find it. Severity:
    real silent-wrong path, but NARROW - requires a near-degenerate config in a
    thin band just outside the eps merge radius. The doc names it ("a missed
    dirty component - silent"); good, but should say WHY (the gate is clean-biased
    by construction, not incidentally).

(ii) GATE DIRTY but B FINDS NOTHING - essentially a NON-ISSUE, and the doc's
    stated mechanism is BACKWARDS. The doc's parenthetical "(2*eps relaxation
    flags near-misses)" is inverted: the 2*eps nudge SUPPRESSES near-misses
    (returns non-intersecting when an eps nudge separates), it does not flag them.
    For the gate to flag a pair it must be a robust >eps-deep crossing between
    non-adjacent triangles - which B ALSO finds (B's shared-vertex skip set is a
    SUBSET of the gate's distance skip set: any topologically-shared vertex is at
    distance 0, so the gate skips it too; therefore B-skipped ⊆ gate-skipped, and
    nothing the gate flags is skipped by B). The only genuine (ii) is a legitimate
    touching contact the gate flags dirty and B re-emits clean = wasteful, not
    wrong. So (ii) is harmless; the doc should fix the inverted mechanism.

COUPLING (unstated, matters): R2(i) and R1 are the SAME blind spot. The weld-merge
in R1 creates the shared-vertex config the gate skips (attack 1a-(1)), so the
re-gate that "catches" R1's failure IS the gate whose clean-bias R2 flags. A
weld-merge fold lands in exactly R2(i)'s blind band. The doc treats R1 and R2 as
independent; they are one hole.


================================================================================
## ATTACK 3 - DISTILLATION FIDELITY : SURVIVE (one minor citation nit)
================================================================================
Spot-checked 8 "proven by pointer" claims against cited notebooks - all trace:

  1. w_S=2 at Havoc lump, two independent algos bit-for-bit  -> v5-verify-probe
     ATTACK 1 (SA solid-angle + 9-dir exact ray, residuals ~1e-100). TRUE.
  2. Static filter sound over ~1.6B predicates, zero certified flips -> v5b-s12:
     722M orient3d (289M vf + 433M ee) + 867M orient2d on GT7081 = ~1.59B; "ZERO
     sign-flips at ratio>1". TRUE (the ~1.6B is an accurate GT7081 aggregate; the
     "-s12" citation is a touch loose since the number is a roll-up, but honest).
  3. Zero genuine multi-sheet edges on the corpus -> v5b-s34 Finding 1: Havoc
     297 / siA 26078 / siB 26078 lines, ALL 2-distinct-plane, 0 with >2 planes;
     seam counts cross-checked bit-for-bit vs the pipeline. TRUE.
  4. Edge-edge z-order 100% certified, fallback only on genuine ties -> v5b-r4:
     D3 z-order n=5376/31656 cert=100%, 0 fallback; siA enumeration 100% certified
     across 94582. TRUE (fallback lives in D1, not the z-order decision).
  5. w=2 alternation residual absorbed by threshold read (measured non-alt SSEE)
     -> v5b-r3 item(i): measured raw SSEE (not SESE); threshold read discards the
     2 spurious interior sheets, thresholds soup winding directly. TRUE.
  6. No ambiguity witness exists; topology a pure function of the soup ->
     a0-verify-witness: "w_S is a PURE FUNCTION of S", "witness CANNOT exist in
     the union regime", constructive witness-attempt-that-fails. TRUE.
  7. Candidate A resolves the decomposable corpus; fold dropped (owner) ->
     v5b-fold + commit 71770901 "owner decision - drop the A fold, keep the gate";
     GT7081 "~1e-10 rel, NOT bit-exact" (bit-exact retracted). TRUE.
  8. Zero silent wrong-resolves on the union corpus -> ExactArrangement3D.md
     section "Zero silent wrong-resolves" (line 1160) exists. TRUE.

SYMBOL MAP (grep, all exist): Canonicalize (overlap3.cpp:218), Winding03/_
(boolean3.cpp:388/462), Shadows (shared.h:121), IsSelfIntersecting
(properties.cpp:138), IsManifold/Is2Manifold gate (overlap3.cpp:1189), Manifold::
Decompose 3D-by-connectivity (constructors.cpp:455 - the doc's "existing Decompose
primitive" is accurate, not the 2D CrossSection::Decompose). NO orient2d/orient3d
in the tree - consistent with the doc listing the static filter as B's own, not as
reuse (honest), though it is unflagged net-new code.

OPEN LIST COMPLETE vs the pre-distillation KEEP list (simplify-ab): THE BUILD /
COMPONENT-LOCAL SEED / SINGLE GLOBAL SoS / NEGATIVE WINDING / ONCE-ONLY all
present; RE-GATE + STATIC FILTER folded into pipeline/mechanism (specified, not
open). Only LOCALIZER + exterior SELECTIVE WELD dissolved (justified) and the A
FOLD dropped (owner). Nothing silently vanished - EXCEPT the localizer's cost
residue (attack 1b), which is a lost CAVEAT, not a lost item.

NO orphaned machinery references: grep confirms zero "fold"/"localizer"/"DCEL" as
live design elements in the NEW doc (the one "DCEL" is the halfedge-not-DCEL trap
warning, intentional).


================================================================================
## ATTACK 4 - STYLE / TERMINOLOGY : SURVIVE (clean)
================================================================================
- DCEL: grep docs/ - ZERO in every history doc (ExactArrangement3D.md sweep of 11
  sites confirmed complete; MaintainedEmission3D / V4ImplPlan / SweepEmit3D /
  SweepPlane3D / Boolean2 all 0). Only site is Regularize3D.md:199, the
  "halfedge-not-DCEL are the standing traps" warning = intentional, correct.
- ASCII: LC_ALL=C non-ASCII scan of Regularize3D.md = PURE ASCII. No em/en-dash,
  no smart quotes, no ellipsis glyph.
- AI-tell vocab (leverage/robust/intricate/tapestry/delve/realm/seamless...):
  none found. House voice holds.
- Magnitudes-only: no bare exact fixture counts; only file:line citations and the
  "~1.6B" magnitude (qualified). Respected.
- Repo terminology: halfedge throughout; eps = kPrecision*bBox.Scale() machine
  weld radius vs tol = user tol, stated correctly (R4's own trap avoided);
  signed-multiplicity fold described as the operand-of-origin DISCARD, not an
  operand label (correct - matches overlap3.h:40-43 CanonicalFace.mult).


================================================================================
## RECOMMENDED DOC CHANGES (NEED-CHANGE, not BREAK)
================================================================================
1. R1: replace "the re-gate catches the failure but does not repair it" with the
   honest split - validity catches non-manifold weld outcomes; IsSelfIntersecting
   is BLIND to self-folds at the merge (the merge creates the shared-vertex config
   it skips). Note wall-A C3 already killed the constant-radius bounded weld
   (x1.8..4.0 -> degenerate collapse). Add B's once-only / radial-not-projection
   defense as the rebuttal, and name the residual (two distinct points rounding
   within eps) as open.
2. Localizer / B-mechanism or open list: state the whole-component winding cost
   residue (O(#cell-components * ntri) seed casts, ~seconds/query-scale on a large
   single component like GT7081 per v5-verify-probe) + the new-surface-area
   adaptive kernel (no orient2d/orient3d in the tree today).
3. R2: state the gate is systematically CLEAN-biased; keep (i) as the sharp
   silent-miss; fix (ii)'s inverted "2*eps relaxation flags near-misses" (it
   SUPPRESSES them) and note (ii) is near-empty (B-skipped ⊆ gate-skipped);
   cross-reference that R2(i) and R1 are one blind spot.

No src/test touched (reads only). This notebook is the sole write; no commit.
