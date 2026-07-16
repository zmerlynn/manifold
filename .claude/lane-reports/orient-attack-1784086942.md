# ATTACK: THE ESCALATING ORIENT3D - necessity, hard

Owner-commissioned. Attack the exact orient3d kernel (Orient3DExactSign + SoS
cascade) - re-prove its necessity per fire class or eliminate it. Canonical
READ-ONLY; experiments in /tmp/orient-attack/. No commits. Zero-oracle-wrong.

Branch explore/sweep-plane-3d-v5 HEAD 3776a29b.

## The kernel (src/overlap3.cpp)

- `sos::ExactOrient3D` (577): the e^0 exact orient3d sign, TOTAL (adaptive-width
  two's-complement accumulator, no window-fail). ~100 lines integer.
- `sos::SoSOrient3D` (611): Edelsbrunner-Mucke symbolic perturbation keyed by
  global vertex index; never 0 for distinct idx. 192 monomials.
- `Orient3DFilterSign` (688): Shewchuk static error-bound filter. Returns +/-1
  certified or 0 uncertain. THE FAST PATH. Kernel is only touched on filter-0.
- `Orient3DExactSign` (722): thin wrapper -> sos::ExactOrient3D.
- `Orient3DSoS` (734): filter -> ExactSign -> SoS cascade. NEVER 0.

## The 3 production callers (+ test probe + banked radial)

1. TIE CASCADE `Orient3DSoS` (734): filter-first (`if (s!=0) return s;`), then
   ExactSign, then SoS. Called from EdgePiercesTriSoS (D1/D3 predicates).
2. `EdgePiercesTriSoS` (772): edge-in-plane guard uses ExactSign directly
   (777-778: `fu==0 && fv==0 && ExactSign==0`), then Orient3DSoS x5.
3. WINDING PROBE `WindingAt` (1028): 5 filter-0 escalations to ExactSign
   (1048,1049,1055,1056,1057). reg3d-c2bx: GT7081 fires=5521, zeros=0.

## Baseline census (from reg3d-c2bx-verify notebook, to reproduce)
- GT7081: WindingAt fires=5521, zeros=0 (every graze decidably NONZERO).
- siA=0, siB=0, openscad=0, EntangledBars=0, CoplanarFold Mult2/3=0.
- BridgedCaps: fires=96, zeros=96 (ALL genuine exact-zero ties -> nullopt/reseed).

## Plan
1. CENSUS every filter-0 escalation across corpus + fixtures, by caller x cause.
2. WindingAt fires (GT7081 5521): AVOID not decide? (a) margin-max ray; (b) path reroute.
3. Tie-cascade fires: structural ties (bitwise-dup coords)? fast path?
4. EdgePiercesTriSoS fires: same.
5. Per-class verdict + re-price structural-tie / fixed-point quantization.

## FIRE CENSUS (instrumented /tmp/orient-attack, per-fixture single-test runs)

Sites: WA da/db/o1/o2/o3 = the 5 WindingAt filter-0 escalations; EG reach =
EdgePiercesTriSoS edge-in-plane guard reaches ExactSign (fu==0&&fv==0), inplane =
both ExactSign==0; SOS filt0 = Orient3DSoS filter-0 (ExactSign fires), ex_nz =
ExactSign nonzero (NEAR-TANGENT), ex0 = ExactSign==0 -> SoS (GENUINE TIE),
struct1 = ex0 repeated-point, struct2 = ex0 shared-axis-coord.

fixture           | WA da(nz/zero) | EGreach=inpl | SOSfilt0 ex_nz ex0 | struct%
------------------|----------------|--------------|--------------------|--------
siA               | 0              | 0            | 4    4   0         | -
siB               | 0              | 0            | 4    4   0         | -
GT7863            | 0              | 1284=1284    | 2435 0   2435      | 79 (48+1874)
Constructed/CFPF  | 0              | 0            | 56   0   56        | 79 (32+12)
CoplanarFold M2   | 0              | 0            | 96   0   96        | 83
CoplanarFold TJct | 24 (0/24)      | 0            | 144  0   144       | 83
CoplanarFold M3   | 0              | 0            | 144  0   144       | 83
EntangledBars     | 0              | 0            | 96   0   96        | 100
EntBarsRotated    | 0              | 0            | 96   0   96        | 67
BarsCrossZ        | 48 (0/48)      | 0            | 96   0   96        | 100
CapSeaming(FC)    | 48 (0/48)      | 0            | 160  0   160       | 90
BridgedCaps       | 192 (0/192)    | 580=580      | 664  0   664       | 100
NearCoplanarFold  | 0              | 0            | 32   8   24        | 67
openscad (FC)     | 0              | 5510=5510    | 14452 4  14448     | 78 (7976+3298)
GT7081            | PENDING (5521 nz expected)   | ...                | ...

### Three hard findings (pre-GT7081)
F1. WA escalation: da-ONLY, and on EVERY fixture except GT7081 all da fires are
    ZERO (genuine tie -> nullopt -> RobustWinding re-seeds). There the escalation
    is a pure NO-OP (0 -> nullopt == filter-0 -> nullopt). So WindingAt's exact
    caller is load-bearing on GT7081 ALONE. da is SEED-INDEPENDENT by
    construction (Orient3DFilterSign(a,b,c,p), no seed arg) => candidate 2(a)
    margin-max-ray and 2(b) path-reroute CANNOT clear a da fire. Refuted by code.
F2. EG guard: reach == inplane on EVERY carrier (openscad 5510=5510, GT7863
    1284=1284, BridgedCaps 580=580). The exact guard NEVER overturns the filter's
    "both near plane"; it always confirms in-plane. => filter-only (fu==0&&fv==0)
    is corpus-identical; the exact guard is an unexercised soundness backstop
    (would catch a filter-0 edge that is exactly OFF-plane; never occurs here).
F3. Tie cascade: ex_nz (exact SIGN load-bearing) is TINY (0 axis-aligned, 4-8
    curved/self-int). The mass (ex0) is genuine ties routing to SoS, 67-100%
    STRUCTURAL (shared-axis-coord / repeated point => orient3d==0 provable by
    comparison, no arithmetic). Non-structural remainder 0-33% (rotated bars,
    GT7863 twin, siA/siB, openscad's ~22%).

## GT7081 census (full test = 2 RemoveOverlaps3D calls, eps + eps*0.5)
WA da=11729 ALL nz (0 zero) - the 5521-per-call class, all near-tangent NONZERO.
EG reach=132=inplane=132. SOS filt0=5100 ex_nz=1104 ex0=3996 struct=23% (544+366).
GT7081 is the SOLE carrier that stresses the exact SIGN: 11729 (WA) + 1104
(tie-cascade ex_nz) near-tangent-nonzero decisions. Its 0.002deg twins.

## KEY STRUCTURAL INSIGHT (code, verified): SoS K=0 == ExactOrient3D
SoSOrient3D enumerates 24*8 monomials; the mask=0 (no perturbed factor) subset IS
the 24 real terms of ExactOrient3D, grouped as K=0 (the first group). So when e^0
is nonzero, SoS returns that same sign. => the `ex=ExactSign; if(ex!=0)return ex`
line in Orient3DSoS (caller 1) is a PERF SHORTCUT, provably redundant with SoS,
NOT a correctness caller. Confirmed bitwise by NO_SOS_SHORT prototype.

## MUTATION RESULTS (prototype evidence, /tmp/orient-attack, env-gated)
- BASE: light Overlap3 30/30 green; GT7081 resolves; openscad fail-closed.
- STRUCT_FAST (skip e^0 accum on structural ties): light 30/30 green (BITWISE:
  structClass!=0 => e^0==0 provable). SHRINKS e^0 evals to the non-structural
  remainder; SoS unchanged, no caller removed.
- NO_EG_EXACT (EG guard filter-only): light 30/30 green. Corpus-identical
  (reach==inplane everywhere). Removes caller 2's exact calls; loses an
  unexercised soundness backstop (a filter-0 edge exactly OFF-plane would be
  wrongly suppressed as in-plane - never occurs on corpus).
- NO_WA_ESC (remove WindingAt escalation) - DECISIVE:
    * light Overlap3 30/30 green (escalation is a no-op everywhere else).
    * GT7081 REGRESSES -> FAIL-CLOSED "seam sub-face arrangement not exactly
      resolvable ... filter-uncertain" (DirtyComponentUnresolved). da fires all
      revert to zero->nullopt->reseed->classify fails. So caller 3 is
      load-bearing for GT7081's resolve ALONE.
- NO_SOS_SHORT (drop the redundant ExactSign shortcut, SoS-direct): light 30/30
  green (BITWISE by K=0 subsumption).
- GT7081 combined reducibles (STRUCT_FAST+NO_EG_EXACT+NO_SOS_SHORT, WA kept):
  PASSED (still resolves oracle-true). With NO_SOS_SHORT the 5100 filt0 all route
  through SoS (ex_nz=0->ex0=5100) and resolve IDENTICALLY - empirical proof SoS
  K=0 decides the ex_nz cases bitwise. WA kept fires 11729.
- openscad combined reducibles: PASSED (stays correctly fail-closed).
- Canonical src/overlap3.cpp: 0 diff (untouched). All experiments in /tmp copy.

## VERDICT SUMMARY (per fire class)
- WINDING-PROBE (GT7081 11729): NECESSARY-BECAUSE GT7081 resolve; candidate 2
  (avoid via ray/seed) REFUTED (da seed-invariant). Only WIND re-arch avoids.
- TIE-CASCADE ex0 (genuine ties, the mass): NECESSARY (SoS). 67-100% structural,
  e^0 half SHRINKABLE (STRUCT_FAST bitwise-green) but SoS irreducible.
- TIE-CASCADE ex_nz (near-tangent nonzero): the ExactSign is SHRUNK-TO-nothing
  (SoS K=0 subsumes, NO_SOS_SHORT bitwise-green).
- EG-GUARD (in-plane detect): SHRUNK-TO filter-only (reach==inplane, corpus-safe).
- Bottom line: ExactOrient3D callers 3 -> 1 at ZERO corpus cost (EG filter-only +
  drop SoS shortcut). 1 -> 0 costs GT7081. SoS kernel irreducible.

## CANDIDATE 2 (avoid GT7081 winding fires rather than decide) - REFUTED
da = Orient3DFilterSign(a,b,c,p): NO seed argument. All 11729 fires are da.
2(a) margin-max ray and 2(b) path-reroute both change only seed/db/o1/o2/o3, never
da. RobustWinding already re-seeds on nullopt and cannot clear a da graze (da is
seed-invariant). Moving the probe p = the O4 reprobe, adjudicated insufficient
(reg3d-c2b: the 0.002deg face is near-tangent over the whole cell extent). The
ONLY genuine avoid path is the coupled-integer-flood WIND re-architecture (removes
ray-casting entirely) - named research axis, not a ray/seed tweak.

## PER-CALLER VERDICT (callers of Orient3DExactSign)
1. Orient3DSoS shortcut (738): SHRUNK-TO-nothing (REDUCIBLE). SoS K=0 group ==
   ExactOrient3D, so the `if(ex!=0)return ex` shortcut is provably redundant;
   NO_SOS_SHORT bitwise-green (light 30/30). Perf-only (saves 192-vs-24 monomials
   on ex_nz, negligible). The SoS it shortcuts to is NECESSARY.
2. EdgePiercesTriSoS guard (777-778): SHRUNK-TO filter-only. reach==inplane on
   EVERY carrier -> NO_EG_EXACT corpus-identical (light 30/30). Exact calls =
   unexercised soundness backstop (SoS can't replace: SoS never returns 0, and the
   guard needs e^0==0 detection). Filter-0 detection suffices on corpus.
3. WindingAt (1048-1057): NECESSARY-BECAUSE-GT7081. NO_WA_ESC regresses GT7081 to
   fail-closed; green everywhere else. IRREDUCIBLE for constructed probe points
   (no vertex index -> SoS impossible). GT7081 is the sole load-bearing carrier.

## KERNEL NECESSITY
- SoSOrient3D (perturbation): NECESSARY / IRREDUCIBLE. Decides coincident-coplanar
  ties (the regularization core) on general-position mesh vertices. Structural
  fast-path shrinks the e^0 half only; SoS itself has no comparison/filter/grid
  substitute for the non-structural remainder (rotated/curved/self-int ties).
- ExactOrient3D (total e^0): 2 of 3 callers reducible at ZERO corpus cost; the
  3rd (WindingAt) is its only index-free irreducible use. If WindingAt is dropped
  or WIND-rearchitected, ExactOrient3D-as-a-function folds into SoS's K=0 group.

## CALLER COUNT 3 -> 2 -> 1 -> 0 (of ExactOrient3D)
- 3 -> 2: EG guard filter-only. FREE on corpus (bitwise-green). Loses an
  unexercised soundness margin (a filter-0 edge exactly off-plane wrongly
  suppressed - never occurs on corpus). Clean minimize win.
- 2 -> 1: drop the Orient3DSoS ExactSign shortcut (folds into SoS K=0). FREE
  (bitwise-green). Leaves WindingAt as the sole ExactOrient3D caller.
- 1 -> 0: remove WindingAt exact. COSTS GT7081's oracle-true resolve (regress to
  fail-closed) OR the coupled-integer-flood WIND re-architecture (research). NOT
  free. This is where the exact e^0 genuinely earns its keep on the corpus.
- SoS kernel itself: cannot reach 0 callers. The tie cascade is the operator.

## RE-PRICING the two revisit candidates against the data
- STRUCTURAL-TIE PROOFS: 67-100% of tie-cascade genuine ties are structural
  (shared-axis-coord / repeated point). Prototype (STRUCT_FAST) bitwise-green,
  shrinks the e^0 integer accumulation to the non-structural remainder. But
  removes NO caller and NO SoS (SoS breaks every tie regardless), and perf is not
  the bottleneck (ray-cast winding is, seconds/query). LOW VALUE: shrink, not
  eliminate.
- FIXED-POINT QUANTIZATION: refuted as a decider. (i) Structural ties are
  grid-invariant exact coincidences (two boxes sharing a face) = intended
  geometry; a grid cannot remove them. (ii) The near-tangent class (GT7081
  0.002deg) gets WORSE under snapping: coarse -> exactly-coplanar ties (stage-7
  emission wall); fine -> filter still can't decide. The one thing a bounded-grid
  CONTRACT buys is replacing the adaptive-width accumulator + totality proof with
  a fixed-width predicate (real corpus e^0 spread is only 9 bits; the 112-limb
  width exists purely for a wild-magnitude totality proof) - but that is exactly
  the int256 window reg3d-s6r deleted, re-introduced under an input contract the
  library (double, eps-valid) does not have. Lateral move, not a win.
