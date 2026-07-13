# v4 IMPL crucible round-1 - ADVERSARIAL LANE (R-BC + R-CLOSURE)

Branch explore/sweep-plane-3d-v4, HEAD 2d404d57. Reasoning lane, paper-execution.
Artifacts read: docs/V4ImplPlan.md (full; stages B/C, arch delta, closure),
docs/MaintainedEmission3D.md (kill table, probes 1-4), walla-arc, v4impl-stage0
(skeleton), v4-round3 (probe 4). Source: src/overlap3.h, src/overlap3.cpp
(Canonicalize/mult), test/overlap3_test.cpp. No src/test delta (git diff empty).

VERDICT: NEED-CHANGE.

## Two structural facts established from source (not prose)

FACT 1 - operand-of-origin does NOT exist in overlap3's input. Canonicalize
(overlap3.cpp:218-305) folds the input `Manifold::Impl& in` by SIGNED MULTIPLICITY
(`it->second.mult += (par==it->second.repPar)?1:-1`; mult==0 records dropped). The
arrangement (ArrangementGeometry, overlap3.h:172) carries: verts, faces
(CanonicalFace{id,verts,normal,mult}), seams (Seam{faceId0,faceId1,vertId0,vertId1}),
face2Group, criticalXs (x-only). It carries FACE-INTERSECTION INCIDENCE and
MULTIPLICITY. It does NOT carry an operand label (A vs B). The a+b oracle is composed
SEPARATELY (test Gate5 comment "SAME operands composed and oracled"); overlap3's input
is the pre-merged self-intersecting mesh. So the escape kill's named disambiguator
"which lump is real OPERAND vs coincidence" (MaintainedEmission3D:472) has NO field to
read - operand identity was discarded upstream of overlap3 at the boolean combine.

FACT 2 - the plan's IdentityTable carries NONE of the incidence its own prose invokes.
Struct (V4ImplPlan.md:71-81): Junction{pos, members(arr.vert ids), spread(bool)};
IdentityTable{junctions, vert2junction, thinnedCriticalXs}. No incident seams, no
incident faces, no mult, no operand. The prose that justifies closure -
"which faces and seams actually meet at a junction ... which operand each face came
from" (V4ImplPlan.md:38, 472) - names data the specified structure does not hold. The
incidence IS recoverable from the arrangement (Seam.vertId0/1 == member -> faceId0/1 ->
CanonicalFace.mult), but the plan neither carries it nor specifies the lookup or the
decision rule that consumes it.

## ATTACK 1 - R-CLOSURE, the paper execution on the Havoc escape-kill case

Walking stage C's mechanism (V4ImplPlan.md:203-226) through the exact collapsed
junction the escape kill (probe 3) lost the 336-unit lump at:

1. Pass 1 builds the table on Havoc's ~212 arr.verts (skeleton: junctions=211, one
   collapse of a 2-vert sub-eps-x cluster dx=0.445 eps). The escape kill found a
   collapsed cluster "sits at a coordinate that is EXACTLY a right-operand vertex" -
   that input vertex is an arr.vert, hence a member of some junction J_lump.

2. What the table KNOWS about J_lump: pos (lex-min member, order-free), members
   (arr.vert ids), spread (classifier = maxDx<=0.7eps AND maxTr>=maxDx). That is ALL.

3. What the killed re-mesh did NOT know: it worked DOWNSTREAM on emitted cap-plane
   images (~1.7 eps), positions-only; it had to INFER clustering from those positions.

4. What the table adds over the re-mesh:
   - "which cluster is one junction vs several": YES, genuinely new - the spread
     classifier + members, order-free from input. This is real upstream information.
   - "which lump is real operand vs coincidence": NO. The table does not carry it
     (FACT 2); operand identity does not exist to carry (FACT 1). Distinguishing a
     genuine-feature vertex from a coincidental crossing-image would require reading
     the arrangement's incidence (does J_lump's member bound real faces with mult?),
     which the plan neither carries nor specifies a rule over.

5. HOW does stage C use this to decide the lump SURVIVES? THE PLAN DOES NOT SAY. Stage
   C's mechanism is "unify cap planes for junction-imaging PAIRS, re-triangulate holes,
   join strips by name." There is NO step keyed on incidence that says "J_lump bounds a
   real thin-operand feature -> preserve it / re-emit the wedge rather than collapse."
   The lump-survival decision rule is ABSENT. The classifier that IS specified is
   BIT-IDENTICAL to the escape kill's (both: maxDx<=dxThresh*eps AND maxTr>=maxDx), so
   it flags and collapses the SAME clusters, and (escape kill, invariant across every
   parameter) the collapse LOSES the 336 lump.

6. Even GRANTING a future incidence rule that marks J_lump "distinct/preserve": the
   escape kill's INSEPARABILITY (proven, MaintainedEmission3D:460-464) says every
   config that keeps a lump-region cluster distinct REOPENS the fan. So the rule yields
   FAIL-CLOSED (NonManifoldEmission), not oracle-true closure. The information lets you
   DIAGNOSE the tension and refuse correctly; it does not RESOLVE the tension. The
   tension is topological, not informational - incidence is necessary-but-insufficient.

ADJUDICATION: the deciding mechanism is HAND-WAVED, not specified. Per the round's
rule, this is the central NEED-CHANGE. The plan must EITHER (a) specify the decision
rule concretely - what incidence is computed (from members via Seam/face/mult), and
the exact rule that turns it into "close vs preserve+failclosed" - OR (b) schedule the
deciding experiment EARLY, at stage A against the existing skeleton, not discover at
stage C that the information is necessary but not sufficient. The escape-kill
inseparability makes (b) the honest move: the likely result is "Havoc stays fail-closed",
settled cheaply at A.

FIRST PROBE TO RUN (stage A, runnable against the skeleton): on Havoc's arr.verts,
build the table AND attach arrangement incidence to each junction - for each member,
its incident seams (Seam.vertId0/1==member), the faces those seams join (faceId0/1),
each face's mult. For J_lump (the junction at the right-operand vertex) vs a plain
collapsible twin (the sub-eps-x 2-vert cluster the skeleton already collapses clean),
ask: (i) is there a COMPUTABLE incidence signature separating "real feature, must
survive" from "one junction, safe to collapse"? (ii) does that signature admit a
re-emission that BOTH closes the fan AND preserves the lump, or ONLY a fail-closed
refusal? The escape kill predicts (ii) = fail-closed. Running this at stage A refutes
or confirms the stage-C premise before A/B/C are built on top of it.

## ATTACK 2 - Havoc's NEVER-BOUND problem (weld images + far tracks)

R2 proved Havoc's blocking junction images are WELD images (name NO canonical arr.vert)
plus FAR-track interpolations (defining canonical ~1e9 eps away in x) - both
EMISSION-CONSTRUCTED in pass 2. The table's vert2junction maps arr.vert -> junction;
it CANNOT map a weld image or an interpolation (neither is an arr.vert). The plan's
only answer is section 1.4 "BuildSlabs section verts BIND into the global identity
(a junction-imaging vert takes the junction's canonical)" - which NAMES the requirement
but not the BINDING PREDICATE. The two natural predicates both fail here:
- positions-based (section vert at yz near a junction's canonical) = the Voronoi/naked
  snap, killed candidate #2;
- provenance-based (section vert traceable to an arr.vert incident to the junction) =
  exactly what CANNOT reach a weld image, which by construction is a pure 2D-crossing
  section vert with no arr.vert provenance (the never-bound case R2 found).
The plan claims BuildSlabs binds BEFORE track extension (dodging the far-track
amplification) - plausible for far tracks. But the WELD-image pole (untraceable 2D
crossings that land near J but are not J's provenance) is exactly the spread-vs-distinct
ambiguity at the section level, and the plan specifies no rule that gives such a vert a
junction id. GAP: "section verts BIND into the global identity" is the sentence an
implementer improvises around; R2 already showed the natural improvisation (ball-bind)
does not reach Havoc's blocker.

## ATTACK 3 - R-BC, the stage split

(a) B independently green + meaningful? Stage B == OV3_IDBIND alone. v4-round3 CONTROL:
"IDBIND alone ball=2 (thr off): FATAL sheet contact, capArr=735" - IDENTICAL to
baseline. B changes ZERO observable outcome on any carrier; its only deliverable is an
internal invariant (Pin_CapImage_CanonicalTransverse). The plan is HONEST about this
(V4ImplPlan.md:197-199). It is a real INFRASTRUCTURE stage but a PROGRESS-VACUOUS
checkpoint: an agent lands it, greens an internal pin, moves no carrier - and it is a
permanent re-implementation of a killed-probe HALF (OV3_IDBIND).

(b) Transient regression at B? B moves EVERY junction-imaging vert corpus-wide from
self-located to canonical yz. The design's own variant-ii analysis names residual-1
(MaintainedEmission3D:205-215): moving V's image to canonical yz perturbs the cap
arrangement and "if the perturbation crosses a neighboring edge the retained region
flips." So B risks flipping a currently-RESOLVING carrier. B's fence is only "no fence
regression" + "diff confined to the transverse bind" - WEAKER than stage A's
"emission bitwise-unchanged", and the "confined" claim is ASSERTED with no skeleton
evidence (the skeleton measured order-invariance/density, never emission-consumption on
resolving carriers). CHANGE: B must PIN bitwise/oracle-invariance on every resolving
carrier, not just "tests still pass".

(c) C scope creep-proof? C's DONE (V4ImplPlan.md:226): "isolated twin oracle-true +
corpus re-adjudicated with zero oracle-wrong (SOME MAY newly resolve; others stay
recorded contracts)". "Some may" is satisfied by "none". So C passes with: isolated
synthetic twin closes (the escape kill ALREADY closes it) + every real carrier stays
fail-closed + zero oracle-wrong (trivially, nothing new resolves). Given probe 3 (Havoc
inseparable) and probe 4 (GT7863 dead zone), the two hand-analyzed PAIRS will NOT flip.
So C's realistic DONE = "isolated twin + zero new real resolves" = the SAME fail-closed
corpus as v3, reached via B+C. C's DONE names no carrier that must change state - it is
flip-OPTIONAL, the opposite of creep-proof. CHANGE: C must name a carrier it predicts
flips (with mechanism), or state plainly the honest expected deliverable is
"principled fail-closed corpus + isolated-twin closure", so C is not mistaken for
real-carrier progress.

VERDICT on the split: the DECOMPOSITION is anatomically real (transverse collapses at
B, x-residual at C = the GT7863 trace) - NOT fictitious. But NEITHER B nor C has a DONE
that requires a real carrier to change state, so the split is honest yet
PROGRESS-VACUOUS end-to-end.

## ATTACK 4 - the kill-table fence

Stage B (cap reads canonical yz) vs kills #2 (Voronoi snap) / #4 (provenance-label
merge): the CANONICAL-CHOICE difference is real (order-free from-input, C1) - not
vocabulary. But the BINDING PREDICATE (which section vert binds) is UNSPECIFIED
(attack 2); if it is yz-proximity it IS kill #2. Difference stated for the canonical,
NOT for the binding - and the binding is where probe 1 died. FLAG (partial).

Stage C (cap-plane unification for junction-imaging pairs) vs kills #5 (cap-plane
unification alone), #6 (collapse+re-mesh = escape kill), #7 (bind+unification = probe
4). Confirmed from v4-round3: probe 4 = OV3_IDBIND (== stage B) + OV3_CAPMERGE
(== stage C's cap-plane unification), run on GT7863, OUTCOME (iii) dead zone. So
STAGES B+C COMPOSED ARE THE KILLED PROBE-4 COMPOSITION, on the exact carrier whose
anatomy the plan uses to justify the B/C split. Claimed differences:
  - global order-free canonical vs probe 4's per-track terminating arr.vert: touches
    C1 (order-invariance), does NOT touch the dead zone. For GT7863's closure = only
    vocabulary.
  - name-based unification (unify precisely the twin's two cap planes) vs probe 4's
    WIDTH-threshold demotion: a GENUINE mechanistic difference. BUT it does not touch
    the load-bearing obstruction: to remove the ~1.12-eps x-residual the twin's two cap
    planes bracket, ANY unification must resolve the REAL M4 cancelled-edge macro
    content between x1 and x2 - collapse it -> SubEpsFeature backstop; keep it ->
    inconsistent macro face. The NAME says "these two images are V"; it does not say
    what to do with the macro content between them. Stage C INHERITS probe 4's dead
    zone. FLAG: difference is real but insufficient; the plan asserts "cap planes unify
    into one plane" without specifying the fate of the bracketed macro content.

Stage C vs kill #6 (escape kill): stage C's mechanism IS collapse+re-mesh; the only
thing that makes its OUTCOME differ (fail-closed, not oracle-wrong) is the
zero-oracle-wrong GATE (pin 5), NOT a new closing mechanism. That gate is legitimate
and is why this is NEED-CHANGE not BREAK - but it means stage C does not CLOSE Havoc,
it correctly REFUSES. The difference from kill #6 is a gate, not a mechanism.

## SYNTHESIS

The zero-oracle-wrong gate genuinely fences the escape kill's failure mode: Havoc and
GT7863 fail-closed correctly, not oracle-wrong. That is sound and is why the verdict is
NEED-CHANGE, not BREAK. But (1) the closure mechanism is HAND-WAVED - the IdentityTable
carries no incidence, operand identity is absent from overlap3's input entirely, and no
lump-survival decision rule exists; (2) the deciding experiment is deferred to the most
expensive stage when it is runnable at stage A on the skeleton; (3) stages B+C
reconstruct killed probe 4 with the load-bearing difference (fate of the macro content
bracketed by the twin's cap planes) unstated; (4) B's fence is weaker than A's against
the design's own named residual-1; (5) C's DONE requires no real carrier to move, so
B+C can land producing the SAME fail-closed corpus as v3.

REQUIRED CHANGES:
1. Carry incidence (incident seams/faces + mult) in the identity table OR specify the
   members->arrangement->incidence lookup; state explicitly that operand-of-origin is
   NOT available and does not become available (FACT 1), so the disambiguator is
   incidence+mult inference, not an operand label.
2. Move the R-CLOSURE deciding probe to stage A (the incidence-signature separation +
   closing-vs-failclosed test on Havoc's J_lump, runnable against the skeleton).
3. Specify the section-vert binding predicate for weld images / far tracks (attack 2),
   or record that Havoc's never-bound images make it structurally fail-closed and stop
   implying the isolated-twin win generalizes.
4. Specify how stage C's name-based cap unification resolves the macro content bracketed
   between a twin's two cap planes (the probe-4 obstruction), or accept GT7863
   fail-closed and stop using its anatomy to imply B/C closes it.
5. Strengthen B's fence to bitwise/oracle-invariance on every resolving carrier
   (residual-1), matching stage A's discipline.
6. Make C's DONE name expected carrier state-changes, or state honestly that the
   expected deliverable is principled fail-closed + isolated-twin closure.

Fence: git diff src/ test/ empty; no machinery touched. This notebook only.
