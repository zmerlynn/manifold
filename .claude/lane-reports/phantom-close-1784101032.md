# TARGETED landing: phantom-guard uncertainty must FAIL-CLOSE (Regularize3D)

Lane: complete the openscad phantom-seam guard so a filter-UNCERTAIN verdict
escalates to the total SoS predicate instead of silently concluding "no seam".
SOLE canonical editor, branch explore/sweep-plane-3d-v5, HEAD 7915174f (unpushed).
Do NOT push. ASCII only.

## The flag (land2-verify-1784100476.md, item iii)

The just-landed phantom-seam guard (RecordSeams, nPts!=2 branch) skips a pair as a
MEASURE-ZERO CONTACT when `!cleanPierce(i,j) && !cleanPierce(j,i)`. `cleanPierce`
is FILTER-ONLY: it looks for an owner edge whose two endpoints are
Orient3DFilterSign-CERTIFIED strictly off tgt's plane on opposite sides AND
EdgePiercesTri==1. A near-tangent GENUINE crossing the filter cannot certify
(endpoint filter-0, or EdgePiercesTri==-1 edge-edge z-order refusal) makes
cleanPierce return false -> the pair is SKIPPED (silent drop), not fail-closed.
That is a lossy off-corpus fallback, banned by the standing owner rule. Contrast:
the orient landing (orient-land-1784094685.md) KEPT the EdgePiercesTriSoS exact
guard for exactly this regime.

## The completion (one predicate, one implementation - no new FORM)

In cleanPierce, when the filter REFUSES (EdgePiercesTri == -1: near-tangent
endpoint sign or an edge-edge z-order the static bound cannot certify), escalate
to the existing TOTAL predicate EdgePiercesTriSoS (the same Orient3DSoS path
recordEdge / offVertexPierces already use on -1). Arms:
- certified no straddle / certified no pierce (filter r==0) -> that edge no-crosses
  (unchanged).
- certified pierce (filter r==1) -> genuine crossing found -> nPts!=2 so FAIL
  CLOSED (A.ok=false), unchanged.
- UNCERTAIN (filter r==-1) -> SoS DECIDES: pierce(1) -> genuine -> fail closed;
  no-pierce(0) -> that edge no-crosses. Never a guess, never a silent skip.

Pattern is byte-for-byte the offVertexPierces recovery lambda (lines ~1239-1244):
`int r = EdgePiercesTri(...); if (r==-1) r = EdgePiercesTriSoS(...); if (r==1)...`.

## Rails (plan)

- On corpus this MUST be a bitwise NO-OP: the filter certifies every corpus case,
  so r==-1 never occurs in cleanPierce -> escalation fires 0 times -> the new
  branch is never executed -> output bit-identical by construction.
- From-scratch build. Instrument the escalation branch (temporary fprintf) to
  COUNT fires across full Overlap3 suite + GT7081 + openscad. Expect 0. Report,
  then remove.
- FNV vs pre-edit HEAD build on resolving carriers (baseline = clean HEAD build/).
- openscad still fail-closed at F4; GT7081 still resolves.
- Mutation note (invert escalation arm: uncertain -> skip) - which pin catches it.
- docs/Regularize3D.md ledger: one sentence (soundness backstop, EG-guard class).

## Baseline

- HEAD 7915174f; `git diff HEAD -- src/ test/` empty (working tree clean on
  tracked src/test). build/ is a clean HEAD Release build (Ninja, MANIFOLD_PAR=ON,
  MANIFOLD_DEBUG=OFF); binary + libmanifold.so are pre-edit HEAD.
- Baseline FNV (harness reused from simp2-verify, linked vs build/): PokedCube
  e2a7ea215b315e32, EntangledBars 8cc4d3319b06cd46, BarsCrossZ d13e905cab2bf5fa,
  GT7863 6f8e3fcb5cd0162c, siA d844dadb573272bd, siB 65e0a03b829473a6 (match the
  simp2-verify baseline row-for-row -> build/ confirmed clean HEAD).

## FINDING (owner premise corrected): the filter does NOT certify every corpus case

First attempt (Option C): dropped the whole `su==0||sv==0||su==sv` pre-filter and
escalated EVERY filter-refused (r==-1) edge to EdgePiercesTriSoS. From-scratch
build, full suite:
- 92 escalations across 32 UNIQUE openscad pairs (NOT 0). The phantom-guard
  cleanPierce path IS reached on the openscad soup with edges whose endpoint sits
  EXACTLY on the other triangle's plane (vertex-on-face / vertex-on-edge
  T-junctions, coincident-position duplicates) - the filter refuses those signs
  (su==0), so the owner's "filter certifies every corpus case" premise is FALSE
  for this path.
- openscad test FAILED (line 1202): it no longer fails at the F4 winding-probe
  wall. EdgePiercesTriSoS PERTURBS an exactly-on-plane vertex to one side and
  reports a PHANTOM pierce -> cleanPierce true -> A.ok=false -> REOPENS the exact
  F11 seam-truncation wall the oscad-reopen lane proved closed (its rational
  reconstruction: these 308 measure-zero pairs have ZERO clean off-plane interior
  pierce). GT7081 still resolves (103s). 31/32 pass.

ROOT CAUSE: SoS answers "which way if perturbed?", NOT "is there a positive-
measure crossing?". For a measure-zero contact (endpoint exactly ON the plane)
the right answer is "touch, no crossing" - EdgePiercesTriSoS is the WRONG
escalation for the STRADDLE test; it manufactures phantoms. Option C is banned.

## Revised design (Option B): exact-complete the straddle, do NOT perturb it

The flag is about a GENUINE crossing (both endpoints EXACTLY strictly off-plane,
opposite sides) the filter can't certify. Keep the clean-straddle precondition,
but complete a filter-REFUSED endpoint sign with the EXACT predicate
(Orient3DFilterSign then Orient3DExactSign - filter-first, the blessed pattern),
NOT SoS:
- endpoint EXACTLY on-plane (exact==0): a T-junction TOUCH, measure-zero -> skip
  the edge (never perturbed into a phantom). This is the corpus regime -> output
  unchanged.
- both endpoints EXACTLY strictly off-plane, OPPOSITE sides: a genuine crossing;
  the interior pierce test (EdgePiercesTri) decides, and where its edge-edge
  z-order is filter-uncertain the TOTAL EdgePiercesTriSoS decides (a genuine
  strict-off straddle -> SoS never perturbs an on-plane endpoint here). Pierce
  found + nPts!=2 -> fail closed.
Consistent with orient-land item 2 (filter-first exact completion, never a
filter-only guess). Adds one FILTER-FIRST Orient3DExactSign caller -> docs
inventory updated.

Rail becomes: on corpus the exact endpoint completion FIRES (the T-junctions),
but every fire returns exactly-on-plane -> skip == the old filter-skip, so OUTPUT
is bitwise-identical; the genuine-straddle SoS-interior escalation fires 0 and no
pierce verdict flips. Instrument all three to prove it.

## FINDING 2: even Option B (exact endpoint straddle) flips openscad via the INTERIOR SoS

Option B kept EdgePiercesTriSoS for the interior test once an exact strict-off
straddle was established. From-scratch rebuild + full suite:
- EXACTFIRE 120 (endpoint completions, benign), OFFSTRADDLE 0 (no filter-refused
  endpoint became a straddle - the on-plane vertices skip correctly), but
  SOSINT 21: filter-CERTIFIED clean straddles (both endpoints strictly off-plane)
  whose crossing point grazes the triangle BOUNDARY (interior o1/o2/o3
  filter-uncertain) still escalate to EdgePiercesTriSoS.
- openscad FAILED again (line 1202): EdgePiercesTriSoS PERTURBS a
  boundary-grazing crossing (crossing exactly ON the triangle edge) to "inside"
  and reports a phantom pierce. So SoS is the wrong tool for the INTERIOR test
  too, not just the straddle - it manufactures a phantom from any measure-zero
  contact (endpoint on-plane OR crossing on boundary).

## Option D (LANDED): the EXACT, NON-PERTURBING, STRICT-INTERIOR completion of EdgePiercesTri

cleanPierce is the exact completion of EdgePiercesTri: every Orient3DFilterSign 0
is filter-first completed by Orient3DExactSign, and ANY exact-zero (endpoint
exactly on-plane, OR crossing exactly on the triangle boundary) is a MEASURE-ZERO
contact -> NOT a clean pierce -> skip that edge. No SoS, no perturbation. A clean
pierce requires both endpoints EXACTLY strictly off-plane on opposite sides AND
the crossing strictly inside (all three edge-edge orientations one nonzero exact
sign). This exactly matches the oscad-reopen definition ("a clean off-plane edge
piercing the STRICT interior") and its exact proof (none of the 308 measure-zero
pairs has one).

From-scratch rebuild + full suite (instrumented):
- 32/32 PASS. openscad still fail-closed at F4 (test asserts the "seam sub-face
  arrangement not exactly resolvable" detail). GT7081 resolves (99.7s).
- EXACTFIRE 206: the exact predicate DOES fire on corpus (filter-refused signs in
  the phantom-guard path - the owner's "filter certifies every corpus case"
  premise is FALSE for this path). Every fire either confirms exactly-on-plane /
  on-boundary (skip) or agrees with a filter-certified verdict.
- DELTA 0: NO cleanPierce ever returned "pierce" via a path that needed exact
  completion. So every cleanPierce verdict equals the OLD filter-only verdict ->
  the skip/fail-close decision is UNCHANGED on every corpus carrier -> OUTPUT is
  bitwise-identical by construction. (DELTA==0 is the airtight corpus rail; it is
  STRONGER than "0 escalation invocations" because it proves 0 output change over
  the entire suite, not just the 6-fixture FNV.)

Honest correction to the owner's rail: the escalation (exact completion) is NOT
0 on corpus - it fires 206 times. But it is a bitwise OUTPUT no-op: DELTA==0.

## Final rails (de-instrumented from-scratch build-phantom)

- Diff: src/overlap3.cpp only, +32/-10, the cleanPierce block; no instrumentation
  residue (grep clean); ASCII-clean. docs/Regularize3D.md caller inventory + one
  backstop sentence. Notebook.
- FNV vs pre-edit HEAD (build/), 6 carriers (PokedCube/EntangledBars/BarsCrossZ/
  GT7863/siA/siB): byte-for-byte IDENTICAL (diff -q equal).
- Full Overlap3 suite minus GT7081 (ulimit -v 6G): 31/31 PASS.
- openscad isolated (ulimit -v 6G, timeout 900): PASS - fail-closed at F4 ("seam
  sub-face arrangement not exactly resolvable"), no partial output.
- GT7081 isolated (ulimit -v 4G, timeout 900): PASS - RESOLVES (109.7s, no OOM at 4G).
- DELTA==0 across the full instrumented suite (all 32) -> bitwise-identical output
  on every corpus carrier by construction.

## MUTATION NOTE (invert the escalation arm: uncertain -> skip)

Inverting the arm means: on a filter-refused sign, SKIP the edge (treat as
no-pierce) instead of exact-completing it - i.e. revert to the OLD filter-ONLY
cleanPierce. That inverted version IS exactly the pre-edit HEAD (build/), which I
proved corpus-identical (FNV 6/6 byte-identical, full suite 32/32, DELTA==0). So
WHICH PIN CATCHES THE INVERSION ON CORPUS: NONE. The arm is a pure SOUNDNESS
BACKSTOP (EG-guard class) - it changes behavior only off-corpus, where a filter-0
edge that is EXACTLY a strict-interior transversal pierce would be silently
dropped (skip) by the inverted arm instead of exact-caught and failed-closed.

The harness IS a live oracle, proven by the two WRONG designs it bit RED:
- Option C (drop the straddle pre-filter, escalate every r==-1 to EdgePiercesTriSoS):
  openscad test RED (line 1202) - SoS perturbs an on-plane vertex into a phantom
  pierce, reopening F11. 92 escalations / 32 pairs.
- Option B (exact straddle, but EdgePiercesTriSoS for the interior): openscad RED
  again - SoS perturbs a boundary-grazing crossing into a phantom pierce. SOSINT 21.
Both flips were caught by the openscad F4 pin, so the corpus DOES guard the
PERTURBING mis-escalations hard; it just cannot guard the correct arm's inversion
because that inversion is measure-zero-identical on corpus (the flag's regime is
genuinely off-corpus). This is the EG-guard pattern exactly.

## VERDICT: LANDED (Option D). One commit. Do NOT push.

cleanPierce = the EXACT, NON-PERTURBING, STRICT-INTERIOR completion of
EdgePiercesTri (filter-first Orient3DExactSign, any exact-zero = measure-zero =
skip). Closes the land2-verify flag (near-tangent genuine crossing no longer
silently dropped) with zero corpus output change. openscad stays F4, GT7081
resolves, all resolving carriers byte-identical.
</content>
</invoke>
