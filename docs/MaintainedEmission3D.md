# Maintained-order emission for 3D overlap removal (v4 design exploration)

STATUS: DESIGN EXPLORATION, crucible round 3 (two adversarial rounds folded), branch
explore/sweep-plane-3d-v4. This is an ARCHITECTURE study, not a pipeline. The v3
sweep-native emission (docs/SweepEmit3D.md) stays intact as the reference and the
fixture harness; nothing here changes v3's default behavior. The deliverable is a
variant adjudication + probe evidence + an honest per-variant verdict on ONE
question:

  Does restructuring emission so that cross-section state is carried CONTINUOUSLY
  across criticals - so a block-rule collapse of a near-degenerate junction
  propagates downstream BY CONSTRUCTION - DISSOLVE the wall-A coupling, or
  RELOCATE it?

The short answer, sharpened by two adversarial rounds:

- Maintained-order emission DISSOLVES the coupling that is really about per-plane
  INDEPENDENCE (v3's many independent per-plane cap computations), and variant
  (iii) - a single global two-pass collect - is the architecture that does it. It
  dissolves cleanly under two NAMED structural conditions (re-derive the collapse
  from input coordinates via an order-independent event ordering; canonicalize
  each cluster exactly once), and it does NOT relocate: round 2 confirmed the
  global batch collapse is bitwise input-order invariant across hundreds of
  permutation and jitter checks, its cluster grouping coordinate-determined.
- The irreducible core - the near-degenerate cluster's spread-vs-distinct
  RESOLUTION - is a REAL coordinate-determined decision, but on well-posed CLOSED
  input its metric consequence is EPS-BOUNDED (a local non-manifold defect), NOT
  macroscopic. Round 3 traced probe 2's "macroscopic" divergence to an open-chord
  (ill-posed winding) fixture artifact; on closed input the same near-concurrence
  resolves Lipschitz.
- The research remainder is thereby RESHAPED. It is not an open disambiguation
  policy (the consequences are eps-bounded and a fenced classifier is
  constructible) but the IDENTITY-CARRYING PIPELINE ITSELF: variant (iii) is the
  unique holder of the upstream arrangement identity that emission discards, and
  round 3's escape kill proves that identity is exactly what a bounded
  positions-only mechanism lacks - the full bounded composition closes a carrier
  to a VALID manifold that is oracle-WRONG by a macro volume. So maintained-order
  emission (variant iii) is NECESSARY and now PRECISELY SCOPED - a specified
  implementation project with named preconditions, not the vaguer "necessary, not
  sufficient" of round 1.

## The wall (settled; see SweepEmit3D.md WALL-A - do not re-derive)

One 3D junction V is emitted INDEPENDENTLY by two adjacent slabs' cap
computations. V sits sub-eps in x from a cap critical; two adjacent criticals a
sub-eps apart in x each run their OWN 2D cap arrangement, and each emits an image
of V at its own plane's track-EXTENDED position. A steep extension track (order a
hundredfold) amplifies the sub-eps x-gap into a super-eps transverse divergence;
the two images live in SEPARATE cap arrangements and meet only at the assembly
eps-weld, where the gap exceeds eps and they do not merge. The result is a twin
cap sheet plus a micro-edge, which the sheet splitter fails to pair. The density
variant (GT7081, hull, self-intersection pair) is the same near-concurrence at a
scale that refuses at the retained-piece budget before emission.

(SUPERSEDED mechanism note, v4 round-1 density measurement - see
.claude/lane-reports/v4plan-verify-density-*. Where this memo attributes the
density blowup to seam-seam crossing bundles / a combinatorial crossing pile, the
measurement corrects it: on GT7081, selfA, and selfB the dense criticals are
99.99%+ DEGENERATE CONTACTS - near-tangent face pairs, seamLen <= eps - with ZERO
M1 seam-seam triple points. The "crossing bundle" language should read
"degenerate-contact clusters". The mechanism arguments transfer verbatim - contact
x's are consumed x-only exactly as crossing x's were, so the thinning safety
argument is unchanged - but there is no triple-point completion to fear, and folding
the full contact set into junction identities stays bounded at vertex scale.)

Five bounded shortcuts are killed in that memo, all by RELOCATION: weld-radius
bump, naked/Voronoi input snap, anisotropic canonicalize, provenance-label twin
merge, sub-eps cap-plane unification. The named real fix is COORDINATED cap+strip
re-emission around a near-degenerate junction cluster: recognize the cluster as
ONE entity and rebuild its caps and all incident strips together, so both flanks
consume one shared subdivision AND one shared junction vertex.

## The 2D law, and what "maintained status" must mean in 3D

The 2D engine (src/boolean2_sweep.cpp, after Smith ch. 6-8) is the proof the
whole 3D program leans on. It runs a lexicographic sweep whose STATUS is the
y-ORDER of the live edges at the current x. Smith's 7.6.2 block rule, at each
event p, brackets the CONTIGUOUS status block straddling p and splits every
bracketed edge THROUGH p, re-inserting the remainders sorted by gradient. This
collapses a dense near-concurrence to ONE shared vertex and forces any sub-eps
residual crossing onto it. The engine is TWO passes over ONE structure: an
arrangement pass (discovers crossings, applies the block rule) and a winding pass
(emits the retained boundary) over the SAME collected arrangement. The law the 2D
campaign proved: local order-free decisions do not compose on dense
near-concurrences; the maintained order is load-bearing; everything downstream
flows through the one collapsed structure.

The 3D question is what "maintained status" means when the swept object is not a
1D edge order but a 2D CROSS-SECTION - a planar subdivision of closed
non-overlapping polygons (Emmett's #1707 "active cross section of closed,
non-overlapping polygons"). v3 answers this by NOT maintaining anything: it
rebuilds an independent 2D arrangement per slab and an independent cap arrangement
per critical, and reconciles them afterward by position (track extension + assembly
weld). The wall is precisely the reconciliation-by-position failing on a
near-degenerate junction. Maintained-order emission proposes to carry the
cross-section state ACROSS criticals so the reconciliation never happens - the one
structure IS the shared truth.

## The relocation test (the adjudication tool)

From the RSI-#3 memo (memory/rsi-3-completion-research-grade.md). A candidate
DISSOLVES the coupling only if the near-degeneracy decision is made ONCE and
everything downstream flows from it. If the decision REAPPEARS in a different
structure - certificate ordering, incremental-update predicates,
reconciliation-by-position, a matching predicate - the candidate RELOCATES and the
wall stands. Every wall-A kill is a relocation verdict; this document is exactly
that ruthless with its own variants.

## Variant space

The three variants below are ordered by how much of v3 they keep. For each: the
maintained structure, the event set, WHERE the near-degeneracy decision lives, the
relocation verdict as currently arguable, the empirical falsifier, and what of v3
survives. The corpus harness survives in all three - always.

### Variant (i): full kinetic cross-section

MAINTAINED STRUCTURE: the planar subdivision (the 2D arrangement of closed
non-overlapping polygons - the retained solid section of #1707) maintained
CONTINUOUSLY in x. The combinatorics (edge adjacency, region incidence) are
invariant between criticals; each critical is an EVENT that updates them
incrementally. This is a kinetic data structure (KDS) over a planar arrangement.

EVENT SET: exactly v3's criticals - vertex events (a face edge enters or leaves
the section, adding or removing a segment), seam-crossing events (two section
segments start or stop crossing, a section vertex appears or disappears), triple
points. Each event carries a certificate (a predicate whose failure time is the
critical x).

EMISSION FROM THE STRUCTURE: a persistent section edge sweeps a strip over its
lifetime (birth event to death event); a cap/fan is emitted at each topology-change
event as the difference in the retained region. Both flow from the maintained
structure rather than from independent per-plane computations.

WHERE THE DECISION LIVES: in CERTIFICATE ORDERING. Near-concurrent events (the
wall's cluster: criticals sub-eps apart in x) have near-equal certificate failure
times. The KDS processes events in sorted certificate order. THE HAZARD: if the KDS
processes each near-concurrent event as a SEPARATE update in sorted order, it
maintains a distinct intermediate subdivision between c1 and c2 - which is exactly
v3's two adjacent caps, now as two sequential events, emitting the junction twice.
A naive KDS RELOCATES the decision into the event queue's comparator and reproduces
the twin. The block-rule analog that would dissolve it is a DEGENERACY-COLLAPSE
policy: process all certificates failing within a sub-eps x-band as ONE
simultaneous combinatorial update, producing ONE shared section vertex. But the 2D
block rule's "contiguous block straddling p" is defined by the 1D STATUS ORDER; in
a 2D subdivision there is no 1D order - the "block" is the set of section
vertices/edges within an eps-ball of the junction in (y,z), and collapsing them is
a 2D merge that must decide, coherently for the whole cluster, whether it is ONE
junction (spread) or SEVERAL near-coincident junctions (distinct). That is the
spread-vs-distinct disambiguation the WALL-A memo names research-grade.

A SECOND hazard is specific to (i): it replaces the BATCH 2D engine (the proven
correctness anchor, which re-derives every arrangement from INPUT coordinates) with
INCREMENTAL maintenance, where a vertex CONSTRUCTED at event c1 is consumed by a
certificate predicate at event c2. Constructed-coordinate consumption across events
is the classic KDS drift: two parts of the structure can resolve the same
near-degeneracy from disagreeing constructed coordinates. Round 2 built the real
incremental structure this hazard predicts - a minimal kinetic order-maintenance
KDS over a nextafter-spaced near-coincident event stream - and measured it: a naive
constructed-coordinate KDS drifts MACROSCOPICALLY. The spread-vs-distinct grouping
count itself swings with event-processing order (from a handful of shared vertices
to two-to-three times as many, positions running away tens of eps) and is not
well-defined. So this drift is not a rounding nuisance; it is the wall reintroduced
inside the maintained structure. The SAME probe found the escape - and it closes
the variant, below.

RELOCATION VERDICT (round 3 - CLOSED, reduces to iii): the drift is real, but a
from-input RECERTIFICATION rule (re-derive every certificate and snap from ORIGINAL
input coordinates, never from updated state) restores EXACT batch-equivalence,
order-free, on every tested pencil. That is the escape - and from-input recert of
the grouping IS clustering all from-input crossings at once, which IS the batch /
global collapse, which IS variant (iii). So the incremental machinery buys nothing:
naive (i) drifts macroscopically; correct (i) reduces to (iii). Variant (i) is
therefore RETIRED as a distinct option - it is either the drift wall or (iii) wearing
KDS clothes. (Round 1 slightly overstated that (i) "trades away" the batch's
robustness with no recovery; the sharper statement is that (i) must recert from input
to be correct, and that recert collapses it into iii.)

EMPIRICAL FALSIFIER (run - round 2, R5): a minimal incremental 1D-order maintenance
KDS over a constructed near-coincident event stream. RESULT: the naive constructed-
coordinate KDS thrashes (grouping order-dependent, macroscopic), CONFIRMING the
hazard on a real incremental structure - and a from-input recertification recovers
exact batch-equivalence, showing (i) recovers robustness only by becoming (iii).

WHAT OF v3 SURVIVES: CANONICALIZE (input quantization), SEAMS (the criticals ARE
the event set), the corpus harness. DIES: the batch per-slab arrangement, the
per-slab resolver + track extension, the independent per-critical cap arrangements.
The batch 2D engine survives only DEMOTED - from batch correctness anchor to a
per-event local re-arrangement primitive.

### Variant (ii): identity-threaded hybrid

MAINTAINED STRUCTURE: v3's batch per-slab arrangements UNCHANGED (the proven
engine), PLUS a persistent naming of section vertices/edges carried ACROSS slabs -
cross-critical identity. Tracks, seams, and the provenance chains already carry
PARTIAL identity (per-input-edge subdivision); the wall-A trace showed emission is
positions-only end-to-end. Variant (ii) adds a JUNCTION-identity layer: one 3D
junction = one identity = one canonical representative (its arr.vert V and V.yz).

EVENT SET: v3's criticals, unchanged. But each critical's cap, instead of
self-locating each junction by extending tracks to its own plane, BINDS the
junction's section-vertex images to the canonical representative. Both caps place V
at V.yz, so their images differ only sub-eps in x and weld.

WHERE THE DECISION LIVES: in IDENTITY ESTABLISHMENT - deciding which section vertex
in each slab is the image of junction V (a matching), then which canonical position
represents V. The bind itself is downstream of that matching.

RELOCATION VERDICT (as arguable): TWO relocation pressures, and probe 1 resolves the
bounded form empirically.
- (a) Choosing V.yz per plane. For the cap's TRANSVERSE position this genuinely
  dissolves the self-location: V.yz is the SAME for both caps, computed once, so the
  two images are sub-eps-x apart and weld - no track amplification. This is the real
  win and the correct shape for the case-(1) spread twin (sub-eps x). BUT the cap
  ARRANGEMENT at c1 runs over V.yz AND all other cap input at c1; moving V's image
  from track-extended to V.yz perturbs that arrangement, and if the perturbation
  crosses a neighboring edge the retained region flips (WALL-A residual-1,
  cluster-collapse re-emission). So the decision RELOCATES from "twin at the assembly
  weld" to "arrangement topology at the cap plane" unless the move stays within the
  arrangement's tolerance AND crosses nothing. Round 3 per-junction correction: on
  GT7863's actual blocking junction the transverse bind DOES collapse the six
  incident images to ONE canonical (y,z) bitwise - the real win - yet a super-eps
  residual (~one-and-a-tenth eps) survives PURELY in x, because the two images sit
  at two cap PLANES that a section-plane (y,z) bind cannot merge. So the
  "sub-eps-x, welds" optimism is FALSIFIED for the real blocker; it holds only for
  an isolated sub-eps-x twin, which no real carrier's first blocker exhibits.
- (b) Identity establishment across a critical where the arrangement changes
  discontinuously (the M4 cancelled-edge disagreement, measured to strict-FP). Left
  and right sections can disagree macroscopically about a cancelled edge's
  subdivision; matching section vertices across that discontinuity is itself a
  near-degeneracy decision, which RELOCATES the coupling into the matching predicate.
  For a dense cluster carrying multiple true vertices within construction noise,
  "which vertex is V" is ill-posed - the spread-vs-distinct disambiguation again.

So variant (ii) DISSOLVES the coupling for an ISOLATED, cleanly-identified case-(1)
spread twin (one junction, sub-eps move, no neighbor crossing), and RELOCATES for
dense clusters (matching ambiguity) and wherever the bind perturbs the cap topology
(residual-1). "Is this just the collapse decision wearing architecture clothes?" -
YES: reduced to a bounded bind rule it IS the killed candidate A (input snap) with a
provenance selection criterion; the provenance selection fixes the SELECTION
robustness (survives >eps gaps) but not the merge ACTION (residual-1) or the weld
images (which name no canonical). It becomes genuinely new only when it ALSO does
coordinated re-emission - i.e. when it becomes the research remainder.

EMPIRICAL FALSIFIER / PROBE 1 (run; see below, per-junction corrected in round 2):
bind Havocglass8/GT7863 cap-input images to their track's canonical vertex and
measure fan close, manifoldness, oracle volume. RESULT: does NOT close at any bind
radius. Round 2 replaced the aggregate counter with a per-junction trace: Havoc's
blocker is NEVER-BOUND (weld images name no canonical, far tracks have no near
canonical); GT7863's is BOUND-AND-STILL-DIVERGED (yz collapses bitwise, residual
super-eps cap-plane x-gap). The bounded hybrid relocates on both, for these
corrected reasons.

WHAT OF v3 SURVIVES: EVERYTHING (batch engine, canonicalize, seams, caps, strips) -
variant (ii) is ADDITIVE (an identity-establishment pass + a bind rule). This is its
appeal: minimal, testable, non-destructive. It is the cheapest honest FRAMING of the
collapse decision - which is a virtue for exposition even though probe 1 shows the
bounded version does not close the wall.

### Variant (iii): two-pass global identity (the 2D two-pass, lifted)

MAINTAINED STRUCTURE: a GLOBAL 3D arrangement - the full set of junction vertices,
collapsed once, globally, before any emission - analogous to the 2D engine's
CollectArrangement (first pass) feeding a winding/emission SECOND pass over the SAME
collected structure. Not maintained incrementally (unlike i); collected in one batch
pass, then emission reads it.

EVENT SET: none in the KDS sense. The collect pass processes all criticals but
produces a STATIC global vertex set + incidence; the emission pass is a read.

WHERE THE DECISION LIVES: in the GLOBAL COLLECT pass - a 3D block-rule analog
collapses all near-concurrent junctions to shared vertices ONCE, globally, before
any cap or strip. This is the STRONGEST dissolution: no per-plane self-location
exists because positions come from the global vertex set, not from per-slab track
extension; every cap and every strip reads the same collapsed structure.

RELOCATION VERDICT (round 2 - DISSOLVES, does NOT relocate, under two named
conditions): R1 attacked exactly this claim, that (iii)'s global collect is (i)'s
relocation deferred with its cluster grouping secretly order-sensitive. It is not.
Across hundreds of sweep-order checks (exact-degenerate, near-concurrent, and
jittered) the global batch collapse is BITWISE input-order invariant, and the
one-crowd/two-crowd grouping is coordinate-determined and order-free. What MAKES it
stable are two structural conditions, now the DESIGN REQUIREMENTS of the 3D lift:
  (1) the collapse is re-derived from INPUT coordinates each event via an
      order-INDEPENDENT event ordering (exact-lex event set + exact predicates),
      never consuming a constructed coordinate in a discovery-order-dependent way;
  (2) each near-degenerate cluster is CANONICALIZED exactly once (one global batch
      fixes one representative choice; the eps-level ambiguity of that choice then
      has eps-bounded consequence on the closed field).
v3's wall is precisely the violation of (2) - per-slab INDEPENDENT canonicalizations
of one junction, whose two eps-different representatives diverge past eps in the
steep-track neighborhood; variant (i)'s naive KDS violates (1) by consuming
constructed coordinates across events. A single global collect satisfies BOTH
conditions STRUCTURALLY: its input is closed solids ((2)'s well-posedness is free)
and one batch pass IS the from-input, canonicalize-once machine. The honest catch
remains: the concrete 3D global arrangement - all triple points resolved, all faces
cut into a shared subdivision - IS the RSI-#3 arrangement completion. Variant (iii)
does not AVOID the research-grade work; it NAMES it as the foundation, and its cost
is the real fix's cost. The spread-vs-distinct disambiguation lives in the global
collapse policy, made once for the whole cluster - exactly where the memo says it
must live.

EMPIRICAL FALSIFIER: a global-collect-then-emit on Havocglass8 - take the union of
all cap-plane arrangements' junction images, run a global 3D block-rule collapse over
them (with cluster disambiguation), re-emit caps+strips from the collapsed set, and
check the twin becomes ONE vertex WITHOUT opening a macro hole (residual-1) and that
distinct-near-degenerate vertices are NOT over-merged (residual-2). This is near the
full research build; a bounded surrogate is probe 1 extended to an all-criticals
merge.

WHAT OF v3 SURVIVES: CANONICALIZE, SEAMS, the corpus harness. The batch 2D engine
survives as the per-critical collect KERNEL within the global pass (each section is
still a 2D arrangement, but its vertices bind into the global set). DIES: per-slab
independence, track extension, assembly-weld reconciliation.

### On a fourth variant

No fourth variant survives scrutiny as distinct. Candidates considered and folded:
a "status = 1D order of a fixed reference direction" reduces to v3's per-slab
sections (the section IS that status, and it is 2D, not 1D). A "maintain only the
junction skeleton, not the full subdivision" is variant (iii) with a thinner
collected structure - it still needs the global collapse. Round 3 BUILT R6's
concrete proposal - a per-junction local re-mesh triggered only at flagged clusters,
a targeted coordinated re-emission that touches no global arrangement - and ran it
end-to-end (the escape kill below). It does not survive as a bounded escape: it
closes an isolated twin but cannot close a real carrier oracle-true, for the same
missing-identity reason (iii) exists to supply. All roads lead to: make the
near-degeneracy resolution ONCE - globally, batch (iii) - and flow it downstream,
with cluster disambiguation. Variant (i) reduces to (iii) (retired), so (iii) is the
one architecture; the variants differ only in the machinery that carries the one
decision, and in what robustness they trade for it.

## Probe evidence

Env-gated probes, smallest-first, on the PROVEN engine and the real fixtures. All
are reverted from src before this document lands (default v3 behavior unchanged);
they are recorded here by mechanism. Probes 1 and 2 are the round-1 falsifiers,
corrected by the round-2 verification lanes; probe 3 (the escape kill) and probe 4
(the last composition) are the round-3 centerpieces. Exact counts live in the lane
notebooks (.claude/lane-reports/v4-*); this section carries magnitudes and kind.

### Probe 1 (variant ii): identity-bind does not close the fan (per-junction)

Mechanism: in the per-slab resolver's Extend, when a section vertex's track
terminates at a canonical 3D vertex within a ball (in x) of the cap plane, BIND the
extension to that vertex's transverse (y,z) instead of extrapolating along the
(steep) track. Provenance-selected: it uses the track's OWN terminating vertex, so
there is no nearest-neighbor ambiguity - the identity-bind signature, distinct from
the nearest-neighbor Voronoi snap. Measured on the two hand-checkable fan carriers
plus the resolving control, sweeping the ball from a fraction of eps to several eps.

RESULT: the fan does NOT close at ANY ball radius; both carriers stay
NonManifoldEmission ("unresolvable sheet contact"); the control resolves unchanged.

WHAT IT ESTABLISHES: round 1 argued this from an aggregate counter (about a quarter
of all cap-input extensions bound, outcome unchanged); round 2 replaced that with a
per-junction trace of the SPECIFIC blocking fan and found the two poles of failure,
on different carriers.
- Havocglass8's blocker is NEVER-BOUND. Its incident verts are WELD images (which
  name no canonical to bind to - structurally unbindable) plus FAR-track
  interpolations whose defining canonical vertex is order-1e9 eps away in x (a
  nearly-x-parallel track outside any finite ball). No bind reaches this junction -
  not through a fixable selection miss, but structurally.
- GT7863's blocker is BOUND-AND-STILL-DIVERGED. At a wide-enough ball all six
  incident verts bind bitwise to ONE canonical (y,z) - the transverse divergence
  collapses to bit-identity, the real win variant (ii) is credited with - yet a
  residual ~one-and-a-tenth eps micro-edge SURVIVES purely along x, because the two
  images sit at two cap PLANES that a section-plane (y,z) bind cannot merge.
So the bounded bind is INSUFFICIENT on both, for corrected reasons: the SELECTION was
never the blocker (Havoc has no canonical to select; GT7863 selects the right one).
The blockers are the weld images / far tracks (Havoc) and the cap-plane x-separation
(GT7863) - neither of which a positions-only bind touches. This is the killed
candidate A (input snap) re-confirmed via the provenance-identity selection. Variant
(ii)'s bounded form RELOCATES; only its coordinated-re-emission form (the research
remainder) could close the fan.

### Probe 2 (variant i): the collapse is order-free; the resolution is eps-bounded on well-posed input

Mechanism (round 1): a dense near-concurrence in the batch 2D engine - square chords
whose pairwise crossings pile into a sub-eps band (a subset exactly concurrent, an
exact degeneracy) - run through RemoveOverlaps2D under permuted input edge order and
sub-eps coordinate jitter (several magnitudes, two independent sign patterns).
Signature: retained-edge count and retained macro length.

ROUND-1 RESULT (one half CORRECTED in round 3):
- INPUT ORDER is BITWISE invariant. Permuting the edge list (reverse, rotate) leaves
  the arrangement and the retained boundary bit-identical (round 2 reconfirmed with an
  EXACT == compare, not a fixed-precision print). The batch block-rule sweep is fully
  order-free.
- Round 1 also reported the RESOLUTION diverging at MACRO scale (retained boundaries
  differing by more than a chord length under two sub-eps jitter patterns). ROUND 3
  REFUTED that framing: the divergence was an OPEN-CHORD fixture artifact. An open
  curve has undefined winding at its dangling ends, so a sub-eps change in the central
  near-concurrence can flip an ENTIRE chord's retention - an ill-posed winding field,
  not a property intrinsic to any near-concurrence. Refuted independently twice (an R3
  well-posedness lane and a fidelity lane that could not reproduce the divergence on a
  generic near-pencil).

CORRECTED RESULT (round 3, well-posed input). Rebuild the SAME near-concurrence on a
CLOSED (well-posed winding) field and the macro divergence VANISHES. The fixture,
recorded exactly so it is reconstructible (the round-2 lesson - probe fixtures must be):
  - N chords through a shared sub-eps band near the origin, nExact of them exactly
    concurrent (an exact degeneracy); at unit scale N=9, nExact=3, band ~0.02 eps.
  - Close EACH chord into a thin triangle of MACROSCOPIC area and run RemoveOverlaps2D
    as a UNION (Add). Oracle = union area (a Lipschitz function of vertex positions)
    plus retained length.
  - Perturb by sub-eps per-vertex jitter (deterministic sign patterns, magnitudes 1e-6
    eps up to 0.9 eps) and by input-order permutation.
RESULT: union area is essentially fixed and tracks jitter LINEARLY, with area spread on
the order of a TENTH of an eps (about 0.15 eps at unit scale) and retained-length
spread under an eps - EPS-BOUNDED, not macroscopic - up to 0.9-eps jitter. This holds
for a MONOTONE union AND a SIGN-CANCELLING (alternating-orientation, closed) field
alike, so the discriminator is CLOSEDNESS, not monotonicity. The local topology DOES
still change (the used-vertex count shifts; spread-vs-distinct genuinely varies), but
with eps-bounded metric consequence.

WHAT IT ESTABLISHES: the block rule dissolves the INPUT-ORDER coupling completely
(bitwise). The near-degeneracy RESOLUTION is a REAL coordinate-determined decision -
the local topology genuinely changes - but on well-posed CLOSED input its metric cost
is a BOUNDED LOCAL DEFECT (eps-scale in area and length), NOT the macroscopic
divergence round 1 reported. Structural backing: a union's boundary is a Lipschitz
function of vertex positions, so a sub-eps near-concurrence cannot flip a macro
region's inside/outside; it only rearranges the boundary locally by O(eps). The 3D
wall lives on CLOSED input and manifests as exactly that eps-bounded local defect - a
local non-manifold twin. HONEST GAP: round 1's original "square chord" coordinates were
reverted and are not recoverable from the prose, so round 1's specific macro MAGNITUDES
are UNVERIFIED (not falsified - the open-field ill-posedness is independently confirmed
as their likely source). The well-posed refutation above is a fresh, reconstructible
fixture and stands on its own.

### Probe 3 (variant iii / the fourth-variant escape): the bounded composition closes the carrier ORACLE-WRONG

This is round 3's centerpiece and the sixth killed candidate. It builds the ONE
composition prior kills left unprobed - not selection alone (Voronoi snap), not
collapse alone (weld-bump), not binding alone (probe 1), but ALL THREE PLUS local
re-emission: disambiguate -> collapse -> REBUILD incident geometry (targeted
per-junction re-mesh) -> verify. This is exactly R6's proposed bounded fourth
variant, made concrete and run end-to-end on Havocglass8 (oracle = a+b Boolean
union; volume + winding-number reference).

THE COMPOSITION, in three fenced pieces:
- CLASSIFIER (the spread-vs-distinct disambiguation R4 said was asserted but never
  constructed): predicate = sub-eps in x (maxDx <= dxThresh*eps) AND
  transverse-dominated (maxTr >= maxDx). It FIRES on every case-(1) spread twin and
  EXCLUDES the case-(2) DISTINCT pair (~72 eps apart in x): across the whole run the
  maximum collapsed-cluster dx is ~0.64 eps, so the 72-eps pair is never even
  clustered. The MUST-NOT-MERGE fence HOLDS.
- COLLAPSE: snap each flagged cluster to one representative.
- BOUNDED LOCAL RE-MESH: per-plane, remove the coincident cap DOUBLE-SHEET and
  re-triangulate the planar cap hole, KEEPING the load-bearing strips (remapped).
  Per-plane bucketing keeps every hole planar.

RESULT: the composition drives the WHOLE Havocglass8 to a VALID manifold - genus-0
(= oracle genus), winding-correct at essentially every grid point, no open edges, no
fan. But its VOLUME is oracle-WRONG by |dv| ~ 336 volume units, about 14000x the eps
oracle bound. And that error is INVARIANT: identical across ball radius, dx
threshold, representative choice (centroid vs member), and collapse-vs-retriangulate
mode (a plain collapse-and-drop with NO re-triangulation gives the SAME |dv| ~ 336) -
the coordinate-determined-resolution signature.

MECHANISM - a macro-region FLIP, localized and physical. The missing volume is a
~336-unit solid LUMP that vanishes; it is not a re-mesh defect (the cap
re-triangulation and strip drops are volume-neutral to well under a unit). The lump
lies INSIDE the thin right operand's bounding box (a thin wedge), and a collapsed
cluster sits at a coordinate that is EXACTLY a right-operand vertex. So the lost lump
is part of the thin operand feature. Its near-degeneracy with the glass surface is
what PRODUCES the wall-A coincident-sheet fans in the first place; collapsing to close
those fans DESTROYS part of the wedge. Closing the fan and preserving the thin feature
are in genuine tension that emitted positions cannot arbitrate.

INSEPARABILITY (proven, not asserted): sparing individual clusters near the lump one
at a time - every configuration that CLOSES the fan LOSES the 336 lump, and every
configuration that keeps a lump-region cluster distinct REOPENS the fan. There is no
single-cluster sparing that both closes and restores the lump. The macro flip is
COUPLED to the closure.

WHAT IT ESTABLISHES: to close the carrier you MUST collapse the wide near-degenerate
clusters; collapsing them yields a valid manifold that is quantifiably oracle-wrong.
Havoc's near-degenerate region admits exactly two positions-only outcomes - collapse
(valid manifold, ~336 short) or no-collapse (non-manifold wall). The oracle-TRUE
topology is UNREACHABLE from emitted positions; reaching it needs the upstream
arr.vert identities that emission discards (which cluster is one junction vs several,
which lump is real operand vs coincidence) - i.e. exactly the information variant
(iii)'s global identity carries and emission does not. R4's core is therefore NOT
demoted: no bounded positions-only disambiguation closes the carrier oracle-true.

THE POSITIVE RESULTS (do not bury them - they refine round 1, they do not overturn
this verdict):
- The CLASSIFIER is CONSTRUCTIBLE and FENCED. R4 challenged "the disambiguation is
  asserted, never built"; it is now built (sub-eps-x AND transverse-dominated) and it
  cleanly separates every case-(1) twin from the case-(2) distinct pair.
- An ISOLATED case-(1) twin IS boundedly CLOSABLE by the local re-mesh: collapse +
  per-plane planar cap re-triangulation closes it to a clean 2-manifold, volume-
  neutral, no residual hole. This REFUTES round 1's implication that no bounded local
  mechanism closes any twin - the bounded per-junction re-mesh is a real partial
  mechanism.
- It does not FLIP any corpus carrier because real carriers have ENTANGLED clusters:
  the closing collapse is coupled to a macro-region flip (the inseparability above).
  The isolated-twin win does not compose into a carrier win.

HONEST CAVEATS: this is a FAMILY-of-schemes result - the collapse+re-emit family across
all its parameters fails; it is not a proof of impossibility for every conceivable
positions-only scheme (the invariance, inseparability, and thin-feature construction
make that strongly implausible without the global collect, but it is not proven).
And the baseline never resolves, so the pre-collapse surface cannot be inspected to
FULLY separate a collapse-induced residual-1 flip from a sweep inaccuracy the collapse
merely exposes. Either way the bounded positions-only escape does not reach oracle-true
and the correction needs the arrangement info emission discards.

### Probe 4 (the last composition): section-plane bind + cap-plane unification on GT7863

The remaining untested composition, run in round 3. Probe 1 and the wall-A round-2
cap-plane probe were each killed SEPARATELY, each for an objection the OTHER half
appears to dissolve: cap-plane unification died on section-plane self-location (which
the bind fixes bitwise); the bind died on the cap-plane x-residual (which unification
removes). Compose them on GT7863 (a PAIR, so the a+b oracle exists) and ask whether
its junction closes.

MECHANISM: OV3_IDBIND binds each track's terminating canonical vertex's (y,z)
(fixing the transverse self-location); OV3_CAPMERGE demotes a sub-eps-width built
slab so EmitCaps merges the two bounding criticals into ONE cap plane (removing the
cap-plane x-gap). Both env-gated; default-off is a bitwise no-op; reverted after.

RESULT: OUTCOME (iii) - the junction does NOT close, at any (bind ball, merge
threshold). The two objections are BOTH present, on opposite sides of a threshold
with no closing window between them:
- Below the twin's ~1.1-eps cap-plane gap (small merge threshold): the bind collapses
  the six incident images to one canonical (y,z), but the cap-plane x-separation
  survives untouched -> SplitTouchingSheets sheet contact (NonManifoldEmission).
  Ball-invariant: widening the bind does not help, because the x-residual blocks, not
  the yz selection.
- At any threshold large enough to start merging sub-eps runs (well below the twin's
  own gap): the demotion collapses REAL macro content (the M4 cancelled-edge
  disagreement) across a run wider than eps and trips the ComputeCap wide-run
  SubEpsFeature FIDELITY BACKSTOP first (SubEpsFeature) - before any threshold large
  enough to reach the twin's gap.
So the twin's cap-plane gap sits in a DEAD ZONE: unreachable below the backstop
threshold (x-residual survives) and shadowed by the backstop above it. On a real dense
carrier no positions-only width threshold isolates the twin's cap-plane gap from the
arrangement's genuine macro content - the SAME entangled-clusters inseparability probe
3 found on Havoc, here at GT7863's cap-placement stage. The composition INHERITS BOTH
kills rather than escaping either. No tier-1 carve-out for the bound-and-diverged
sub-class; the kill table is complete.

## The killed candidates (summary)

Six bounded, positions-only mechanisms, all killed by RELOCATION or by
oracle-wrongness. The first five are the wall-A memo's (SweepEmit3D.md); the sixth is
round 3's:
1. weld-radius bump - relocates (collapses without re-emission; fan reopens).
2. naked / Voronoi input snap - relocates (selection was never the blocker).
3. anisotropic canonicalize - relocates (scoped to the distinct pair, misses the twin).
4. provenance-label twin merge - relocates (positions-only, names no canonical for
   weld images).
5. sub-eps cap-plane unification - relocates / over-merges (a shared cap x is not a
   shared point; wide-run backstop fires on real content). Its composition WITH the
   identity-bind (probe 4) also fails: dead zone, both objections present.
6. spread-vs-distinct classifier + collapse + bounded local re-mesh (probe 3) - the
   full unprobed composition. CLOSES the carrier to a valid manifold that is
   oracle-WRONG by a macro volume (~14000x the eps bound), via a macro-region flip
   inseparable from closure. The oracle-true topology is unreachable from emitted
   positions.

The pattern across all six: a bounded positions-only mechanism either relocates the
near-degeneracy decision into another structure (1-5) or resolves it to a valid but
oracle-wrong answer (6). Both failures are the same missing information - the upstream
arrangement identity that only variant (iii) carries through emission.

## Convergent verdict

Per variant, the honest current answer:

- Variant (i) FULL KINETIC: RETIRED as a distinct option. Round 2 built the real
  incremental KDS and measured it: naive constructed-coordinate maintenance drifts
  macroscopically (the grouping count is not well-defined), and the only correct form
  - from-input recertification - IS the batch global collapse, i.e. variant (iii). So
  (i) is either the drift wall or (iii) wearing KDS clothes; the incremental machinery
  buys nothing.

- Variant (ii) IDENTITY-THREADED HYBRID: its bounded form is KILLED (probe 1,
  per-junction: Havoc never-bound, GT7863 bound-and-still-diverged; probe 4: the
  bind+unification composition sits in a dead zone). Its lasting value is as the
  cheapest honest FRAMING - "one junction, one identity, one representative" - and as
  the additive layer a coordinated-re-emission fix would sit on. As a standalone
  bounded mechanism it is the killed input snap.

- Variant (iii) TWO-PASS GLOBAL IDENTITY: THE architecture. It DISSOLVES the
  plane-independence coupling and does NOT relocate - confirmed by R1 under two named
  structural conditions (from-input re-derivation; canonicalize each cluster once),
  which are the design requirements of the 3D lift. And it is the UNIQUE holder of the
  information probe 3 proved necessary: the upstream arrangement identity that emission
  discards is exactly what a bounded positions-only mechanism lacks, so only the global
  identity collect can reach the oracle-true topology. It IS the RSI-#3 arrangement
  completion, paid in full - not a shortcut around the research-grade object, the
  correct name for it.

THE CROSS-VARIANT CONCLUSION, sharpened. Maintained-order emission dissolves the
coupling that is genuinely about per-plane INDEPENDENCE - reconciliation-by-position
across independent computations - and that dissolution is real, worth the architecture,
and now has NAMED preconditions (R1's two conditions). The near-degenerate cluster's
RESOLUTION is a real coordinate-determined decision, but round 3 corrected its cost:
on well-posed closed input the consequence is EPS-BOUNDED (a local non-manifold defect,
probe 2), not macroscopic. So the research remainder is RESHAPED. It is NOT an open
disambiguation policy - the consequences are eps-bounded and a fenced spread-vs-distinct
classifier is constructible (probe 3). It IS the IDENTITY-CARRYING PIPELINE itself:
probe 3 shows a bounded positions-only composition closes the carrier to a VALID
manifold that is oracle-WRONG by a macro volume, because the oracle-true topology needs
the upstream identity emission throws away. Variant (iii) is the one architecture that
carries that identity. So the Smith lift to 3D is genuine AND variant (iii) is
NECESSARY and now PRECISELY SCOPED: a specified implementation project (global identity
collect + a fenced cluster classifier, with R1's two preconditions), not the vaguer
"necessary, not sufficient" research-grade unknown of round 1. The hard object is named,
bounded, and buildable - it is the pipeline, not an open theory.

## What the adversarial rounds found (R1-R6, run)

The six risks below were the round-1 attack surface. Rounds 2 and 3 ran every one; the
outcomes are folded into the sections above and recorded here as the audit trail. None
overturned the convergent verdict; two sharpened it materially (R3, R4/R6).

- R1 (the load-bearing attack: is (iii)'s global collect just (i)'s relocation
  deferred?) - SURVIVES. The global batch collapse is BITWISE input-order invariant
  across hundreds of sweep-order checks, its cluster grouping coordinate-determined and
  order-free. (iii) DISSOLVES, does not relocate, under two now-named structural
  conditions (from-input re-derivation; canonicalize each cluster once) - the design
  requirements of the 3D lift.

- R2 (probe 1's kill was an aggregate counter, not a per-junction trace) - UPHELD, kill
  re-grounded. Per-junction: Havoc's blocker is NEVER-BOUND (weld images + far tracks,
  structural); GT7863's is BOUND-AND-STILL-DIVERGED (yz collapses bitwise, residual
  super-eps cap-plane x-gap). The kill stands on both, on better evidence.

- R3 (probe 2's macro divergence may be an open-chord ill-posedness artifact) - CONFIRMED
  as an artifact. On well-posed CLOSED input the same near-concurrence resolves Lipschitz
  (eps-bounded area/length spread up to 0.9-eps jitter). Probe 2's "macroscopic" WORDING
  is corrected; the resolution is real but its cost is a bounded local defect.

- R4 (the spread-vs-distinct disambiguation is asserted, never constructed - a bounded
  policy would overturn the verdict) - CONSTRUCTED and FENCED, verdict NOT overturned.
  The classifier exists (sub-eps-x AND transverse-dominated) and separates every case-(1)
  twin from the case-(2) distinct pair; but the full bounded composition (probe 3) closes
  Havoc to a valid manifold that is oracle-WRONG by a macro volume. The core is not
  demoted - it is missing upstream identity, not a missing bounded predicate.

- R5 (the incremental-drift claim rests on a batch surrogate, not a real KDS) - the real
  KDS was BUILT. Naive constructed-coordinate maintenance drifts macroscopically
  (grouping count order-dependent); from-input recertification recovers exact
  batch-equivalence - and that recert IS the global collapse (variant i reduces to iii).

- R6 ("no fourth variant survives" - propose a per-junction local re-mesh at flagged
  clusters, argue it is bounded) - the fourth variant was BUILT and run end-to-end (probe
  3). It closes an ISOLATED twin (a real partial mechanism, refuting round 1's blanket
  claim) but cannot close a real carrier oracle-true: carrier closure requires collapses
  that flip a macro region, inseparable from the closure, unreachable from positions.

NEXT (implementation, not adjudication): variant (iii) is now a specified project -
a global identity collect (RSI-#3 arrangement completion) satisfying R1's two
preconditions, plus the fenced cluster classifier of probe 3 lifted to 3D. The
adjudication is closed; the remainder is engineering the identity-carrying pipeline.

## Stage-A0 addendum (measured after this doc converged)

The implementation crucible's stage A0 (V4ImplPlan.md, RECORDED RESULT +
verification round) answered the sufficiency question this doc left open, and
the answer sharpens the convergent verdict:

- The identity table's incidence channels are DEGENERATE on the hard carriers:
  signed multiplicity is constant positive in the union-composed regime, the
  Canonicalize fold is a pure bijection there (nothing cancels - measured
  population-wide, calibrated against positive controls), and every other
  incidence signature is identical between must-survive and must-collapse
  clusters.  NOT-SEPARABLE at eps quantization, pre-fold and post-fold.
- No ambiguity witness exists: the oracle union is a pure function of the
  oriented face soup.  The hard carriers are NOT ill-posed; the deciding
  winding structure lives at sub-eps scale.  THE BARRIER IS PRECISION, NOT
  INFORMATION.
- Consequence for this doc's verdict: variant (iii) at eps precision delivers
  the independence-coupling dissolution, the density collapse, and the
  consistency pins - but NOT the hard-carrier resolutions; those require the
  exact-arithmetic arrangement completion (the RSI-#3 object), which is the
  named upgrade path, not a defect of the architecture.  Fail-closed remains
  the honest terminal state for those carriers at this precision model.
