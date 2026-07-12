# Maintained-order emission for 3D overlap removal (v4 design exploration)

STATUS: DESIGN EXPLORATION under crucible review, branch
explore/sweep-plane-3d-v4. This is an ARCHITECTURE study, not a pipeline. The v3
sweep-native emission (docs/SweepEmit3D.md) stays intact as the reference and the
fixture harness; nothing here changes v3's default behavior. The deliverable is a
variant adjudication + probe evidence + an honest per-variant verdict on ONE
question:

  Does restructuring emission so that cross-section state is carried CONTINUOUSLY
  across criticals - so a block-rule collapse of a near-degenerate junction
  propagates downstream BY CONSTRUCTION - DISSOLVE the wall-A coupling, or
  RELOCATE it?

The short answer this document argues, with probe evidence: maintained-order
emission DISSOLVES the coupling that is really about INPUT/PLANE independence
(v3's many independent per-plane cap computations), and this is worth doing; but
it RELOCATES the irreducible core - the near-degenerate cluster's
spread-vs-distinct RESOLUTION - into whatever structure makes that one decision.
The resolution is a real, macroscopic, coordinate-determined choice (probe 2), so
no architecture makes it free; the architectures differ only in whether they make
it ONCE (dissolution of the plane-independence wall) or drift it across events
(a new, worse hazard). The honest verdict is therefore
DISSOLVES-UNDER-PRECONDITIONS, and the precondition is the same research-grade
cluster disambiguation the WALL-A memo and the June RSI-#3 memo already name.

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
near-degeneracy from disagreeing constructed coordinates. Probe 2 below shows that
disagreeing coordinate resolutions of one degeneracy differ MACROSCOPICALLY - so
this drift is not a rounding nuisance, it is the wall reintroduced inside the
maintained structure.

RELOCATION VERDICT (as arguable): RELOCATES unless the KDS carries a correct 3D
block-rule / simultaneous-event collapse WITH cluster disambiguation - at which
point it DISSOLVES the plane-independence coupling but INHERITS the disambiguation
as its degeneracy policy's precondition, and additionally must defeat incremental
constructed-coordinate drift that the batch engine avoids for free. Net: the
weakest variant on robustness - it gives up the batch engine's from-input
re-derivation, which probe 2 shows is the very thing that makes the collapse stable.

EMPIRICAL FALSIFIER: a toy kinetic-maintenance probe on a constructed near-coincident
bundle - does a block-rule collapse keep the maintained structure order/certificate
stable, or does sub-eps certificate reordering thrash the emission? (Probe 2 runs
the batch surrogate of this and finds: order-free, but resolution is
perturbation-direction determined at macro scale - so an incremental version that
perturbs certificates via constructed drift would thrash.)

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
  arrangement's tolerance AND crosses nothing.
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

EMPIRICAL FALSIFIER / PROBE 1 (run; see below): bind Havocglass8/GT7863 cap-input
images to their track's canonical vertex and measure fan close, manifoldness, oracle
volume. RESULT: does NOT close, even binding a quarter of all extensions - the
bounded hybrid relocates.

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

RELOCATION VERDICT (as arguable): DISSOLVES cleanly IF the global collect exists and
is correct - and here is the honest catch: the concrete 3D global arrangement, with
all triple points resolved and all faces cut into a shared subdivision, IS the
RSI-#3 arrangement completion. Variant (iii) does not AVOID the research-grade work;
it NAMES it as the foundation and is honest that v3's wall exists precisely because
v3 SKIPS the global collect (doing per-slab local arrangements instead). There is no
cheaper 3D global collect than either (i)'s incremental maintenance or the full
batch arrangement - the 2D collect is cheap only because its status is a 1D order.
So (iii) does not relocate; it costs what the real fix costs. The spread-vs-distinct
disambiguation lives in the global collapse policy, made once for the whole cluster -
which is exactly where the memo says it must live.

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
collected structure - it still needs the global collapse. All roads lead to: make
the near-degeneracy resolution ONCE (globally, batch (iii), or incrementally with a
correct block-rule collapse (i)) and flow it downstream, with cluster
disambiguation. The variants differ only in the machinery that carries the one
decision, and in what robustness they trade for it.

## Probe evidence

Two env-gated probes, smallest-first, on the PROVEN engine and the real fixture.
Both are reverted from src before this document lands (default v3 behavior
unchanged); they are recorded here by mechanism.

### Probe 1 (variant ii): identity-bind does not close the fan

Mechanism: in the per-slab resolver's Extend, when a section vertex's track
terminates at a canonical 3D vertex within a ball (in x) of the cap plane, BIND the
extension to that vertex's transverse position instead of extrapolating along the
(steep) track. Provenance-selected: it uses the track's OWN terminating vertex, so
there is no nearest-neighbor ambiguity - this is the identity-bind signature the
verification lane flagged as un-probed under its true name, distinct from the
nearest-neighbor Voronoi snap. Measured on the two hand-checkable fan carriers plus
the resolving control, sweeping the ball from a fraction of eps to several eps.

RESULT: the fan does NOT close at ANY ball radius; both carriers stay
NonManifoldEmission ("unresolvable sheet contact"); the control resolves unchanged.
The bind FIRED HEAVILY - about a quarter of all cap-input extensions snapped to a
canonical vertex, and widening the ball several-fold barely changed that count and
changed the outcome not at all.

WHAT IT ESTABLISHES: binding cap-input images to their canonical 3D vertices - even
provenance-selected, even aggressively - is INSUFFICIENT. Two reasons, both already
in the memo and now confirmed from the identity angle: (a) forced-through WELD images
name no canonical to bind to, so they stay self-located; (b) the cap arrangement runs
over ALL input, so a bind that fires does not reach the weld images and/or perturbs
the arrangement elsewhere (residual-1). This is the killed candidate A (input snap)
re-confirmed via the provenance-identity selection - the SELECTION was never the
blocker; the MERGE ACTION and the weld images are. Variant (ii)'s bounded form
RELOCATES; only its coordinated-re-emission form (the research remainder) could close
the fan.

### Probe 2 (variant i): the collapse is order-free; the resolution is not free

Mechanism: a dense near-concurrence in the batch 2D engine - square chords whose
pairwise crossings pile into a sub-eps band (a subset exactly concurrent, an exact
degeneracy) - run through RemoveOverlaps2D under (a) permuted input edge order and
(b) sub-eps coordinate jitter at several magnitudes and two independent sign
patterns. Signature: retained-edge count and retained macro length.

RESULT, two clean findings:
- INPUT ORDER is BITWISE invariant: permuting the edge list (reverse, rotate) leaves
  the arrangement and the retained boundary bit-identical. The batch block-rule sweep
  is fully order-free.
- The near-degeneracy RESOLUTION is PERTURBATION-DIRECTION determined, at MACRO scale:
  two different sub-eps jitter patterns give retained boundaries differing by more
  than a chord length, each STABLE across jitter magnitude (a fraction of eps up to a
  tenth of eps). The exact-degenerate base is a third distinct, self-consistent answer.

WHAT IT ESTABLISHES: the block rule dissolves the INPUT-ORDER coupling completely,
but it does not (cannot) make the near-degeneracy RESOLUTION coordinate-free - the
resolution is a real decision with macroscopic consequences, fixed by the sub-eps
coordinate details (= the certificate order). This is the wall's mechanism in
miniature: independent resolutions of ONE degeneracy, fed slightly different
coordinates, diverge macroscopically. It also shows WHY the batch engine is robust
where an incremental KDS would not be: the batch re-derives every arrangement from
INPUT coordinates, so its resolution is single-valued and order-free; an incremental
KDS (variant i) consumes CONSTRUCTED coordinates across events, which reintroduces
exactly the cross-place resolution divergence this probe measures. The batch
from-input re-derivation is the property that makes the collapse stable, and variant
(i) trades it away.

## Convergent verdict

Per variant, the honest current answer:

- Variant (i) FULL KINETIC: RELOCATES (to certificate ordering) and additionally
  gives up the batch engine's from-input re-derivation, which probe 2 shows is what
  makes the collapse stable. Not recommended: it takes the plane-independence wall
  and swaps it for an incremental-drift wall that is at least as hard, plus it still
  needs the same cluster disambiguation. DISSOLVES only under a strictly harder
  precondition than (iii).

- Variant (ii) IDENTITY-THREADED HYBRID: DISSOLVES-UNDER-PRECONDITIONS for an
  isolated, cleanly-identified case-(1) spread twin; RELOCATES in its bounded form
  (probe 1: does not close the fan) and for dense clusters (matching ambiguity). Its
  lasting value is as the cheapest honest FRAMING - "one junction, one identity, one
  representative" - and as the additive layer a coordinated-re-emission fix would sit
  on. As a standalone bounded mechanism it is the killed input snap.

- Variant (iii) TWO-PASS GLOBAL IDENTITY: DISSOLVES the plane-independence coupling
  cleanly and does NOT relocate - but only because it IS the RSI-#3 arrangement
  completion, paid in full. It is the honest architecture: it says the wall exists
  because v3 skips the global collect, and the fix is to do it, once, globally, with
  cluster disambiguation. This is the same research-grade object the WALL-A and
  RSI-#3 memos name; maintained-order emission is a correct NAME for it, not a
  shortcut around it.

THE CROSS-VARIANT CONCLUSION: maintained-order emission dissolves the coupling that
is genuinely about per-plane INDEPENDENCE - reconciliation-by-position across
independent computations - and that dissolution is real and worth the architecture.
But probe 2 shows the near-degenerate cluster's RESOLUTION is a genuine, macroscopic,
coordinate-determined decision; no architecture makes it free. The architectures only
differ in whether they make it ONCE (variant iii, dissolving the plane-independence
wall and exposing the resolution in its cleanest single-decision form) or drift it
across events (variant i, a worse wall). The irreducible core - deciding, coherently
for a whole cluster, spread-vs-distinct BEFORE re-emission - is RELOCATED, not
removed, by every variant; it is the research remainder regardless of architecture.
So the Smith lift to 3D is genuine (the maintained structure is the right shape and
variant iii is the honest architecture), but it is NECESSARY, not SUFFICIENT: it
dissolves the easy coupling and exposes the hard one, which stays research-grade.

## RISKS (attack surface for the adversarial round)

A skeptic should try to break the following:

- R1 (the strongest attack): is variant (iii)'s "global collect DISSOLVES cleanly"
  claim just variant (i)'s relocation deferred? Attack: the global batch collect must
  still ORDER its collapse events; show a construction where the global collapse's
  cluster grouping is itself certificate-order sensitive at macro scale (probe 2 lifted
  to the global 3D collect), proving (iii) relocates into the global collapse's
  grouping predicate rather than dissolving.

- R2: probe 1's kill scope. It fired the bind on the two hand-checkable carriers plus
  the control; it did NOT separately confirm that the SPECIFIC failing junction's images
  were among those bound (the counter is aggregate). Attack: instrument the first
  failing fan's vertices and show whether they were bound-and-still-diverged (residual-1)
  or never-bound-because-weld (the weld gap) - the kill is a population argument, not a
  per-junction trace.

- R3: probe 2's fixture is open chords (an ill-posed winding field), and the
  macroscopic resolution difference could be attributed to that ill-posedness rather
  than to certificate order per se. Attack: reproduce the perturbation-direction
  macro-divergence on a WELL-POSED fixture (a union of near-degenerately overlapping
  closed polygons with a stable area oracle), to show the resolution sensitivity is
  intrinsic to the near-concurrence, not to the open-chord winding.

- R4: the spread-vs-distinct disambiguation is asserted as the irreducible core but
  never constructed. Attack: exhibit a bounded, correctness-fenced disambiguation
  policy (e.g. a Voronoi-ratio cluster test lifted from the identity-extension arc)
  that closes at least one fan carrier oracle-true without over-merging a distinct
  pair - which, if it exists, would DEMOTE the core from research-grade to bounded and
  break the whole convergent verdict.

- R5: the incremental-drift claim against variant (i) rests on probe 2's batch
  surrogate, not on a real incremental KDS. Attack: build a minimal incremental 1D-order
  maintenance across constructed events and show it EITHER thrashes (confirming the
  hazard) OR stays stable via a from-input recertification (showing (i) can recover the
  batch's robustness and the drift objection is weaker than stated).

- R6: "no fourth variant survives" is an exhaustiveness claim. Attack: propose a
  maintained structure neither the full subdivision (i) nor the global collect (iii)
  captures - for instance, a per-junction local re-mesh triggered only at flagged
  clusters (a targeted coordinated re-emission), and argue it is bounded rather than
  research-grade because it touches only flagged clusters, not the global arrangement.

Round 2 should attack R1 and R4 first: R1 is the load-bearing claim that (iii)
dissolves rather than relocates, and R4 is the single result that would overturn the
convergent verdict if a bounded disambiguation exists.
