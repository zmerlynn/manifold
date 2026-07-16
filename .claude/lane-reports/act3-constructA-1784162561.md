# THIRD-ACT lane A: two constructions of Smith-style topological robustness for SELF-overlap

Branch explore/sweep-plane-3d-v5. Canonical READ-ONLY (this notebook is the
deliverable). No experiments needed - both formulations resolve by code-read +
theory against the grounded corpus evidence. ASCII; file:line for code claims;
[SPEC] marks speculation.

Mission (owner-commissioned, third act): the exact kernel (Orient3DExactSign +
SoS) is the proof-of-concept; the destination is eliminating it via a
Smith-Dodgson third act (approximate arithmetic, topological guarantee,
eps-geometry). I own TWO formulations distinct from the sibling sd-self lane
((i) direct seam-graph telescoping, (ii) convention+cocycle-repair - NOT
touched here):

  1. TWO-COPY BIPARTITIZATION - reduce self-overlap to Smith's bipartite boolean.
  2. LAYERED-CONSISTENCY - assemble the proven-free pieces, find the residue.

Grounding read: src/boolean3.cpp (Shadow/expandP/Winding03), docs/Boolean2.md
(Smith ch.7 block rule), docs/Regularize3D.md (witness thm, coupled winding,
cocycle), orient-attack-1784086942.md (WindingAt is the sole corpus-irreducible
exact winding caller; da seed-invariant; only a flood re-arch avoids it),
nomerge-1784160552.md (identity is level-0 free, aliasGap=0), attack22-c3-
1784146879.md (radial order level-0 free from input normals).


===============================================================================
FORMULATION 1: TWO-COPY BIPARTITIZATION
===============================================================================

## 1.1 What Smith's bipartite structure actually requires (code-read)

boolean3 keys its symbolic perturbation PER OPERAND, not per face or per pair.
`expandP_ = (op == OpType::Add)` is a SINGLE bool for the whole operation
(boolean3.cpp:475), threaded through the shadow cascade as
`withSign(expandP, a0xp) - b1exp` (boolean3.cpp:63-66): every A-feature is nudged
by +expandP, every B-feature by the base. Shadows(p,q,dir) breaks the p==q tie by
sign(dir) (shared.h:121-122). So the operand label is a property of the FACE
(via its operand), uniform across ALL that face's crossings.

The crossing enumeration is bipartite by CONSTRUCTION: Intersect12<true> finds
P-edge x Q-face, Intersect12<false> finds Q-edge x P-face (boolean3.cpp:508,514).
There is no P-edge x P-face pass - a clean operand has no self-crossings. That
absence IS the bipartite hypothesis. Winding03 then unites vertices joined by an
UNBROKEN edge (an edge in no crossing pair, boolean3.cpp:407-412), ray-casts ONE
seed per resulting component (Kernel02, :439-448), and floods within
(:451-456). Crucially the seed ray asks "winding of this P-vertex INSIDE the
Q-SOLID" - Q must be a closed oriented surface bounding a region.

Two non-negotiables fall out: (N1) each operand is a CLOSED surface; (N2) each
face carries ONE operand identity, driving ONE perturbation direction.

## 1.2 The min(i,j) per-pair convention: acyclic, but the wrong type

The proposal assigns the P-role to min(i,j) at each self-crossing. Test the
prompt's hypothesized failure (role cycles along a seam loop): a cycle needs
i<j, j<k, k<i - IMPOSSIBLE, since "<" on integer indices is a total order. The
min-index convention is ALWAYS acyclic. So role cycles are NOT the obstruction.

The real defect is a type mismatch. The role is a property of the PAIR {i,j},
not of the face. Face i is P when paired with j>i and Q when paired with h<i, so
it plays BOTH roles. This violates N2 directly: at a vertex shared between the
i-j seam and the h-i seam, the per-operand perturbation would nudge face i in
two contradictory directions at once. And it violates N1: "P = faces that are
the min of their pair" is not a fixed face-set at all (membership depends on the
partner), so there is no P-surface, closed or otherwise, hence no "inside Q" for
the winding seed. The convention dodges cyclicity precisely by refusing to be a
face-partition - and that refusal forfeits the two operands.

## 1.3 The sharper test: does a GLOBAL operand 2-coloring exist?

The honest strengthening: forget per-pair; ask for a fixed face->{P,Q} coloring
such that every self-crossing is P-Q (then each face has one operand identity,
satisfying N2). This is exactly: is the SELF-CROSSING GRAPH bipartite? (vertices
= faces, edges = crossing pairs from RecordSeams, overlap3.cpp:1383).

CONCRETE OBSTRUCTION. Three faces i,j,k that pairwise cross - a triple
self-incidence, three sheets meeting near a common line - form an odd 3-cycle
i-j-k in the crossing graph. An odd cycle is NOT 2-colorable, so NO global P/Q
coloring exists. The two-copy reduction is obstructed at the triple point.

This is not a fringe case: it is EXACTLY the campaign's unbuilt >2-sheet radial
junction (docs/Regularize3D.md:344-348, Cluster 3). The structure that defeats
Smith's bipartite reduction is the identical structure the exact pipeline gates
behind the kernel tripwire. Two-copy bipartitization and the radial branch are
the SAME obstruction viewed from opposite sides: "3 sheets at an edge" is
non-bipartite as a crossing graph and >2-sheet as a radial fan.

## 1.4 Even WHERE the graph is bipartite, the pieces are not operands

Suppose the crossing graph happens to be bipartite (no triple incidences) and we
2-color into P_set, Q_set. Neither is a closed surface: P_set is an arbitrary
subset of S's triangles, and S's shared manifold edges are split between the
colors, so P_set has boundary. N1 fails again. Bipartiteness kills only P-P/Q-Q
crossings; it says nothing about closure. So the reduction is DOUBLY obstructed:
non-bipartite at triple points (1.3), and open-pieced even when bipartite (1.4).

## 1.5 The winding algebra: the identity exists, but it is trivial

"Recover w_S from P-Q crossing parities (w_S vs the w_P,w_Q algebra)":

Winding is additive over ANY face partition - a ray counts each face's signed
crossing independently. So for any coloring,
    w_S(x) = sum over faces delta_f(x) = w_{P_set}(x) + w_{Q_set}(x).
This identity holds trivially and is precisely what WindingAt already computes:
w += delta per triangle, delta = sign(dot(seed-p, n)) (overlap3.cpp:1116,1236).

But w_{P_set} is NOT the boolean's w_P. The boolean's w_P is the integer cover of
a CLOSED SOLID (how many times P's region contains x); w_{P_set} is the winding
of an OPEN surface, not integer-valued as a solid cover and not thresholdable by
a Smith fill rule (Add: w>0, Intersect: w>1, docs/Boolean2.md:119-126). So the
Smith solid-winding algebra (w_P, w_Q as operand covers) does NOT EXIST for the
self case; only the trivial per-face additivity does, and that adds nothing the
resolver's per-face delta sum does not already do.

## 1.6 The salvage (the constructive best-of)

The bipartite idea has ONE surviving fragment, and it is already landed. The
value of an operand label was never the surface split - it was a GLOBAL TOTAL
ORDER that breaks perturbation ties consistently (Shadows' dir, keyed by operand
in the boolean). The self-overlap analog replaces "operand" with "global vertex
index": Edelsbrunner-Mucke SoS keyed by global vertex index, never 0 for
distinct indices (orient-attack:16; docs/Regularize3D.md:233-238, stage 6). The
min(i,j) intuition is CORRECT at the tiebreak layer (a global index total order
decides consistently) and WRONG at the surface layer (no two operands). The
campaign already banked the correct half: SoS-by-index IS "min-index lifted from
operand-role to vertex-role," and it is total by construction (no window-fail).

## 1.7 Formulation 1 verdict

OBSTRUCTED (doubly). (1.3) The global operand 2-coloring exists iff the
self-crossing graph is bipartite; a triple self-incidence is an odd 3-cycle =
the campaign's >2-sheet radial junction, so no coloring exists there. (1.4) Even
when bipartite, the color classes are open surfaces, not Smith operands. (1.5)
The winding "algebra" collapses to trivial per-face additivity (= WindingAt),
not the solid-winding w_P+w_Q. (1.6) The only transferable piece - the global
total-order tiebreak - is already the landed SoS-by-index. Contributes NO new
residue: its obstruction IS the already-named radial branch.


===============================================================================
FORMULATION 2: LAYERED-CONSISTENCY
===============================================================================

Assemble the campaign's proven-free pieces into a full no-exact design and find
what is LEFT that still needs a metric SIGN approximate arithmetic can get wrong
non-realizably.

## 2.1 The assembled design (each layer with its code anchor + proof pointer)

(a) IDENTITY = provenance + equality weld. Proven level-0: aliasGap=0 (provenance
    never absent), every aliasing pair decidable by an input-coordinate incidence
    predicate, not distance (nomerge-1784160552.md:90-110, 155-161). Sites:
    RecordSeams dedup is per-pair provenance-consistent (overlap3.cpp:1383+); the
    overlay/weld metric merges are REPLACEABLE by level-0 incidence.
(b) PER-SEAM CANONICAL 1D ORDER computed once in 3D. Total orders are always
    realizable (SoS-breakable); all incident faces consume the SAME order because
    a triple point is pre-split ONCE and threaded into every incident face's
    overlay bit-identically (overlap3.cpp:1947-1954, Intersect3Planes symmetric
    key :1957-1972).
(c) PER-FACE 2D OVERLAY = Smith's block rule with (b) as input. This is literally
    the built path: EmitSeamedFace calls RemoveOverlaps2D per crossed face
    (overlap3.cpp:2492, 2737), the boolean2 sweep + block rule (docs/Boolean2.md:
    38-52, section 7.6.2), on EXACT axis-drop 2D coords (no projection rounding,
    docs/Regularize3D.md:654-667, L1 resolved).
(d) WINDING BY FLOOD, no seed rays. Seed one extremal cell (w jumps 0->1 at the
    lexicographically extremal vertex, combinatorial), propagate integer deltas
    +-1 across the arrangement's cell-adjacency graph. This REPLACES the ray-cast
    WindingAt (overlap3.cpp:1226-1238) whose da graze is the sole corpus-
    irreducible exact fire (orient-attack:107-113,133-140).
(e) SoS ONLY for exact ties (structural coincidences), realizable by index
    (Formulation 1.6).

## 2.2 Suspect 1 - cross-face block-grouping agreement at triple points

Smith's block rule (docs/Boolean2.md:44, section 7.6.2) groups near-concurrent
events by eps-proximity ON ONE SWEEP LINE - a PER-FACE, PER-SWEEP, METRIC
decision. Three faces sharing a triple point each run their OWN RemoveOverlaps2D
in their OWN axis-drop projection, so the SAME two 3D crossings can be within-eps
on face A's projection and beyond-eps on face B's (different dominant axis,
different 2D metric). If A groups {s1,s2} and B splits them, the shared seam s1's
sub-edge subdivision disagrees and the assembly weld fans open.

Is there a CANONICAL GLOBAL blocking? As a metric operation across faces: NO -
eps-proximity is projection-dependent. It EXISTS only if every triple point is
pre-constructed exactly (once-only) and handed to each face's overlay as an
EXPLICIT shared vertex, so the block rule never has to eps-rediscover it. That is
exactly what B1 does when the crossing is level-0 enumerable (overlap3.cpp:
1947-1954): pre-split, symmetric key, bit-identical to all three faces ->
grouping is canonical BY CONSTRUCTION, no per-face eps decision survives. So:

  - Triple points the enumeration RESOLVES (level-0 on input coords): pre-split,
    canonical, COMBINATORIAL-SAFE. This is the whole resolving corpus (zero
    genuine multi-sheet edges, docs/Regularize3D.md:185).
  - Triple points the enumeration REFUSES (near-tangent, phantom-seam guard,
    overlap3.cpp:1714-1744): NOT pre-split, left to the per-face eps block rule,
    which groups differently across faces / collapses a cluster spanning two
    winding regions into one cell. This is the openscad E1 wall EXACTLY (near-
    tangent triple points from ~1e-4 to ~2eps that RemoveOverlaps2D's eps-merge
    collapses, docs/Regularize3D.md:561-573).

VERDICT: GENUINE METRIC SIGN, but it is not a NEW residue - it is E1. The
"canonical global blocking" the prompt asks about is the exact per-face
near-tangent 2D arrangement over CONSTRUCTED crossings (the priced terminal): the
one tripwire crossing is orient2d on a constructed near-tangent point
(docs/Regularize3D.md:581-590; attack22-c3:104-119).

## 2.3 Suspect 2 - the flood's delta signs at near-tangent crossings

The re-architecture (d). Stepping from cell A across arrangement-face f to cell B,
delta = +1 if the traversal enters f's +n_f half, else -1. Decompose:

  - WHICH face is crossed: NOT re-derived by a sign test. The cell-adjacency
    graph RECORDS that A and B share f (combinatorial). Contrast the ray-cast,
    which asks da = Orient3DFilterSign(a,b,c,p) "is p above face t's plane" for
    EVERY t (overlap3.cpp:1100) - that is the near-tangent graze that fires exact
    on GT7081's 0.002deg twins (orient-attack:88-89), and it is SEED-INVARIANT so
    no ray/seed tweak clears it (orient-attack:133-140). The flood NEVER computes
    a da: it only crosses faces the adjacency already says it crosses.
  - THE SIGN once f is known: sign of n_f's component along the A->B step. n_f is
    an INPUT triangle normal (exact, never near-zero - near-tangency is a small
    ANGLE BETWEEN two faces, not a small normal OF either). The +n_f-side label of
    f's two incident cells is computed ONCE at build time (a level-0 input-normal
    sign, never near-zero) and stored. So the delta is a STORED COMBINATORIAL BIT
    plus an exact input-normal read.

VERDICT: COMBINATORIAL-FROM-PROVENANCE (safe). No metric threshold, no
near-tangent fire. This is the mechanism by which the flood removes the sole
corpus-irreducible exact winding caller (WindingAt/GT7081) - it converts "which
side of a near-tangent face is the probe on" (metric) into "cross a face the
adjacency already names" (combinatorial). Confirms orient-attack's "only a
coupled-integer-flood WIND re-architecture avoids them" and shows WHY.

COUPLING (the honest cost). The flood is well-defined only if the cell complex is
CLOSED - path-independence (the cocycle, docs/Regularize3D.md, coupled winding
deltas + path-independence) holds iff every arrangement face borders exactly two
cells. At the E1 open-boundary fan the complex is NOT closed, two paths to one
cell give different w, and the flood is ill-defined. So the delta SIGN is safe,
but the flood's WELL-DEFINEDNESS is gated on the SAME arrangement completeness as
E1. The cocycle is the formal statement of that gate. You do not get the flood's
exact-freedom until you have paid E1's completion; once paid, the winding half is
exact-free.

## 2.4 Suspect 3 - fold / coplanar cluster membership

  - EXACT coplanar fold: two faces coplanar iff their planes are identical =
    Orient3DExactSign of one face's verts against the other's plane == 0 (level-0
    on input coords, realizable by SoS). COMBINATORIAL-SAFE.
  - NEAR-coplanar membership ("gap below eps"): a genuine METRIC threshold
    (SnapNearCoplanarClusters, the smaller directional vertex-plane max < eps,
    docs/Regularize3D.md:210-218). BUT it is an INPUT-conditioning decision, not
    an arrangement-topology sign: it either snaps to EXACTLY coplanar (then the
    exact fold runs, realizable) or FAILS CLOSED at the global-planarity guard
    (:214-218, S3 census). Never a silent wrong sign on the arrangement.

VERDICT: exact membership combinatorial-safe; near-coplanar membership metric but
CONTAINED (snap-to-exact or fail-closed guard - a decision-complete refusal,
docs/Regularize3D.md:604-609, not a residue that emits wrong geometry).

## 2.5 Formulation 2 verdict - ONE genuine residue

Assembling (a)-(e), the residue-per-junction:
  identity (a)          -> level-0 free (nomerge)                    SAFE
  1D seam order (b)      -> total order, always realizable           SAFE
  block overlay (c)      -> Smith's proven second act                SAFE given (b)
  flood delta (d)        -> combinatorial (Suspect 2)                SAFE
  SoS ties (e)           -> realizable by index (Formulation 1.6)    SAFE
  radial order at nodes  -> level-0 from input normals (attack22-c3) SAFE
  block GROUPING at      -> metric IFF near-tangent triples unsplit  = E1
    triple points (S1)      (canonical when pre-split once-only)
  coplanar membership(S3)-> exact safe / near contained             fail-closed

The design has EXACTLY ONE genuine metric-sign residue: the EXACT PER-FACE
NEAR-TANGENT 2D ARRANGEMENT over CONSTRUCTED crossing points (E1 / Cluster 1).
Its tripwire is orient2d on a constructed near-tangent crossing (off the level-0
input-coordinate discipline - the coordinates are built intersection points, not
input). Everything else is combinatorial, realizable-order, or contained
fail-closed. And the flood (d) STRICTLY IMPROVES on the current pipeline: it
removes the winding-probe exact caller (WindingAt), leaving the E1 constructed
sign + the realizable index-SoS as the only exact surface - but only after the
same arrangement completion E1 requires (the cocycle coupling, 2.3).


===============================================================================
MIGRATION-PATH SKETCH (exact kernel as scaffolding + differential oracle)
===============================================================================

The exact pipeline stays as the DIFFERENTIAL ORACLE throughout: every no-exact
step is validated byte-or-oracle against the exact resolve on the corpus before
the exact path is retired for that step. Order (least to most coupled):

  STEP A [survives as-is]. Enumeration level-0 predicates + SoS-by-index for
    structural ties (Formulation 1.6). Already the landed design; keep. The
    index-SoS is the bipartite-role salvage and is realizable, so it is NOT a
    kernel-elimination target - it is the intended terminal for exact TIES.

  STEP B [replace: winding]. Swap ray-cast WindingAt (overlap3.cpp:1226-1252)
    for the combinatorial flood over the cell-adjacency graph (Suspect 2). This
    is the single highest-value swap: it deletes the sole corpus-irreducible
    exact winding caller (WindingAt/GT7081, orient-attack:151-153). GATE: the
    flood needs the complete cell complex (cocycle, 2.3), so B cannot land before
    the arrangement is closed. Until then, run flood and ray-cast in parallel and
    diff (the exact ray-cast is the oracle); they must agree on every closed
    corpus component. [SPEC] the flood likely lands first on the already-closed
    corpus (siA/siB/GT7863/Bridged/Bars) and stays fail-closed-or-oracle on
    openscad, matching today's boundary.

  STEP C [survives, hardened]. Identity: replace the overlay/weld METRIC merges
    (nomerge sites 2/4) with level-0 incidence + equality weld (nomerge:112-153).
    Free on resolving carriers (weldConflate=0 by once-only, nomerge:65-73);
    net-new work is once-only construction plumbing, not a decision. Oracle: the
    exact resolve's vertex identity.

  STEP D [the terminal - the only genuine kernel target]. E1 / Suspect 1: the
    exact per-face near-tangent 2D arrangement over constructed crossings. This
    is where the exact kernel is genuinely ELIMINATED vs merely INDEX-REALIZED:
    carry each near-tangent crossing SYMBOLICALLY as its plane-triple {f,g,h} and
    evaluate orient2d via degree-bounded coefficients on the EXISTING
    Orient3DExactSign accumulator - ONE new predicate FORM, zero vendored code
    (~150-250 LOC, attack22-c3:113-119; docs/Regularize3D.md:581-590). Completing
    this closes the cell complex, which UNBLOCKS the Step B flood on openscad.
    Oracle throughout: the exact pipeline's own emission.

Survives unchanged: decompose/gate/compose, the coplanar EXACT fold, the block
rule (boolean2), radial order (level-0), SoS-by-index for ties. Replaced: ray-
cast winding -> flood (B), metric identity merges -> level-0 incidence (C).
The exact e^0 orient3d as a FUNCTION never fully disappears - it survives as the
tie-cascade's realizable index-SoS (Step A) and as the {f,g,h} FORM's evaluator
(Step D). "Kernel-free" is achievable for the WINDING and IDENTITY halves; the
ENUMERATION half keeps a realizable exact tie-breaker by design, not by defeat.


===============================================================================
RESIDUE SUMMARY
===============================================================================

Formulation 1 (two-copy bipartitization): OBSTRUCTED (doubly). Contributes 0 new
  residue - its obstruction (odd-cycle self-crossing graph at a triple
  incidence) IS the already-named >2-sheet radial branch. Salvage: the global
  total-order tiebreak survives as the landed SoS-by-index.

Formulation 2 (layered-consistency): SOUND design, 1 genuine metric-sign residue.
  RESIDUE-1 = E1 / exact per-face near-tangent 2D arrangement over CONSTRUCTED
  crossings (block-grouping at unsplit near-tangent triples; orient2d on a
  constructed point). All other layer-junctions resolve SAFE: identity (level-0),
  seam order (realizable), block overlay (Smith proven), flood delta
  (combinatorial), radial order (level-0), SoS ties (realizable), coplanar
  membership (exact-safe / contained fail-closed).

CONVERGENCE: both formulations, from opposite directions, land on the SAME single
wall the exact pipeline already names - the constructed near-tangent 2D
arrangement (E1 / Cluster 1). Formulation 2's flood additionally REMOVES the sole
corpus-irreducible exact WINDING fire, but is coupled (via the cocycle) to the
very arrangement-completeness E1 demands. Net: the honest no-exact destination is
kernel-free winding + kernel-free identity, gated on ONE terminal (E1) whose
exact surface is a single {f,g,h} orient2d FORM, not a general kernel.
