# A0 VERIFY / AMBIGUITY-WITNESS lane notebook

Branch: explore/sweep-plane-3d-v4 @ dd4e7119 (canonical READ-ONLY;
work in /tmp/a0-verify-witness copy).
Started: 2026-07-13.
Mission: the OPPOSITE direction from the other A0-verify lane. That lane asks
"did the probe miss retained data." THIS lane asks whether the INPUT ITSELF
determines the oracle answer at all. Construct an AMBIGUITY WITNESS: two operand
pairs (A,B) and (A',B') with byte-IDENTICAL composed soups but MACROSCOPICALLY
different true Boolean unions. If it exists -> fail-closed upgrades to "NO
algorithm consuming this input can match the oracle without guessing" (proof-
shaped closure). If every construction fails -> characterize the invariant that
forces uniqueness (the future research direction), RSI-#3 style.

## Step 0: grounding (the pipeline's actual input + fill rule)

- ComposeImpl (test/overlap3_test.cpp:47) = pure ORIENTED-FACE CONCATENATION:
  it appends operand B's vertProperties + (triVerts+base) onto A's, one Impl.
  No boolean, no reorientation. The pipeline (RemoveOverlaps3D) sees exactly
  this soup and nothing about the A/B split (confirmed: Canonicalize folds by
  signed mult, discards operand label - closure lane FACT 1).
- FILL RULE. The pipeline emits the boundary where the winding number crosses
  the >= 1 threshold. Confirmed verbatim by the NestedCubes fixture comment
  (test/overlap3_test.cpp:205-212): "The inner shell is entirely inside:
  winding=2 on both sides -> not a fill transition -> not emitted." So the
  pipeline computes {p : w(p) >= 1} for the soup, w = signed winding.
- ORACLE. The tests use `a + b` = Boolean3 union (e.g. :710). For two
  positively-oriented closed solids this is the set-union {p in A or p in B}.

## The proof sketch that frames the whole lane (INPUT-DETERMINES candidate)

Let S be a fixed ORIENTED face soup (a multiset of oriented triangles). Define
the signed winding field w_S(p) = sum over faces f in S of the signed solid-
angle / ray-crossing contribution of f at p. w_S is a PURE FUNCTION of S
(nothing else): the winding at a point is determined by the oriented faces.

Union regime (the A0 targets: every canonical face mult +1, all outward). If S =
compose_union(A,B) = faces(A) + faces(B) with A,B positively-oriented closed
manifolds, then w_S = winding_A + winding_B, and BOTH terms are >= 0 everywhere.
Hence
    A union B = {p : winding_A(p) >= 1 OR winding_B(p) >= 1}
              = {p : winding_A(p) + winding_B(p) >= 1}       (both terms >= 0)
              = {p : w_S(p) >= 1}.
Every quantity on the last line depends ONLY on S. So for ANY decomposition of a
fixed non-negative-winding soup into positively-oriented closed manifolds, the
set-union is exactly {w_S >= 1} - a PURE FUNCTION OF THE SOUP. Partition-
equivalence is FORCED. The witness CANNOT exist in the union regime.

The ONLY escape is an operand whose winding goes NEGATIVE somewhere (an inward /
"hole" shell). But a negative winding requires INWARD-oriented faces in S -> a
DIFFERENT oriented soup. So the escape is excluded by the "same soup" premise.
(Empirically pinned as Case C below: flipping an operand's orientation changes
the canonical soup multiset AND the union.)

Cross-operation corollary. The pipeline computes {w_S >= 1} and does not know
the "operation". If one soup S were reachable as BOTH a union-compose and a
subtraction-compose, both would still have oracle result {w_S >= 1} - the SAME
set - because subtraction bakes B's flip INTO S (w_S = winding_A - winding_B is
already the subtraction field). There is no conflict; the soup's {w_S >= 1} is
simultaneously correct for every operation that produced it. So no cross-op
witness either. The only way to break determinism would be a DIFFERENT fill rule
(even-odd, >=2); the pipeline and the positive-fill Boolean both use >= 1.

## Step 1: the probe (test/overlap3_test.cpp Overlap3.WitnessProbe, env-gated
OV3_WITNESS, THROWAWAY). Builds operands from raw exact-coordinate triangle
primitives (per-sign fixed diagonals so a coincident wall is 2-manifold),
deduping verts EXACTLY; compares composed soups via an orientation-preserving
canonical triangle multiset with EXACT double compare (house rule); oracle =
library `a + b`; also runs the composed soup through RemoveOverlaps3D. Fresh
cmake vbuild (Release, MANIFOLD_TEST=ON) - built clean, ran in 3ms.

### Case A - COMPONENT-RELABEL (two VALID non-degenerate splits, IDENTICAL soup)
Soup = two disjoint unit cubes. Split1 = (P, Q). Split2 = (Compose{P,Q} as ONE
2-component manifold, EMPTY). Both operands valid (Status ok).
  MEASURED: composed soups BYTE-IDENTICAL (canonical multiset EXACT-EQUAL=1).
  oracle union(P,Q) = 2.0000  ==  union(PQ, Empty) = 2.0000.
This is the constructive witness-ATTEMPT-that-fails: two genuinely different
operand stories, byte-identical soup, IDENTICAL union. The invariant holds
against the real library. (The only non-degenerate multi-split freedom is moving
whole connected components between operands - provably union-invariant.)

### Case B - COINCIDENT-WALL repartition (the crux the design worried about)
Split1 = L[0,1]^3 + R[1,2]x[0,1]^2 sharing an exact wall at x=1 (24 tris,
including the buried coincident +x/-x double-wall). Split2 = Big[0,2] (wall
interior, no x=1 faces) + DoubleWall (the coincident sheet as its own operand).
  MEASURED: L,R,Big all valid (vol 1,1,2). DoubleWall built ok but vol=0, and
  Manifold's ingest CLEANED IT AWAY: split2 soup = 20 tris (the 4 x=1 wall tris
  are the whole difference; verified all-x=1-plane tris present in split1,
  ABSENT in split2). So the two soups are NOT byte-identical (24 vs 20).
  BUT oracle union(L,R) = 2.0000  ==  union(Big,DoubleWall) = 2.0000, and the
  PIPELINE resolves BOTH to vol 2.0000.
FINDING: I could not make the coincident-wall soups byte-identical, and the
reason REINFORCES input-determines: (a) the clean-operand partition of a
coincident-wall soup is essentially UNIQUE (a +x wall at x=1 must bound a solid
to its left, forcing L=[0,1]; likewise R) - so there is no alternative clean
split to begin with; (b) the ONLY achievable variation is burying-vs-splitting
the coincident double-wall, which is WINDING-NEUTRAL (the two opposite coincident
faces cancel in w everywhere off the measure-zero wall), and Manifold cleans it;
(c) union and pipeline output are invariant (2.0) across that variation. The
oracle tracks {w>=1} EXACTLY even through coincident faces - no special-casing
deviation.

### Case C - the ORIENTATION ESCAPE (why the soup pins the union)
Body[0,3]x[0,3]x[0,1] + sliver[1,2]x[1,2] standing z in [1,2]. Sliver emitted
OUTWARD (survives as a protrusion) vs INWARD (all faces flipped).
  MEASURED: sliverOut valid vol=+1.0; sliverIn valid vol=-1.0 (Manifold ACCEPTS
  the inward shell as a closed 2-manifold). Soups DIFFER (soupEQ(out,in)=0 - the
  flip changes every triangle's winding key). Unions DIFFER:
  union(body,sliverOUT) = 10.0  vs  union(body,sliverIN) = 8.0.
This is the proof's escape clause made concrete: the ONLY way to change the union
is to change the soup, and negative winding requires inward faces = a different
oriented soup. Within a FIXED soup there is no orientation freedom left to
exploit. (Note: a negatively-oriented operand is where a SUBTRACTION-composed
input's -1 mults would live - out of the A0 all-+1 union scope, and still
soup-determined by the same argument with the flip baked into S.)

### Pipeline symmetry
RemoveOverlaps3D on split1 (24-tri, buried wall) and split2 (20-tri, clean box)
both -> resolve vol=2.000000. Same {w>=1} answer; the pipeline is invariant to
the winding-neutral soup difference, as its fill rule predicts.

## VERDICT: INPUT-DETERMINES (no ambiguity witness exists)

For the union regime A0 measured (and, with the flip baked into S, every regime
the pipeline sees), the oracle result is a PURE FUNCTION of the oriented face
soup: result = {p : w_S(p) >= 1}, w_S determined by S alone. Proved (necessity
of orientation, non-negativity of union windings) and confirmed against the real
library oracle across three construction families - component-relabel (identical
soup -> identical union), coincident-wall (partition rigid; only variation
winding-neutral), orientation-flip (changes the union ONLY by changing the soup).
Every witness construction collapses to the invariant. No two operand splits with
the same soup have different unions.

CONSEQUENCE for the fail-closed adjudication (the point of the lane). The
NOT-SEPARABLE fail-closed does NOT upgrade to a proof-shaped "no algorithm can
match the oracle" impossibility. This lane proves the OPPOSITE: because the input
DETERMINES the answer, an EXACT-arithmetic resolver CAN match the oracle
(evaluate {w_S >= 1} exactly). A0's NOT-SEPARABLE is a statement about the
RETAINED, EPS-QUANTIZED identity table (arr.verts/seams/mult), which quantizes
away the sub-eps coordinate distinctions that carry w_S near the degeneracy. The
escape kill's lost 336-lump is therefore recoverable IN PRINCIPLE from the exact
soup (the lump has real thickness -> winding 1 in S); the pipeline loses it only
because eps-tolerant collapse discards the distinguishing sub-eps geometry. So
the barrier is PRECISION, not INFORMATION.

This CONFIRMS the docs' named direction rather than closing it off: the research
frontier is EXACT arrangement completion (RSI-#3 / MaintainedEmission variant
iii), which is precisely the one channel that carries the input's full winding
information through emission. The honest fail-closed posture for the CURRENT
eps-machinery stands; it is not a fundamental information barrier.

SCOPE / caveats. (1) Manifold's ingest cleans winding-neutral coincident content
from a single operand (Case B), so a byte-identical-soup pair can only be built
where the differing content is winding-neutral anyway - which is exactly why no
witness survives. (2) Subtraction inputs carry -1 mults (Case C's inward shell);
still soup-determined, out of A0's union scope. (3) The proof assumes the >= 1
positive-fill rule, which both the pipeline (NestedCubes comment) and the
library's union use; a different fill rule (even-odd, >=2) is not in play.

## WRAP / fence
- Probe REVERTED in the /tmp copy (git checkout test/overlap3_test.cpp); the copy
  is the only thing touched. Canonical repo untouched except this notebook.
- vbuild is throwaway under /tmp/a0-verify-witness.
- No git commit, no push. Notebook-only deliverable.
