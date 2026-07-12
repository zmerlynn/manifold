# WALL-A tier-3 memo - adversarial verification lane (epistemics)

Artifact: docs/SweepEmit3D.md WALL-A section (## WALL-A, lines 1255-1395) +
notebook .claude/lane-reports/walla-arc-1783880487.md. Repo HEAD b68c6ea4
(explore/sweep-plane-3d-v3). Job: find the hole in a tier-3 "irreducible
coupling" memo - exhaustiveness of candidate space + conclusiveness of kills.

Method: read memo + notebook + cited sections (PROVENANCE CHAINS, 3D-IDENTITY
EXTENSION, CHAIN-PLANE RULE, STRICT-FP, WALL-B) + RSI-#3 template
(memory/rsi-3-completion-research-grade.md; the .claude/plans path in the brief
does not exist, memory holds the template). Traced the emission path in
src/overlap3.cpp / overlap3_sweep.cpp / boolean2.h. Re-read the arc's OWN
surviving drivers+logs in this session's scratchpad (git diff src/ empty =
canonical HEAD, instrumentation reverted, confirmed by grep: no OV3_THIN/FANDBG/
getenv in src). No heavy re-runs needed - the strongest hole is a code-trace +
a contradiction against the arc's own fan dump.

## VERDICT: NEED-CHANGE

The irreducible-coupling CONCLUSION survives (no bounded fix actually works: the
C3 any-axis weld-bump to 4x eps did not resolve Havocglass8; the fan's
twin-merge genuinely fails). But the memo under-earns tier-3 on TWO points: an
axis-inverted anatomy that mis-targets a kill, and an unnamed candidate. Neither
is a BREAK.

## STRONGEST HOLE: the anatomy inverts the twin's separation axis, mis-targeting C4

MEMO CLAIM (anatomy part 1, docs line 1283-1286): the fan twin is "the
3D-IDENTITY EXTENSION wall's near-degenerate arr.verts cluster (one 3D junction
= two verts kept distinct by Canonicalize because they are more than eps apart
in x while sub-eps in the transverse plane)".

ARC'S OWN FAN DUMP (scratchpad/fan_dump_havoc.log, Havocglass8, current HEAD)
shows the INVERSE for the actual failing twin. The TIE micro-edge v4-v509:
  v4  = (12942.044899999999, -549.12713600000006, 2009.0866700000008)
  v509= (12942.044899999757, -549.12713599066421, 2009.0866700370959)
  dx = 2.42e-10 = 0.0107 eps   (SUB-eps in x)
  dy = 9.34e-9  = 0.415 eps
  dz = 3.71e-8  = 1.648 eps    (SUPER-eps, dominates)
  |3D| = 1.698 eps  (dump: len/eps=1.698)
So the twin is SUB-eps in x, SUPER-eps in the transverse (yz) plane - the exact
opposite of the memo's ">eps in x, sub-eps transverse".

CORROBORATION from the cited 3D-IDENTITY EXTENSION section + identity-extension
notebook: there are TWO distinct near-degeneracies, and the memo grafts one onto
the other:
  - CASE (2) DISTINCT pair (identity-ext nb lines 161-166): V1,V2 = two genuine
    canonical verts 72.8 eps apart in X, sub-eps in yz; "Canonicalize did NOT
    merge them (>eps in x)". THIS is what C4 was killed on.
  - CASE (1) SPREAD (identity-ext nb 192-193, 234, 326-327): ONE vertex V whose
    cap critical sits 0.3-4 eps (SUB-eps) from V.x; a steep track amplifies that
    sub-eps x-offset into a large yz-divergence. THIS is the fan twin.
The fan (case 1) twin's two images are NOT "two arr.verts kept distinct by
Canonicalize" - they are two EMITTED cap-plane images of ONE arr.vert V, from
two cap planes 0.0107 eps apart in x (v4.x != v509.x; TriangulateCap emits every
cap vert at exactly xCap, so distinct x = distinct critical = distinct cap),
kept distinct by the ASSEMBLY eps-weld (BuildImpl, radius eps), not by
Canonicalize. capEps=8*eps never saw them together because they are in SEPARATE
cap arrangements; they meet only at the assembly weld where 1.7 eps > eps.

CONSEQUENCE (the kill is mis-aimed). C4 (ANISOTROPIC CANONICALIZE / yz-merge) is
KILLED in the memo (line 1327-1332) entirely on case-(2) grounds: "the
near-degenerate arr.verts PAIR is a real MACRO separation in x (order tens of
eps)... merging them in the transverse plane collapses real x-extent". For the
ACTUAL fan twin (case 1) there is NO macro x-extent to collapse - x-separation
is 0.0107 eps. A reader who identifies the fan twin as sub-eps-in-x and proposes
a yz-merge (or a near-x critical-plane merge letting one capEps=8*eps
arrangement absorb the 1.7-eps transverse spread) finds C4's stated kill FALSE
for their configuration. The candidate is real: merging the two cap planes
0.0107 eps apart would put both extensions in ONE arrangement and capEps=8*eps
(> 1.7 eps) would merge them - closing THIS fan locally. It is neither in the
killed list under its true signature nor covered by C2 (which touches only
crossings, never vertices - and the notebook confirms crossing-thinning leaves
the fan untouched, so the fan-driving criticals are VERTICES, exactly what C2
excludes).

The correct kill EXISTS but is not routed here: the fan twin-merge is killed by
residual-1 cluster-collapse (3D-IDENTITY EXTENSION section) + the C3 any-axis
weld-bump evidence (a yz-merge is a subset of the 4x-eps weld-bump that failed
on Havocglass8). So NEED-CHANGE, not BREAK: separate the case-(1) spread twin
(sub-eps in x) from the case-(2) distinct pair (>eps in x) in the anatomy; scope
C4's kill to case (2); route the case-(1) fan yz/critical merge kill explicitly
to residual-1 / C3.

## ATTACK 1 - unprobed candidate hunt (per (a)/(b)/(c)/(d))

(a) PROVENANCE-IDENTITY TWIN MERGING - CODE TRACE, then RELOCATES; UNNAMED.
  Traced: OutTri3D = {vec3 v[3]} (overlap3.cpp:603) - positions only, no
  identity. edgeSubdiv fed to strip chains = vector<vector<vec2>> (positions;
  boolean2.h:225-236, ComputeCap:853). Cap arrangement RemoveOverlaps2D
  (ComputeCap:855) passes NO edgeClass/classSubdiv - RemoveOverlaps2D's public
  signature does not even take them (they live on internal SweepWinding,
  boolean2.h:207-212). EdgeM.srcId is a FACE id and is UNCONSUMED in 3D (WALL-B;
  SlabResolver keys by POSITION). BuildSlabs' section call
  (overlap3_sweep.cpp:230) also passes no classSubdiv. SlabResolver::Track =
  {bool weld; vec3 a,b} (positions). BuildImpl weld = eps hash grid over vec3
  (positions, radius eps). => provenance is POSITIONS-ONLY end-to-end; no
  cross-cap 3D-junction identity exists in the emission path. So the candidate
  is "thread arr.vert(V) identity through slab->cap->edgeSubdiv->OutTri3D->weld,
  merge twins by label-equality." VERDICT: RELOCATES (RSI-#3). Unlike
  negEdges/classSubdiv (which fixed mis-PAIRING / mis-SUBDIVISION with ZERO vertex
  move), closing the fan requires collapsing a POSITION divergence (1.7 eps) -
  the label changes only the SELECTION criterion (robust to the >eps gap), but
  the merge ACTION still moves a vertex 1.7 eps => degenerate/doubled cap sheets
  => residual-1 re-emission. The memo's "no output weld can merge the images
  without collapsing load-bearing pieces" (line 1377-1379) substantively covers
  a label-weld (same action, different selection), so the substance is present -
  but the candidate is NOT NAMED and the positions-only code trace is absent.
  For a tier-3 memo written right after the provenance-chains machinery landed
  ("the necessary subdivision foundation"), a reader's most natural next idea is
  "merge by provenance identity." Leaving it implicit is an exhaustiveness gap.
  MINOR NEED-CHANGE (name it + trace positions-only + state the relocation).

(b) WHY ADJACENT CAPS DISAGREE - addressed via the Voronoi-snap kill. The fan
  dump proves the disagreement mechanism: v4/v509 at x 0.0107 eps apart = two
  caps at two criticals, each extending V independently. The constructional
  analog ("both caps consume one shared V", the M4 one-arrangement lift) is the
  Voronoi cap-input snap (candidate A), validated-safe but NON-CLOSING
  (residuals 1+2) - the memo cites this correctly. But see the STRONGEST HOLE:
  the memo mis-describes WHICH separation keeps them apart.

(c) BUNDLE-COLLAPSE vs THINNING - partial. C2 thins CROSSINGS to eps-spaced
  representatives (the design keeps reps >= eps apart so no feature vanishes -
  which subsumes collapse-to-one for sub-eps crossing bundles). The 2D
  block-rule COLLAPSE-TO-ONE of a >eps VERTEX cluster is the fan root and is
  folded into the re-emission remainder (correct), but "collapse (vs thin)" is
  not named as its own candidate; the reasoning (>eps move loses a feature bound)
  is present implicitly. Minor.

(d) NEAR-X CRITICAL-PLANE MERGE - the candidate that falls out of the CORRECTED
  signature (see STRONGEST HOLE). Genuinely not in the killed list under its own
  signature; relocates to residual-2 (capEps=8*eps over-merges genuinely
  distinct near-degenerate verts) + residual-1. NEED-CHANGE to name+route.

## ATTACK 2 - kill conclusiveness

- kAngleTie (C1): population-corroborated. fan_light.log shows minTieGap/kTie=0
  and nearTieFans(<=100kTie) tiny (0-10) across all 4 cases; the MECHANISM (tie =
  angularly-coincident twin sheet) is directly observed in the dump (t429/t434 at
  ang=4.510891518 EXACTLY, gap=0). Solid, though "gap=0" is a min-statistic
  generalized to all tie fans.
- weld-bump (C3): SINGLE-SITE (Havocglass8 only, "does NOT resolve at ANY
  multiplier up to 4x"). Mechanism corroborated to population-of-2 by the
  identity-ext input-snap on GenericTwin7863 ("fixes the micro-edge but opens a
  macro hole elsewhere"). Adequate but thin; the "reopen ELSEWHERE" at the
  OUTPUT level is asserted, not separately dumped.
- C4 (yz-Canonicalize): MIS-TARGETED (STRONGEST HOLE) - killed on the case-(2)
  distinct pair, not the case-(1) fan twin.

## ATTACK 3 - is the coupling proven irreducible + costed?

The RELOCATION pattern is shown for C3 (->degenerate collapse) and C4
(->x-extent collapse) and C1 (structural). But "every bounded shortcut reduces
to the same problem" (line 1374) is asserted over a set of 3 killed + 1
validated that is not proven exhaustive: the provenance-identity route (a) and
the near-x critical merge (d) are missing/mis-routed. The real fix is costed by
SHAPE ("coordinated cap+strip re-emission; recognize the cluster as one entity;
one shared subdivision AND one shared junction vertex; precondition =
spread-vs-distinct disambiguation") but thinner operationally than the RSI-#3
exemplar (which named the DCEL + radial face-side ordering). Acceptable for a
research-remainder memo; the coupling is well-motivated, not merely asserted.

## ATTACK 4 - not-landing adjudication + "one wall"

C2 not-landing: SOUND and fairly presented. vrf_selfA/B.log: tol=off FATAL
budget; tol=0.5/1/2/4 all RESOLVE valid, volumes BITWISE-identical across tol
(self_A 0.142662155468, self_B 0.142664431424), A~=B to 2e-6 = real corroboration
but NOT an external a+b oracle. Both directions costed (landable IF owner accepts
un-oracle'd single-mesh MustResolve + hull budget->fan slowdown + core-stage
change). Consistent with the absolute zero-oracle-wrong posture and the
Voronoi-snap precedent.

"ONE wall, scale continuum": DEFENSIBLE but partly interpretive. The discriminator
(maxSpread/eps: controls <0.04, all failing >1.5, GT7081 126) is clean and
continuous. But equating the A2 crossing-bundle spread with the A1 cap-plane
endpoint twin as "the same near-concurrence at different stages" is an
interpretation from CO-OCCURRENCE + shared >eps signature, not a demonstrated
single mechanism; the notebook itself carves hull/self_A OUT as a distinct
size sub-class. Honest as "one wall + a size sub-class".

## ATTACK 5 - template conformance

Magnitudes-only / by-kind: yes. ASCII: yes. Supersession notes: present and
careful (STRICT-FP supersedes earlier eps-gate sections). Machinery delta: ZERO
confirmed - git diff src/ empty, no OV3_THIN/FANDBG/getenv in src, git status
clean. Overclaim check: the ">eps in x" anatomy line is the one factual
overreach (STRONGEST HOLE); the STRICT-FP "bitwise" was already self-corrected
to "only Cray bitwise". No other overclaims spotted.

## Bottom line
NEED-CHANGE. Fix the inverted twin-separation axis in anatomy part 1, re-scope
the C4 kill to the case-(2) distinct pair, and explicitly route the case-(1)
fan twin-merge family (yz-Canonicalize / near-x critical-plane merge /
provenance-identity label-weld) to its true kill (residual-1 cluster-collapse +
C3 any-axis weld-bump). The irreducible-coupling verdict itself stands.
