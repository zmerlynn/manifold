# reg3d OPENSCAD REOPEN - exact second-endpoint reconstruction (the terminal skipped this)

Lane: REOPEN, OPENSCAD. Branch explore/sweep-plane-3d-v5 (HEAD 3776a29b). Canonical
READ-ONLY; ALL work in rsync copy /tmp/oscad-reopen (fresh cmake -B vbuild Release -j8).
Deliverable = patch + scorecard (bake-off; orchestrator lands). Notebook here only,
no commits. ASCII; halfedge vocab; zero-oracle-wrong ABSOLUTE.

## THE MEASUREMENT THE TERMINAL SKIPPED
The terminal (reg3d-oscad) classified 259/323 truncated seams as "no fold-owned
completing endpoint" and declared a plane-based-rep PROOF wall. But its classifier
only searched for an INTERIOR-PIERCE dropped candidate. It never reconstructed the
seam's TRUE second endpoint. This lane does: for every one of the 323 nPts==1
truncations, reconstruct the EXACT tri-tri intersection segment in rationals
(doubles are dyadic -> Fraction is exact; hex-float probe dump = lossless).

## BUCKET TABLE (exact tri-tri intersection of all 323 truncated pairs)
  256  POINT_contact   exact intersection is a SINGLE POINT (not a segment)
        128  coincident-position vertex, distinct global id (unwelded CSG duplicate)
        128  vertex-on-edge / vertex-on-face T-junction (minVtxDist median 0.8, NOT
             near-coincident) - a vertex lands exactly on the other tri's edge/face
   52  EDGE_on_edge    segment lies on a shared BOUNDARY (collinear edge overlap);
        the two tris share a collinear edge that partially overlaps, fold apart -
        measure-zero contact, no interior crossing (31 share a coincident vtx, 21 not)
   15  XSEG_transversal  GENUINE crossing: segment interior passes through BOTH tri
        strict interiors.  12 share a coincident vtx whose MISSING endpoint is that
        vertex (position-wjump shape); 3 are sub-eps SHORT (seglen<=eps).

WHERE the true missing endpoint lives + WHY production fails to produce it:
- 308 pairs (256 POINT + 52 EDGE_on_edge): there is NO genuine transversal crossing;
  the exact intersection is a point or a boundary segment.  The filter/SoS recorded a
  PHANTOM single endpoint (a vertex-on-face touch decided a pierce, or a collinear-
  edge endpoint), nPts==1.  The index-keyed sharesVert skip MISSES the coincident-
  duplicate family (the terminal measured "nSh=0" by INDEX and wrongly concluded "not
  a shares-vertex instance") - it IS the third incidence family, keyed on POSITION.
- 15 pairs: genuine crossing; the missing endpoint is (12) the coincident vertex the
  index-wjump misses, or (3) a sub-eps interior pierce below the weld radius.

## SOUNDNESS (exact): of the 308 measure-zero pairs, ZERO have a clean off-plane edge
piercing the other's STRICT interior.  Skipping them drops NO real crossing
(verify_skip.py: bad=0).  Only the 15 XSEG are genuine.

## FIX (decision-completions, wjump/7863c1 precedent, no new predicate FORM)
1. POSITION-coincidence generalization of the shares-vertex skip/recovery: a pair
   meeting only at a coincident-position vertex (distinct id) is self-adjacent-at-a-
   corner; apply the SAME wjump gate (skip unless a genuine off-vertex pierce, in
   which case recover [coincidentV, offVertexP]).  Handles the 171 shared cases (128
   skip + 31 skip + 12 recover).
2. PHANTOM GUARD at nPts!=2: a lone/absent endpoint is a genuine seam only where a
   clean off-plane edge pierces the other's strict interior (self-certifying
   transversal).  Else the pair is a measure-zero contact -> NO seam (skip), not
   fail-closed.  Handles the 128 T-junctions + 21 collinear-edge (nosh).  Leaves the
   3 sub-eps XSEG as honest residue (witness theorem: sub-eps seam = collapse).


## RESULT (LANDED in /tmp/oscad-reopen; canonical byte-untouched)
The F11 seam-truncation wall the terminal called a plane-based-rep PROOF is CLOSED.
323 truncations = 308 measure-zero phantom contacts (skip, exact-proven sound) + 12
position-wjump crossings (recover) + 3 sub-eps collapses. RecordSeams completes;
openscad now fails one wall DEEPER at F4 (GT7081-class winding-probe filter-precision,
kernel-tripwire-gated). Still fail-closed, zero wrong resolve.

VALIDATION:
- Openscad pin flipped F11 "degenerate incidence" -> F4 "seam sub-face arrangement
  not exactly resolvable". FULL Overlap3 suite 32/32.
- BITWISE (FNV) vs COMMITTED-canonical baseline: IDENTICAL across all 52 resolving
  carrier outputs. Fix touches openscad only.
- MUTATIONS (each reverts openscad to F11): MUT-A index-coincidence 12 FAIL; MUT-B
  no-guard 152 FAIL; MUT-C no-collapse 3 FAIL. All 3 mechanisms load-bearing.

FIX (3 decision-completions in RecordSeams, no new predicate FORM):
1. shares-vertex skip/recovery keyed on POSITION coincidence (unwelded duplicates
   = the shares-vertex family the index-keyed skip missed; terminal's nSh=0 was by
   INDEX). 2. phantom-seam guard: no clean interior pierce -> measure-zero -> skip.
   3. sub-eps collapse: recorded endpoints merged under eps (nPtsPre>=2 -> nPts==1)
   -> witness-theorem collapse -> skip.

CONTAMINATION NOTE: canonical /home/zml working tree was concurrently modified by
another lane (a BroadphasePairs refactor, uncommitted). This lane rebased on the
COMMITTED HEAD 3776a29b (git show HEAD) - immune to it. Orchestrator reconciles.

PATHS: patch /tmp/oscad-reopen/oscad-reopen.patch ; scorecard /tmp/oscad-reopen/scorecard.md
