# VERIFY LANE: push-gate for the input-exact stack + fold member-guard fix

Branch explore/sweep-plane-3d-v5. Range 7f8aec4e (pushed) .. 2d7d42c6 (local HEAD):
- 1568e9e0 = input-exact/coverage stack (all levers default-inert + P3 cap provenance)
- 2d7d42c6 = fold member-guard fix (+12 LOC ACTIVE default-path change, openscad 21->18)

From-scratch Release builds, both SHAs, worktrees /tmp/plv (HEAD) + /tmp/plv-base
(base), Ninja, -DMANIFOLD_PAR OFF+ON. FNV/determinism measured with an env-gated,
print-only hook at the two resolver returns (RemoveOverlaps3D + ResolveComponentDirect),
byte-identical in both trees, no-op unless ZZZ_FNV set (hashes vertPos_ + halfedge Start
ids + sizes = the BitIdenticalMesh basis). Canonical untouched (worktrees are throwaway).

## VERDICT: PASS (push-safe)

Diff scope is exactly src/overlap3.cpp + the two notebooks; no test/header/other-src
change (git diff --stat). Every claim reproduced.

### 1. Suite / heavies (HEAD, both PAR)
- Overlap3 32/32 PASS at HEAD PAR=OFF and PAR=ON (base also 32/32).
- GT7081 resolves (nopar+par, and under ulimit -v 4000000 timeout 900).
- openscad fail-closed: HEAD large-comp openEdges=18, base=21 (both confirmed);
  fail reason unchanged ("unresolvable sheet contact", NonManifoldEmission) - the
  openscad test asserts the reason and passes at both SHAs.

### 2. FNV base vs HEAD (full Overlap3 suite, per-test attributed, 26 resolving carriers)
- Byte-identical on every resolving carrier EXCEPT EntangledBarsRotated.
- EntangledBarsRotated: base h=a53e916c745b6337, HEAD h=0e95fb8f0a5bdb06, BOTH
  nv=32 nh=180 (same vertex set, same 60 tris). Test passes at HEAD asserting its
  own oracle: volume == INDEPENDENT (A+B) boolean union, GWN over 20k pts (checked
  >3000), one Decompose solid, tol-invariant, not self-intersecting. Nothing beyond
  the re-triangulation changed -> rebless valid.

### 3. Fold fix (2d7d42c6)
- Diff read: the onCount>=2 skip is now gated behind `if (kMemberGuard)`
  (FOLD_MEMBER_GUARD lever); default path splits at EVERY interior registry junction
  unconditionally (once-only). Rationale comment matches the measured bug (cap
  sub-face spanning the shared line unsplit; 21-18 = the 3 cap-wall T-junction opens).
- Mutation FOLD_MEMBER_GUARD=1 restores openscad opens 18->21.

### 4. Stack (1568e9e0) default-INERT
- ATTRIBUTION: HEAD+FOLD_MEMBER_GUARD reproduces the base EntangledBarsRotated hash
  EXACTLY (a53e916c745b6337) and base openscad opens (21), with the full stack
  present -> the stack contributes ZERO bytes; the fold guard is the sole change.
- Levers EX2_*/GRAZE_*/STITCH_*/SYM_*/IX_OFF/EX2_CLUSTER/EX2_ALL/EX2_OFF all env-
  gated (default off); BuildSymbolicProvenance early-returns when no wedge face is
  flagged; F4B_DUMP inert when unset. The empirical FNV identity is the proof
  (F4B_TRIPLES=0 equivalent).
- sos:: kernel: the diff of the sos namespace is PURE ADDITIONS (0 base-only lines);
  Orient3DExactSign / EdgePiercesTriSoS / SumSignN<8> / TermN<NMag> byte-identical.
- No new predicate FORM in the core kernel: P3 SegPlane{HPoint,EHPoint,BigHPoint} is
  a depth-1 point-pair-meet-plane on the existing HPoint/Poly form (X,Y,Z deg ~3, W
  deg 2, << deg-20), full filter/exact/big triad, dispatched consistently on
  VProv.segPlane (trivial/segPlane/cramer) in all three bases. The re-landed input-
  exact overlay adds a separate default-inert exact-arith BigOrient2D (the "landed
  form" the P3 rides) - NOT a deepening of the untouched core kernel.
- fail-closed / Fatal() arms: TEXT-identical base vs HEAD. Contract-eps snap untouched.
- inventories synced (3-basis triad complete); entire committed diff ASCII; no
  FNV/memcmp diagnostic in committed code.

### 5. Determinism (4t=PAR build, 1t=nopar build)
- siA: par==nopar==(par run B) h=81f1da4e0d8ff9ec (nv=8474 nh=50844).
- GT7081: par==nopar==(par run B) h=e54127a46f6b652a (nv=15754 nh=94380).
- openscad census: large-comp openEdges=18 on par runs A and B (== nopar).

Notebook only; DID NOT push.
