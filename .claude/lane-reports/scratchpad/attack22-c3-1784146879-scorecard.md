# C3 lane scorecard (attack22-c3-1784146879)

QUESTION: is the branch-node radial order decidable at level 0 from input-plane
provenance?  (the last of the 22, Regularize3D)

## HEADLINE

- LEVEL-0-ORDER: 13/13 branch nodes (10/10 multi-halfedge fans) orderable at level 0
  from input-plane normals + the emitted fwd bit, with ZERO genuine angular ties and
  ZERO exact sheet overlaps (rigorous exact comparator).  The radial order is NOT the
  tripwire.  (Under the task's literal "same-plane = tie" criterion the count reads
  0/13, but those same-plane pairs are benign antipodal through-planes labeled by fwd.)
- TRAJECTORY: 22 -> 22.  No winding-only completion closes any open edge.
- ORACLE VERDICT: residue is a genuine near-tangent REPRESENTABILITY wall (missing
  sub-face EXTENT), NOT an ordering wall.  Nothing shipped (a completion here distorts).

## STEP-BY-STEP

MEASURE (step 1):
  - reproduced baseline EXACTLY: 22 open edges, 43 dangling HEs, 20 verts, deg hist
    {2:7,4:7,6:4,8:1,12:1}, balanced 1-cycle, 13 branch nodes.
  - per-fan (exact rationals): 9/10 fans have two HEs from the SAME input face; 1/10
    all-distinct.
  - rigorous steelman: those 9 are near-ANTIPODAL through-planes; oriented input-normal
    order == true constructed order on 10/10 fans, 0 genuine ties.  fwd flag consistently
    labels the antipodal pair (oppFwd && oppOracleSign on every same-plane pair).
  => order level-0 decidable; the pairing is NOT the tripwire.

CLOSE WHAT CLOSES (step 2):
  - every open edge has a fwd/bwd deficit (sum 23) => sheets genuinely MISSING, not
    mis-ordered.
  - 5 components: comp0 (9 edges) lies in NO single input plane (dev 8.5e-7 >> eps) ->
    multi-plane synthesis, gated; comp1-4 (13 edges) single-plane candidates.
  - production-winding probe of comp1-4: every component has faces the both-sides
    {w>=1} rule REJECTS.  comp2's loop-triangle validates at its centroid (coarse
    22->19) but a finer read drops a third -> region not uniformly boundary
    (reproduces f4-r5 "finer read is worse" from the output side).
  => trajectory 22 -> 22; nothing shippable without distortion.

RAILS (step 3):
  - no env: full Overlap3 suite 32/32 PASSED (byte-clean, openscad fail-closed).
  - C3_COMPLETE=1: full Overlap3 suite 32/32 PASSED (probe never mutates emitted soup;
    fires only on the fail-closed path; resolving carriers unchanged by construction).
  - mutation: instrumentation is load-bearing for NOTHING (off->22, on->22); that IS
    the terminal.
  - kernel: zero new predicate FORM; probe reuses RobustWinding (blessed caller).
  - zero-oracle-wrong: no geometry emitted; rejected candidates prove the wall.

TERMINAL / PRICE (step 4):
  - the terminal is the exact per-face near-tangent 2D arrangement (missing sub-face
    EXTENT) = r6's ~600-1300 LOC; the radial rule is FREE (level-0), removed from the
    cost.
  - the tripwire predicate (orient2d on constructed near-tangent crossing): design-b's
    plane-based symbolic-rational {f,g,h} on the EXISTING adaptive-width accumulator
    FITS - ONE new FORM, ZERO vendored code, ~150-250 LOC (design-b B2); vs vendored
    Shewchuk predicates.c ~4-5k LOC.  design-b fits.

## FILES

- notebook: .claude/lane-reports/attack22-c3-1784146879.md
- patch (measurement, env-gated, NOT for landing): attack22-c3-1784146879.patch
- offline: /tmp/attack22-c3/{analyze_c3,steel2,pernode,coplanar_test,loops}.py
