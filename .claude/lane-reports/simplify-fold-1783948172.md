# Simplify A+B distillation - FOLD lane (verify -> doc)

Branch explore/sweep-plane-3d-v5. Folds the three NEED-CHANGE findings from the
adversarial verify lane (.claude/lane-reports/simplify-verify-1783947535.md) into
docs/Regularize3D.md. Verdict there was NEED-CHANGE, not BREAK: the design
survives, the fail-closed posture holds; these are doc-level distillation gaps
(an overclaim, a hidden cost, an inverted mechanism). Docs-only diff, ASCII.

## Folds landed

1. R1 (internal weld) - REMOVED the overclaim "the re-gate catches the failure
   but does not repair it." The re-gate is structurally BLIND to the weld-merge
   self-fold it was cited to backstop: when the assembly weld merges two
   exact-distinct B verts, the straddling triangles now share that vertex
   position, so IsSelfIntersecting's shares-vertex skip drops the pair - the fold
   is invisible. Only the validity half (IsManifold/Is2Manifold) catches a weld
   outcome, and only the non-manifold one (pinch/tear), never a fold. Also stated
   that a constant-radius bounded weld is NOT the escape (wall-A weld-bump kill:
   dead at every multiplier - too small re-manufactures the twin, too big
   collapses slivers to holes). R1's honest rebuttal is now B's STRUCTURAL
   defense, not the gate: once-only construction + radial (not projection,
   not self-location) assembly. Named the surviving residual (two genuinely
   distinct arrangement points rounding within eps) as open.

2. LOCALIZER dissolution - ADDED the hidden cost consequence to B's WINDING
   bullet. Each seed cast is an O(ntri) winding query scoped to the component, so
   ONE large component pays the whole-component cost - seconds-per-query scale on
   the largest corpus component (GT7081), no adaptive orient3d kernel in the tree
   today = net-new surface area. Cited v5-verify-probe-1783913337.md.

3. R2 (gate-notion mismatch) - CORRECTED the inverted parenthetical: the 2*eps
   shares-vertex relaxation SUPPRESSES near-miss detection (returns
   non-intersecting when an eps nudge separates), it does not flag it. Stated the
   gate is systematically CLEAN-biased. Kept gate-clean/B-dirty (a near-degenerate
   band just outside the eps weld) as THE silent-miss direction; noted the reverse
   (gate-dirty, B-clean) is near-empty (B-skipped subset of gate-skipped).
   Cross-referenced that R2(i) and R1 are the SAME blind spot - the weld-merge
   fold lands in exactly the shared-vertex config the gate skips.

## Fence

Docs-only: `git diff --stat` = docs/Regularize3D.md only (no src/, no test/).
Non-ASCII scan of the edited doc: clean (pure ASCII). Verify lane's ATTACK 3/4
(fidelity, style) were SURVIVE and are untouched here.
