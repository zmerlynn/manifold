# NEXT SESSION BRIEF: M4-literal, main-agent-authored

THE TASK (user-approved): implement E' "one arrangement per
critical, three consumers" LITERALLY, as main-agent surgical work
(not a subagent lane - three lanes have escaped this mechanism).
Then burn down the judge's 13, then re-run both closing audits.

STATE: branch explore/sweep-plane-3d-v3 at 3d0046ec (pushed), all
green: 33+1 Overlap3 gates, 111 2D fences, 579 full suite. The
spec is docs/SweepEmit3D.md, section E' (the [R2-fold] acyclic
dataflow paragraph is the contract). The judge's 13 items:
.claude/lane-reports/v3r3-simplicity.txt; the fidelity findings:
v3r3-fidelity.txt.

WHAT M4-LITERAL REQUIRES (and what exists):
- Today: ComputeCap runs TWO RemoveOverlaps2D-class passes over
  shared rawVerts (same input order -> same snaps, empirically
  consistent) and strips SELF-EXTEND their c-side boundaries; the
  cap and strip subdivisions coincide only by shared-snap luck.
- Spec: ONE arrangement pass per critical over both extended
  limits (L @ +1, R @ -1); its output subdivision serves (a)
  cap_plus loops (region diff > 0), (b) cap_minus loops, and (c)
  REPLACES both adjacent slabs' strip-edge c-side intervals - the
  strips must consume the subdivision, which likely means the
  strip-emission pass runs AFTER the cap pass per critical, or
  strips are emitted with provisional c-side boundaries that the
  cap pass then refines before BuildImpl.
- WHY LANES BAILED (their notebooks, same dir): the strip pass and
  cap pass are structured as independent loops over slabs vs
  criticals; the three-consumer rule inverts that ordering. Read
  v3-round3-*.md and v3-final-*.md notebooks for the two concrete
  bail points before writing code.

RITUAL AFTER THE CODE: scripts/format.sh; full suites; commit;
push; re-run BOTH audits (prompts pattern in
/tmp/codex-branch-review/prompts3/, or reconstruct: the simplicity
judge prompt sed'd to docs/SweepEmit3D.md + the new commit; a
fidelity prompt verifying M4 + the 13). The judge reaching SURVIVE
closes the arc; NEED-CHANGE with a shorter list gets recorded
honestly in docs/SweepEmit3D.md's close section, then the user
reviews.

HOUSE: lab notebook for MY OWN work too
(.claude/lane-reports/m4-main-agent-<ts>.md); ASCII; la::;
DEBUG_ASSERT; eps metric table only; fences frozen; the fence rule
(no test weakening, ever) applies to me most of all.
