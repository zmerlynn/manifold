# Post-flip cleanup pass - landing record

Owner-queued heavy cleanup on the single coordinated-emission engine tree
(after the act-3 flip deleted the per-face world). Codex gpt-5.6-sol executed
in a disposable copy; coordinator (this session) reviewed + re-verified on
canonical.

## Removed (all verified dead at current HEAD)
- Comment rot: references to deleted per-face functions (EmitSeamedFace,
  FoldCoplanarClusters, per-face RemoveOverlaps2D) rewritten to the single-engine
  reality.
- Dead pipeline superseded by the coordinated registry engine: EnumerateTriplePoints,
  EnumerateWedgeSplits, BuildJunctionRegistry (~621 lines of pass bodies + fields).
- Dead functions: ExtractCells (defined, 0 call sites - confirmed by grep: only a
  definition + one comment mention), WindingAtCands, unused sheet-normal plumbing,
  diagnostic-only counters.
- Instrumentation prune: getenv sites 34 -> 2, fprintf sites 73 -> 1; 25 obsolete
  dev-lever names + all file/timing/targeted dumps removed.

## Kept (load-bearing - NOT removed)
- IX_OFF: default-ON behavior kill-switch.
- E1_FLOODDIFF: the flood-vs-exact shadow validator (smoke-tested).
- The mutation levers the pins depend on.
- Six unused helpers inside the frozen sos:: kernel: FLAGGED for owner approval,
  left untouched (kernel is frozen).

## Verification (coordinator, on canonical, from scratch, PAR=ON)
- Kernel sos:: / Orient3DExactSign core byte-identical (md5-compared).
- ExtractCells etc. confirmed 0 call sites before accepting the removal.
- Build clean; fast corpus 29/29; openscad pin PASS; GT7081 individual OK;
  GT7081 JOINED (Corpus_GenericTwin7081_Resolves) OK 254s - the decisive
  engine+partition exercise Codex had skipped, run here.
- Net: src/overlap3.cpp -~1674 lines.

## Follow-up
- Owner decision: remove the six unused sos:: helpers? (needs approval - frozen region.)
- Owner decision: squash the superseded retry (dfbd6950/592aa834) before any PR?
