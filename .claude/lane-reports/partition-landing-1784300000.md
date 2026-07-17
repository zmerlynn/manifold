# Reliability-partitioned winding field - landing record

Closes the GT7081-joined STOP with the FINISHED local-criterion form,
replacing the rejected field-retry (dfbd6950) per owner law (fallbacks are
fundamental smells).

## Design (as landed)
- RELIABILITY = DECIDEDNESS. Each flood edge carries `decided`: true iff every
  predicate contributing to its winding delta was decisive (exact result or
  filter-certified sign); false for near-tangent incomplete chord transport,
  sub-ULP band endpoints, non-conserved parity, or uncertified nonzero handoff
  corrections. Zero-delta clean handoffs are decided by equal exact offsets.
- BFS floods over decided (`flOk`) edges only -> the cells partition into
  reliable-connected subgraphs.
- The hull-vertex seed anchors the outer subgraph; every other reliable subgraph
  receives exactly ONE absolute anchor = one exact winding probe at a
  well-conditioned interior cell. A subgraph with no probeable cell fails closed
  (honest, not a fallback).
- DELETED: the whole-field nuke and the `keepField` retry. The driver calls
  emission once. A conflict among decided edges is an invariant violation ->
  debug-assert + fail-closed (never a regime switch).
- The ladder non-boundary acceptance heuristic is gated to sub-weld cells
  (extW <= 0.99*eps), so it can never override an exact combinatorial jump
  outside that regime (the 872-drop bug).

## Execution / verification
- Implemented by Codex gpt-5.6-sol (xhigh) in a disposable copy /tmp/partition-codex.
- Coordinator review (this session, on canonical): diff localized (~+110 net
  logic lines in overlap3.cpp; test change is comment-only, assertions
  byte-identical - test NOT weakened); no kernel (sos::/Orient3DExactSign)
  changes; no disguised retry/fallback; grep-clean of keepField/field-nuke.
- From-scratch Release (PAR=ON) canonical re-verify: fast corpus 29/29;
  openscad pin PASS (12.3s); GT7081 individual OK; GT7081-joined
  Corpus_GenericTwin7081_Resolves OK (251.6s, RED->GREEN), volume oracle passed.

## Load-bearing anchor
The joined test was RED at ef8c7cbf and GREEN after this change alone, with
its assertions unchanged - the partition is the load-bearing fix.

## Pending
- Stage-2 adversarial verify of the WHOLE act-3 arc (flip + flood + perf +
  this closure) before push.
- Owner-queued heavy cleanup crucible on the single-implementation tree.
- Optional: squash the superseded retry (dfbd6950/592aa834) before the final push.
