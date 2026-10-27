# 14. Edge queries pair edges near the overlap

Date: 2026-10-27

Status: Accepted

## Context

The SAT edge query tries every edge of A against every edge of B. For two 96 edge shapes
that is over 9000 pairs, and it was most of the narrowphase for hulls and cylinders.

When two shapes overlap, the edge pair that decides the contact meets within the overlap
depth of both shapes. An edge axis only beats the best face axis with less overlap than
that face, so the depth is at most the face's.

When two shapes are apart, this doesn't hold. The edge pair that proves they are apart
can sit far outside the overlap of their bounds. Pairing only nearby edges then misses it,
and the pair comes back touching, sometimes deeply.

Considered:

- Pairing every edge, as before. Correct and slow.
- GJK distance for every pair before the SAT. Adds 200 to 900 ns to every pair, including
  boxes, where the full edge query is cheap.
- Pairing only nearby edges, and confirming with a GJK distance check when any edge was
  skipped.

## Decision

Edges are paired only when they reach into the overlap of the two bounds, grown by the
larger of the margin and the best face's overlap. When this skips any edge, a GJK distance
check runs, and the pair has no contact if the shapes are further apart than the margin.

## Consequences

- Rocks against rocks are about 4 times faster, cylinders against cylinders about 3.
- Over 400,000 random pairs, contacts agree with pairing every edge. The only differences
  are pairs at the margin that are actually further apart, which no longer get a contact.
- Skipping the pruning for small shapes, so boxes avoid the GJK check, was slower
  overall. Every pair is pruned.
