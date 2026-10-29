# 30. Bounds tree leaves are bodies

Date: 2026-10-29

Status: Accepted. Replaces [decision 16](0016-bounds-tree-with-persistent-pairs.md).

## Context

Every shape had a leaf in the bounds tree, and the world kept the pairs of shapes whose
leaves overlapped. A traced cutout is a body of several convex pieces, so a pile of them
kept a pair for every two pieces whose leaves overlapped, and moved a leaf for every piece.
In Cutouts the broadphase took 0.13 ms of a 0.46 ms step, more than the solver.

Options:

- Keep shape leaves, and make pair finding cheaper. The pair count stays.
- One leaf per body, and pair the shapes of two bodies only when their bodies' leaves
  overlap.

## Decision

Each body with a shape other than a half-space has one leaf, covering the bounds of those
shapes. The world keeps the sorted pairs of bodies whose leaves overlap, found and dropped
as decision 16 did for shapes. For a pair where a body moves, each body lists its shapes
near the other body's bounds grown by both bodies' travel, and the shapes of the two lists
whose grown bounds overlap and whose filters agree become candidates. Half-spaces stay in
their own list, and are paired with a body's shapes only when the body's grown bounds
reach the plane.

A body with more than 64 shapes keeps them in a tree of its own, so listing the shapes near
a box, a ray or a sweep doesn't test them all. Smaller bodies list theirs by testing each
one. Queries go through the body tree and then each body's shapes the same way.

## Consequences

- In Cutouts the kept pairs fall from 751 to 330, and with cheaper bounds for boxes and
  prisms the broadphase from 0.13 to 0.016 ms a step. With 9,800 static boxes as one body
  and 40 boxes falling on it, the step costs what it did with shape leaves. As 9,800
  separate bodies it is the same within noise.
- The candidates for the narrowphase are exactly the old ones, so results are bit for bit
  the same. Only the order the tree finds them in changed, and pairs are still sorted.
- A body's leaf moves when any of its shapes leaves it. A body with a few shapes far
  apart has a leaf over the empty space between them, and its pairs list shapes that
  don't touch.
- Adding a shape to a body with its own tree inserts one leaf. Changing a smaller body's
  shapes makes its leaf again.
