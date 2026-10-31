# 35. Box pairs test their axes in closed form

Date: 2026-10-31

Status: Accepted.

## Context

Two boxes went through the general polyhedron SAT: two face queries over vertices, an edge
query over 12 by 12 edges pruned to the overlap of the bounds, and a GJK check whenever that
pruning skipped an edge. Box pairs took 17% of a step in Sweeper, 14% in Dominoes and 6% in
Rain, at 0.2 to 0.8 µs a pair.

A box has three face axes and three edge directions, so the 15 axes of a pair and the
separation along each follow from the half extents and the rotation between the boxes.

Options:

- Test the 15 axes in closed form and clip the contact as before.
- Also clip box faces with a routine of their own. Clipping took about 0.2 µs of a touching
  pair, but the clipper is shared with every polyhedron, so it is better made faster for all.

## Decision

A pair of boxes tests all 15 axes in closed form, in A's frame. An edge pair counts only
where its arcs cross on the Gauss map, as in the general query. The face and edge numbers
come from a table built once from makeBox, so the same face and edge clipping makes the
contact. The edge pair kept within the step is not used.

When the boxes are apart, the gap is checked against the features meeting across it: for a
face, the deepest vertex of the other box; for an edge pair, the closest points of the two
edge lines. If those lie on both boxes, the gap is the distance between them and GJK is
skipped. Otherwise GJK confirms they are within the margin.

## Consequences

- The box queries cost about 0.1 µs a pair, from 0.2. GJK ran on 69k of 260k gapped pairs in
  Sweeper and on 9k of 175k in Rain. Box pairs took 260 ms of Sweeper's 12 seeds of 600 frames,
  from 410, and 170 ms of Rain's, from 375.
- Whole steps on one core, 12 seeds of 600 frames, best of four, before and after:
  - With sleep off: Dominoes 0.109 to 0.102 ms, Sweeper 0.285 to 0.276 ms. Other runs put
    Sweeper 6% to 11% faster. Rain, Meshes, Joints and Pyramid are within noise.
  - With sleep on: Dominoes 0.088 to 0.080 ms, Sweeper 0.256 to 0.239 ms.
- On overlapping pairs the closed form agrees with the general query: in Rain, Sweeper,
  Joints, Meshes, Dominoes and Pyramid, no normal was off by more than 1 mrad, no depth by
  more than 0.1 mm, and no pair switched between face and edge. Pairs with a gap can differ.
  The 0.9 rule lets an edge axis win when it separates a little less than the best face.
  The general query had often pruned that edge pair and took a worse one, while the closed
  form takes the edge pair that separates most.
- Results change wherever boxes touch boxes. Over 12 sets of 32 seeds, frames with overlap
  past 10 / 20 / 50 mm, before and after:

  | Scene | Before | After |
  | --- | --- | --- |
  | Rain | 89235 / 18804 / 1688 | 86408 / 21469 / 1747 |
  | Sweeper | 379 / 22 / 0 | 347 / 18 / 0 |
  | Meshes | 27744 / 3383 / 80 | 27446 / 2971 / 79 |
  | Joints | 309760 / 27666 / 2314 | 309905 / 27649 / 2317 |
  | Dominoes | 3 / 0 / 0 | 2 / 0 / 0 |

  Rain's sets range from 1209 to 2900 frames past 20 mm under both, and the count before
  decision 34 was 22414. Pyramid and Tower have no frame past 10 mm either way. Hulls,
  Cutouts, Flat layers, Tower and Slopes keep their hashes.
- One gapped Pyramid pair out of 1276 that the check kept was one GJK put past the margin,
  by 38 µm.
