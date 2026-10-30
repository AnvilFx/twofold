# 34. Fast pairs keep their edge pair within the step

Date: 2026-10-30

Status: Accepted. Extends [decision 24](0024-fast-solid-pairs-collide-every-substep.md).

## Context

Colliding fast pairs again within the step took 26% to 34% of a step in Rain, Hulls,
Meshes and Sweeper, nearly all of it in the narrowphase. Pairs of polyhedra were most of
that in Hulls and Sweeper, and in Hulls their edge query, its box test and the GJK check
after it cost about twice the two face queries.

A pair collided again had usually not changed its mind. In 70% to 93% of those collisions
the new normal was within 0.8 degrees of the old one.

Options:

- Clip along the feature found last, faces or edges, with no search. Hulls ran 20% to 25%
  faster, but 8% of the reused features weren't the ones a search would pick, and 0.8% put
  the deepest point more than 5 mm off.
- Search the faces again and test only the edge pair found last. Every switch to a face is
  still caught.
- Either one only while the pair has turned little since the step started.

## Decision

A polyhedron manifold records the faces or edges it was clipped from. Within the step, a
fast pair that turned less than 0.1 rad against itself since the step started runs both
face queries as before, but tests only the edge pair it last met on, while that pair's
arcs still cross on the Gauss map. If they don't, the full edge query runs. If it last met
on a face, no edge pair is tested. Collisions at the start of a step search everything.

## Consequences

- Whole steps on one core, 12 seeds of 600 frames, best of three, before and after:
  - With sleep off: Hulls 0.43 to 0.38 ms, Sweeper 0.32 to 0.30 ms. Rain, Meshes and Joints
    are within noise.
  - With sleep on: Hulls 0.32 to 0.27 ms, Sweeper 0.27 to 0.25 ms.
- Results change wherever fast polyhedra touch. Over 12 sets of 32 seeds, frames with
  overlap past 10 / 20 / 50 mm, before and after:

  | Scene | Before | After |
  | --- | --- | --- |
  | Rain | 99178 / 22414 / 1693 | 89235 / 18804 / 1688 |
  | Hulls | 4189 / 1079 / 65 | 5133 / 858 / 46 |
  | Sweeper | 362 / 19 / 0 | 379 / 22 / 0 |
  | Meshes | 26892 / 2556 / 78 | 27744 / 3383 / 80 |
  | Joints | 309742 / 27760 / 2284 | 309760 / 27666 / 2314 |

  Single sets swing more than these totals, and two sets hold the Meshes change.
  Pyramid, Tower, Slopes, Dominoes, Cutouts and Flat layers keep their hashes.
- Rejected:
  - Reusing face features too. Hulls frames past 20 mm went from 1079 to 1269, and one
    Joints frame reached 623 mm.
  - Reusing the edge pair with no limit on the turn. Hulls frames past 20 / 50 mm went to
    1643 / 110.
- A pair whose deciding axis moves to a new edge pair within the step keeps the old one
  until the next step. Under 1% of reused contacts came out more than 1 mm off what a search
  gives, and all of them were this case.
