# 8. Collide fast pairs again mid-step

Date: 2026-10-27

Status: Accepted

## Context

Contact points are computed once per step. A body spinning fast moves its real overlap
away from those points within the step, so the solver pushes on stale points. Two
cylinders spinning at 25 rad/s sank 48.7 mm into each other.

Options tried:

- More substeps. Costs every pair, not just the fast ones.
- Stronger pushout for deep overlap. Clears the worst cases faster, but some piles of hulls
  never go to sleep with it.
- Collide again every substep. Correct, but the narrowphase cost multiplies by the
  substep count.
- Collide only the fast pairs again, once, halfway through.

## Decision

Pairs whose relative motion over the step exceeds the speculative distance are collided
again at substep n / 2. Before that, every active contact is rebased: anchors become the
bodies' current lever arms, the stored separation absorbs the motion so far, and each
body's motion since the base pose resets. New points take impulses from old points within
the speculative distance, and the restitution data of the nearest old point.

## Consequences

Measured over 36 seeds, frames with overlap past 50 / 100 mm:

| Scene | Before | After |
| --- | --- | --- |
| Rain | 907 / 110 | 589 / 53 |
| Hulls | 298 / 20 | 129 / 0 |

- The spinning cylinder case drops from 48.7 mm to 5.9 mm.
- Costs 13 to 20% in scenes with many fast bodies, nothing in scenes without.
- Rebasing is required. Without it the new anchors act like rolling resistance.
- Carrying the restitution data is required. Without it a full-restitution bounce loses
  about 15% of its height.
- Spinning bodies trapped deep under a pile still take 8 to 12 frames to recover. That
  case exists without the refresh and is still open.
