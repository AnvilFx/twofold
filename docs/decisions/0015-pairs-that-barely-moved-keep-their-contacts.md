# 15. Pairs that barely moved keep their contacts

Date: 2026-10-27

Status: Accepted

## Context

Every touching pair was collided from scratch every step. In a resting pile almost
nothing changes from one step to the next, and the narrowphase redoes the same work.

Jolt keeps a pair's contacts while the pair moves less than 1 mm and turns less than
about 2 degrees relative to itself, and only measures the separations again.

## Decision

A contact stores the pose of body B in body A's frame, the normal in A's frame, and each
point's position on both bodies. While the pair stays within `contactReuseDistance` and
`contactReuseAngle` of that pose, the points are kept. Each separation is measured again
from the two stored positions, and the normal turns with body A. Otherwise the pair is
collided again.

Changing a shape's geometry clears the stored pose of its contacts, and so does the
mid-step refresh.

## Consequences

- About 10 to 15% faster on busy scenes.
- At 1 mm, Dominoes spent 926 frames past 5 mm of overlap instead of 500. The default is
  0.5 mm and 0.035 rad, which leave it as it was. Rain and Hulls overlap the same or less.
- A kept contact doesn't notice a new feature coming into reach within those limits. The
  limits are small enough that the speculative margin covers it.
