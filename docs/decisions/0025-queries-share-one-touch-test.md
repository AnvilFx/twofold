# 25. Queries share one touch test

Date: 2026-10-28

Status: Accepted

## Context

Overlap queries and shape casts have to agree with each other. When a cast reports a hit,
an overlap query with the geometry at the hit pose should find that shape, and a cast that
misses should leave a path where no overlap query finds anything.

The step already has the pieces, but none of them answer that exactly:

- The narrowphase reports contacts out to the speculative distance, and on cylinders its
  separations come from the feature axis estimate.
- `timeOfImpact` stops once the gap along the GJK normal is under its target. On a
  cylinder cap the normal tilts, and casts stopped up to 4.4 mm before the shapes touched.
  Changing it would change stepping.
- GJK on a cylinder can run out of iterations and stop above the true distance, by up to
  about 0.15 mm, and by a different amount depending on which shape comes first.

## Decision

- Two shapes touch when their GJK distance is within 0.01 mm, scaled by `lengthUnit`. The
  query geometry is always the first argument.
- When either shape is a cylinder and the distance is just above that, GJK runs again in
  the other order and the smaller distance is kept. Both orders are upper bounds.
- `overlapGeometry` uses that test against each candidate, and against each triangle of a
  mesh under the query's bounds. A mesh has no inside.
- `overlapPoint` tests each shape kind analytically to the same 0.01 mm. `overlapAabb`
  compares bounds only.
- `castGeometry` moves the geometry without turning. Each candidate is swept with its own
  conservative advancement on the same test: the step is the larger of the gap along the
  GJK normal over the closing speed and the GJK distance over the speed, and a step that
  lands inside is bisected back to the touch. `timeOfImpact` is left alone.
- A cast that starts touching reports fraction 0 with `initialOverlap` set. Its point and
  normal are the deepest narrowphase contact at the start, or a zero normal when the
  narrowphase finds none.
- A hit's point and normal come from the narrowphase run 0.5 mm before the hit, since a
  normal across a gap of a few micrometers is mostly rounding. The point is moved from the
  middle of the gap onto the hit shape.
- Candidates are tried in shape index order and the lowest index wins a tie.
  `castRayAll` sorts by fraction, then shape index.

## Consequences

- Over 600 random casts through scenes with every shape kind, the hit shape is found by an
  overlap query at the hit fraction, nothing is found 1 mm before it, and a missed path is
  clear at every sampled pose.
- Casts against cylinder caps stop where the shapes touch instead of millimeters early.
- Stepping is unchanged.
- A cylinder's GJK distance can still come out slightly high in some poses, so a cast can
  stop a fraction of a millimeter late against one.
- A cast closing on a shape at almost no speed along the GJK normal reports a miss, even
  if it would graze within 0.01 mm.
- A cast takes up to 96 GJK calls per candidate, two near a cylinder. Overlap queries take
  one per candidate or triangle.
