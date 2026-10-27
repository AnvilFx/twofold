# 9. Locked directions have no effective mass

Date: 2026-10-27

Status: Accepted

## Context

A planar box sliding at 3 m/s on a half-space slid 27 m. On a box floor it stopped in
0.84 m.

A planar body on a half-space goes through the 3D path. The half-space normal is a
rotated z axis and has a little rounding in it. One friction tangent then lies almost
entirely along the planar body's locked depth axis. Its response is not zero but around
1e-15 times the inverse mass, so its effective mass is huge. The friction impulse along
that tangent uses up the whole friction cone, and the in-plane direction gets almost none.

Options:

- Give half-spaces a profile on the planar path. Fixes half-spaces but not any other pair
  that ends up on the 3D path with a nearly locked direction.
- Project contact tangents into the plane for planar bodies. Needs special cases for pairs
  with different or non-parallel locks.
- Treat any response below a small fraction of the inverse masses as zero.

## Decision

`effectiveMass` returns zero when the response is below 1e-6 times the sum of the two
inverse masses.

## Consequences

- Works for any locked axis on any path, not just planar bodies on half-spaces.
- The planar box now stops in 0.83 m.
- Bodies that used to glide now roll. A traced ring with a bar through it now rolls between
  the bar ends before it sleeps, the same as it does on a box floor.
- A real response this small would need a mass ratio around a million to one, far beyond
  what the solver handles anyway.
