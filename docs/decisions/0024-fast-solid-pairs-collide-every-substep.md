# 24. Fast solid pairs collide every substep

Date: 2026-10-28

Status: Accepted. Extends [decision 8](0008-collide-fast-pairs-again-mid-step.md).

## Context

A thin plate hitting the ground at 100 m/s or more came back up with restitution 0, at 1
to 34 m/s, and spun up. Two things caused it:

- Reduction chose the spread among points near the deepest. The corners of the edge that
  hits first won, the plate pivoted on that edge and flipped up.
- Points are anchored at the middle of the gap. A fast pair's speculative points can be
  meters apart, so friction acted half that distance from where the plate landed and
  turned its spin into upward speed.

Choosing the spread among every point predicted to land and anchoring far points on each
body's surface stopped the bounce. On its own it made Cutouts, Flat layers, Meshes and Rain
worse. The anchors alone caused the Cutouts and Flat layers part, and neither half fixed
the bounce without the other.

Options tried:

- Both changes for every pair. Stops the bounce, worse in several scenes.
- Both changes for fast pairs only. Same scenes worse.
- Colliding fast pairs again every substep, alone. Bounces up to 37 m/s.
- Colliding every fast pair again every substep with both changes. Fixes the bounce,
  better almost everywhere, but planar pairs in Cutouts went from 38 frames past 10 mm to
  766.

## Decision

- For a fast pair, reduction chooses the spread among every point predicted to reach the
  surface.
- Points farther apart than the speculative distance are anchored on each body's surface,
  except on planar pairs.
- Fast pairs without a planar body are collided again before every substep after the
  first, with the margin covering what is left of the step. Planar fast pairs keep the
  single refresh halfway through.

## Consequences

Measured over 32 seeds, frames with overlap past 20 / 50 mm and the worst frame:

| Scene | Before | After |
| --- | --- | --- |
| Rain | 2028 / 260, 268 mm | 1284 / 167, 163 mm |
| Hulls | 306 / 17, 71 mm | 56 / 0, 38 mm |
| Meshes | 447 / 9, 86 mm | 184 / 6, 54 mm |
| Sweeper | 13 / 0, 45 mm | 0 / 0, 15 mm |
| Joints | 2619 / 270, 121 mm | 2073 / 180, 67 mm |

- The plate stays under 0.6 m/s up to 500 m/s. At 2000 m/s it comes back at 5 to 6 m/s.
- Tower, Pyramid and Cutouts are unchanged. Dominoes and Flat layers move by tenths of a
  millimeter either way.
- A few more bodies are still awake at the end of Rain and Meshes.
- Steps cost more where fast solid pairs are common: Meshes 46%, Hulls 33%, Sweeper 24%,
  Rain 13%. Scenes without them cost the same.
