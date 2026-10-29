# 29. Planar pairs are solved in 2D

Date: 2026-10-29

Status: Accepted. Replaces [decision 1](0001-planar-bodies-are-locked-3d-bodies.md).

## Context

Decision 1 made every body a 3D body and planar a set of locks. The 2D path only built the
manifold. Everything around it ran the 3D way: profiles were made again from the shapes'
world poses for every pair every step, contact rows had 3D lever arms and a 3x3 inverse
inertia on each side, and a turned body had its world inertia rebuilt every substep. Cutouts,
a scene of nothing but layers, took 0.44 ms a step, about three times what Jolt takes for it
with its bodies locked to the plane. That is no better than locking axes in a 3D engine.

Options:

- A separate 2D world, bridged to the 3D one for mixed pairs. Two integrators, two sleep
  systems, and every pair between them special.
- A 2D body type with its own state. Every place that reads a body would branch on it.
- Keep one body type, and give pairs whose bodies all move in one plane a 2D collision and
  a 2D solve.

## Decision

Bodies stay 3D bodies with locks, so integration, sleeping, joints and every pair that
mixes layers with free bodies work as before.

A pair whose dynamic bodies are planar with parallel normals collides in 2D. Each shape
keeps its profile in its body's frame, made once for the normal as the body sees it and
made again only when that changes. Each body is placed in the plane once per round of
collisions: its position and turn in the plane's 2D frame and its depth along the normal.
A pair then brings b's profile over by one turn and one shift, rejects pairs whose boxes
are further apart than the margin, and runs the polygon, circle and half-plane tests.

A contact whose moving bodies are all planar in one plane, with its normal in that plane,
is solved in 2D. Each point has a scalar arm per direction, the turn a unit impulse gives
about the plane normal, and each body a scalar inverse inertia about that normal. The rows
work on the bodies' in-plane velocity and spin. Contacts holding no impulse skip warm
starting, and a point with no normal and no friction impulse skips its friction row.
Everything else, including a layer against a free 3D body, uses the 3D rows.

## Consequences

- Cutouts goes from 0.44 to 0.125 ms a step, and Flat layers from 0.067 to 0.044. Run
  against Jolt in the same harness, Cutouts takes 0.13 ms a frame to Jolt's 0.145, and its
  worst overlap is 17 mm to Jolt's 141. Over 48 seed sets Cutouts overlaps a little less
  than before, and Flat layers over 24 is the same within noise.
- A layer falling on a 3D box, or a 3D ball hitting a stack of layers, is still an ordinary
  3D contact. The locks keep the layer in its plane.
- Shapes hold a profile and bodies a pose, which the 3D path never reads.
- A planar body's world inverse inertia is no longer rebuilt as it turns, since turning
  about the plane normal leaves it unchanged.
- Layers in planes that aren't parallel still use the 3D path against each other.
