# 23. Half-spaces are half-planes on the planar path

Date: 2026-10-28

Status: Accepted

## Context

A 2D scene is mostly layers on a floor, and the floor is usually a half-space. `makeProfile`
had no case for half-spaces, so every layer on the floor took the 3D path. There a layer
only a few millimeters thick rests on the floor through SAT and the 3D manifold, which
leaves its contact normals slightly out of the plane on a tilted floor and makes it stop
and settle differently than it does on a static box floor.

The options were:

- Keep the 3D path for half-spaces.
- Cut the half-space with the plane at every depth of the other profile, giving a line
  that moves with depth when the half-space leans.
- Treat it as one line taken at a single depth.

## Decision

A half-space whose plane normal is perpendicular to the planar body's plane normal, within
the alignment tolerance the planar path already uses, becomes a half-plane profile covering
every depth. The line is taken at the middle of the other profile's depth range, where the
contacts are placed anyway. A polygon touches it with its edge facing the line most, up to
two points, and a circle with one point. A half-space leaning further than the tolerance,
or a pair of half-spaces, falls back to 3D.

## Consequences

- Planar boxes, polygons and discs on a half-space floor rest, slide and roll as they do
  on a static box floor, within 0.05 mm, and their normals have no part along the plane
  normal.
- Flat layers and Cutouts step differently. Overlap is about the same, ground contacts
  have about half as many points, and a step costs slightly less.
- Within the tolerance a leaning half-space is off by at most its lean times half the
  layer thickness, which for layers is well under the linear slop.
