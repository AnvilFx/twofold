# 5. Cylinders are prisms in 3D

Date: 2026-10-26

Status: Accepted

## Context

Exact cylinder collision in 3D needs either GJK and EPA for every cylinder pair, or a set
of special cases for cylinder against box, cylinder against cylinder and so on. The rest
of the 3D path is SAT over polyhedra.

In 2D, a circle is cheap and exact, and wheels rolling in a plane are a common case.

## Decision

A cylinder is a 24-sided prism on the 3D path. On the planar path, a cylinder whose axis
is along the plane normal is an exact circle.

## Consequences

- One 3D path for every shape except spheres and half-spaces.
- A cylinder rolling in 3D rolls over 24 flat sides. Close enough for most scenes, and
  visible if you look for it.
- A wheel in a 2D scene rolls exactly.
