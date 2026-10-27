# 19. Meshes collide one triangle at a time

Date: 2026-10-27

Status: Accepted

## Context

Terrain and imported models need to be collision geometry as they are, concave and open.
A body sliding over a mesh must not catch on the edges between triangles, even though
each triangle on its own has edges that stick out.

Options:

- Convex decomposition at load time. Approximate, slow to build, and still leaves seams
  between the pieces.
- Each triangle as a thin hull. Reuses everything, but a box sliding across a flat floor
  hits the side faces of the next triangle.
- Each triangle as a two-sided face, with precomputed flags saying which edges stick out,
  as in Bullet and Jolt.

## Decision

A `TriangleMesh` is built once: vertices at the same position merge, degenerate triangles
drop, and a bounding volume tree goes over the rest. Each triangle edge gets a flag per
side. An edge has a normal of its own on a side only where the neighbor folds away from
that side by more than about 5 degrees, or where it is open or shared by more than two
triangles.

A body collides with each triangle near it using the same SAT, clipping and sphere tests
as two convex shapes. The triangle is a polyhedron with a front and a back face, so both
sides collide. When the contact normal is off the face and the closest feature has no
normal of its own on that side, the contact is made again against the face. If the face
then claims more overlap than the shapes really have, the contact belongs to a neighbor
and is dropped.

Manifolds from triangles facing the same way merge and reduce to four points.

## Consequences

- Boxes slide and balls roll across seams of flat and gently curved meshes without
  catching.
- Meshes are only for static and kinematic bodies. They don't collide with each other or
  with half-spaces.
- Without the overlap check, a body resting on one side of a ridge got a contact from the
  far triangle pointing down into the mesh, up to 0.8 m deep. The check keeps the deepest
  contact on a hilly mesh under 50 mm.
- Every triangle is collided on its own, so a body over many small triangles costs more
  than one over a few large ones.
