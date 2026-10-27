# 1. Planar bodies are locked 3D bodies

Date: 2026-10-26

Status: Accepted

## Context

The engine has to handle scenes where flat layers and free 3D bodies share one world. A
layer falling onto a 3D box has to land on it, and a 3D ball thrown at a stack of layers
has to knock it over. Layers in the same plane should behave the way they would in a 2D
engine: exact circles, no drift out of the plane, no contacts pushing them apart in depth.

Options:

- Run a 2D world and a 3D world side by side and bridge the pairs between them.
- Add a separate 2D body type with its own integration, and special-case every pair that
  mixes the two.
- Keep every body 3D, and make "planar" a set of locks.

## Decision

Every body is a 3D body. `LockPlanar` locks translation along a frame's z axis and
rotation about its x and y axes. Locks are applied by projecting velocities onto the free
axes and inverting inertia only in the free subspace.

A pair whose dynamic bodies are all planar with parallel normals collides through a 2D
path first. Shapes that extrude along the normal become polygons and circles in the plane,
and contacts get an in-plane normal. Anything the 2D path can't handle falls back to the
3D path.

## Consequences

- Mixed pairs need no special code. A planar layer against a free sphere is an ordinary 3D
  contact, and the locks keep the layer in its plane.
- One integrator, one solver, one sleep system.
- A planar body on the 3D path can get contact directions that point along its locked
  axis. The locks keep it in place, but the solver has to cope with those directions. See
  [decision 9](0009-locked-directions-have-no-effective-mass.md).
- Layers in planes that aren't parallel use the 3D path against each other.
