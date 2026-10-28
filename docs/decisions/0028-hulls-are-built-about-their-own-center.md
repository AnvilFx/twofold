# 28. Hulls are built about their own center

Date: 2026-10-28

Status: Accepted

## Context

A hull's face planes came from absolute coordinates. For a cloud 10000 units from the
origin, a float only resolves about 1 mm, and cross products of vertices that far out
lose most of their digits. Face normals were off by up to 0.012 and planes by 3 to 6 mm,
and at 100000 units quickhull dropped points: a 64 vertex sphere came out with 35. A sphere
hull resting on a slab hull, both built that far from their bodies, came to rest at a
different height than the same pair at the origin.

Options:

- Document it and ask for hulls built around the body. Easy to miss, and a level made in
  world coordinates is built exactly that way.
- Double precision in the builder. Every point gets converted, and the polyhedron keeps
  float planes anyway.
- Build about a nearby center and give back the original points.

## Decision

`computeHull` moves the points by the center of their bounds before building and gives the
original points back in the mesh. The tolerance it builds to is 1e-4 of the cloud's size,
with no allowance for distance from the origin. The checks for flat, collinear and single
point clouds, and the slab they get, still use the tolerance that grows with distance.

A polyhedron's face normals come from Newell's method about its centroid instead of the
origin.

## Consequences

- At 10000 and 100000 units, hulls keep every vertex and their normals are within 1e-6 of
  the exact planes, as at the origin. Planes are as accurate as the vertices allow, about 2 mm at
  10000 units. A test drops a sphere hull on a slab hull 10000 units from their bodies and
  checks it rests as it does at the origin.
- Clouds near the origin give the same vertices and the same face counts. In the bench
  scenes, face normals move by up to 5e-6 and offsets by up to 4e-6, and faces may come
  out in another order. The Hulls, Meshes, Rain, Flat layers and Cutouts hashes change
  for that reason.
