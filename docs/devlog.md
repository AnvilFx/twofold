# Devlog

Notes as I go, newest first. The longer reasoning behind some of this is in
[decisions](decisions/).

## 2026-10-27 00:50

Convex hulls from point clouds. Quickhull, stopping at a vertex limit, always taking the
farthest point out first, so a hull with a low limit still keeps its corners.
Coplanar triangles merge back into polygon faces so SAT doesn't test the same plane twice.
Flat or collinear clouds get thickened into a thin slab so they still have volume.

Outlines from alpha images. Marching squares with interpolated crossings, Douglas-Peucker
to cut the point count, then corners moved back to where the straight runs meet, since
tracing rounds them off. Contours are simplified together so a thin gap between two
outlines can't close up. Holes get bridged into their outline, then ear clipping and
Hertel-Mehlhorn to merge triangles back into convex pieces of at most 8 vertices.

## 2026-10-26 23:07

First commit. Spheres, boxes, cylinders, prisms and half-spaces, a soft step solver with
substeps, sleeping islands, locks, and the planar path for layers that share a plane.
Then the debug viewer, with a recorded timeline that rewinds and replays exactly.

The viewer showed problems right away:

- A cylinder on its side wouldn't settle. Reduction kept the four deepest points, and those
  could all sit at one end of the face it rests on. Reduction now ranks by separation
  predicted after the step and picks the spread among points close to the deepest.
- Boxes sliding over the edge of a slope sank in. Nearly parallel edge pairs give an
  unstable axis, so those are left to the face queries now.
- Rain, Sweeper and Flat layers spawned bodies inside each other. Fixed the scenes.
- Moving a floor with `setTransform` left the bodies on it floating. Waking only started
  from dynamic bodies. A static or kinematic body now wakes what it touches.
