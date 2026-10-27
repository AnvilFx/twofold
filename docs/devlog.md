# Devlog

Notes as I go, newest first. The longer reasoning behind some of this is in
[decisions](decisions/).

## 2026-10-27 02:33

Viewer came up tiny on a high DPI monitor. It now scales with the monitor it opens on.

A box lying across the edge of a slab slowly sank into the edge. SAT picked an edge pair
and an edge pair only ever made one contact point, so the box pivoted on it. Edge contacts
now clip the face of B against the face of A, seen along the edge axis, which gives two
points when an edge lies across a face. Also, two small faces deep inside each other can
miss each other entirely when clipped, and the pair lost its contact. The deepest vertex
stands in for that case now.

## 2026-10-27 01:27

GJK distance and conservative advancement for time of impact, for continuous collision.
The simplex solve uses Ericson's closest point on a triangle. Nothing calls it yet.

Continuous collision. Speculative contacts already stop anything moving at the speed the
step started with. What gets through is a body that speeds up within the step: knocked by
another body, or pushed by a huge force. After the solve, any shape that moved more than
the step expected plus half its inner extent gets swept against everything near its path.

The sweep is done twice. First a small sphere at the shape's core. If the core never
reaches the other shape, nothing tunneled, and the body stays where the solve put it. Only
if the core gets through does the full shape get swept to find the first touch. Without the
core check, every graze rewound the body and fast bodies stuttered.

Dynamic bodies that hit each other rewind together. Rewinding one into a body that didn't
rewind just makes a new overlap.

Thin planar layers would always count as fast, since their inner extent is the layer
thickness. For planar bodies the faces along the plane normal are left out of it.

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

Added a Cutouts scene: rings, gears, stars and hooks traced from images, falling on pegs
and into a traced bowl. A ring catches on a peg through its hole, which was the point.

Cutouts turned up a bad case in the 2D polygon test. Two pieces passing corner to corner,
8 cm apart, got contacts saying they were touching. The incident edge lay entirely past
the end of the reference edge, and the clip extrapolated it into points on neither shape.
Separated polygons nearest at two corners now get one point on the line between the
corners. An incident edge past the reference edge uses its deepest vertex instead of
being clipped.

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
