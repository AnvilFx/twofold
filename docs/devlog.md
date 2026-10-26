# Devlog

Notes as I go, newest first. The longer reasoning behind some of this is in
[decisions](decisions/).

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
