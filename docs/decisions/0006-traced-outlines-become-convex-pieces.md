# 6. Traced outlines become convex pieces

Date: 2026-10-27

Status: Accepted

## Context

Flat layers come from images with arbitrary outlines: concave edges, holes, separate
islands. All collision in the engine is between convex shapes.

Options:

- A concave polygon or triangle mesh shape on the planar path.
- Decompose each outline into convex polygons and use them as prisms on one body.

## Decision

Decompose. `alphaPieces` traces the alpha channel with marching squares, simplifies the
contours together so they never touch, bridges holes into their outlines, ear clips, and
merges triangles back into convex pieces with Hertel-Mehlhorn, up to `maxVertices` each.
The pieces become prisms on one body.

## Consequences

- No new shape type and no new collision code. Traced bodies work on both the planar and
  3D paths.
- Holes are real. A ring can hang on a peg through its hole.
- A detailed outline becomes many shapes, so the contact count between two traced bodies
  grows with the number of pieces. `tolerance` and `maxVertices` trade accuracy for piece
  count.
