# 21. Cylinders are exact in 3D

Date: 2026-10-27

Status: Accepted

## Context

As 24-sided prisms, cylinders rolled over flat sides. A cylinder rolling on the ground
bobbed by almost 1% of its radius and lost speed at every side, a third of it within a
second on a mesh. Wheels and barrels are common enough that this shows.

GJK and EPA alone give one point and a normal that wanders for flat contacts, which is
most of what a cylinder does: a cap on a floor, a side on a box face, a can on a can.

## Decision

Cylinders keep their exact shape. A pair with a cylinder tests feature axes built from the
flat and straight parts of both shapes: face normals, cylinder axes, and edge and axis
crosses. Those are exact whenever a flat or straight feature pair touches. Where a rim
meets an edge or another rim, none of them is, so the GJK direction between the shapes and
the one between their rounded cores are tried as well. A GJK direction only wins by a
clear margin, so resting contacts keep their feature axis.

Contacts come from clipping the two features along the chosen normal. A cap is 24 rim
points, starting at the deepest one when it is tilted.

GJK treats a nearly flat tetrahedron face as facing the origin. With curved supports, such
a face could otherwise turn a few millimeters of overlap into a contact 200 mm deep.

## Consequences

- Cylinders roll at a constant height. They keep their speed on boxes and half-spaces,
  and lose under 2% of it over three seconds on a mesh.
- Ray casts hit the true surface.
- Two cylinders spinning fast into each other reach 16 mm of overlap where the prisms
  reached 10 mm. The rim points of a cap miss some of a spinning rim within a step.
- Without rolling resistance, cylinders on the ground roll until something stops them, so
  scenes full of them sleep later.
- The planar path still treats a cylinder along the plane normal as a circle.
