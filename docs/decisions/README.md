# Decisions

Design decisions worth explaining. Each record says what the problem was, what else was
considered, and what follows from the choice. When a decision changes, a new record
replaces the old one instead of editing it.

1. [Planar bodies are locked 3D bodies](0001-planar-bodies-are-locked-3d-bodies.md)
2. [Soft step contact solver](0002-soft-step-contact-solver.md)
3. [Single threaded and deterministic](0003-single-threaded-and-deterministic.md)
4. [Sort and sweep, rebuilt every step](0004-sort-and-sweep-rebuilt-every-step.md), replaced by 16
5. [Cylinders are prisms in 3D](0005-cylinders-are-prisms-in-3d.md), replaced by 21
6. [Traced outlines become convex pieces](0006-traced-outlines-become-convex-pieces.md)
7. [Continuous collision only past the speculative margin](0007-continuous-collision-past-the-speculative-margin.md)
8. [Collide fast pairs again mid-step](0008-collide-fast-pairs-again-mid-step.md)
9. [Locked directions have no effective mass](0009-locked-directions-have-no-effective-mass.md)
10. [Moved anchors wake what they overlap](0010-moved-anchors-wake-what-they-overlap.md)
11. [Impact contacts get a second pass](0011-impact-contacts-get-a-second-pass.md)
12. [Old contact points share their impulse](0012-old-contact-points-share-their-impulse.md)
13. [The relax pass takes friction in reverse](0013-relax-pass-takes-friction-in-reverse.md)
14. [Edge queries pair edges near the overlap](0014-edge-queries-pair-edges-near-the-overlap.md)
15. [Pairs that barely moved keep their contacts](0015-pairs-that-barely-moved-keep-their-contacts.md)
16. [A bounds tree with persistent pairs](0016-bounds-tree-with-persistent-pairs.md)
17. [Soft joints with a settle pass](0017-soft-joints-with-a-settle-pass.md)
18. [Limits and motors move the whole joint](0018-limits-and-motors-move-the-whole-joint.md)
19. [Meshes collide one triangle at a time](0019-meshes-collide-one-triangle-at-a-time.md)
20. [A mesh pair has a contact per facing](0020-a-mesh-pair-has-a-contact-per-facing.md)
21. [Cylinders are exact in 3D](0021-cylinders-are-exact-in-3d.md)
22. [Rolling resistance follows the spin](0022-rolling-resistance-follows-the-spin.md)
23. [Half-spaces are half-planes on the planar path](0023-half-spaces-are-half-planes-on-the-planar-path.md)
24. [Fast solid pairs collide every substep](0024-fast-solid-pairs-collide-every-substep.md)
25. [Queries share one touch test](0025-queries-share-one-touch-test.md)
26. [Events follow the solve, and sensors test overlap](0026-events-follow-the-solve-and-sensors-test-overlap.md)
