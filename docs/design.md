# Design

This is how Twofold works inside. The README covers what it does and how to use it. The
records in [decisions](decisions/) explain why some of these choices were made.

## Goals

- One world for 2D and 3D bodies. A flat layer and a free 3D body collide with each other
  like any two bodies, and two flat layers in the same plane behave like a 2D engine.
- Deterministic. The same calls in the same order give bit-identical results, and a copied
  world steps exactly like the original.
- Small and readable. No dependencies outside the standard library, no threads, no SIMD.
  Speed matters, but not at the cost of being able to follow a step from start to end.
- An API that does not assume a particular host. Plain value types in, handles out.

## Layout

| File | Contents |
| --- | --- |
| `include/twofold/world.h` | The public API: `World`, descriptors, geometry types, queries |
| `include/twofold/math.h` | `Vec2`, `Vec3`, `Quat`, `Mat3`, `Transform` |
| `src/internal.h` | Internal types shared by the source files |
| `src/world.cpp` | Bodies, shapes, the step, broadphase, sleeping |
| `src/solver.cpp` | Contact preparation, warm starting, the contact solve, restitution |
| `src/collide3d.cpp` | SAT for polyhedra, sphere and half-space pairs, manifold reduction |
| `src/collide2d.cpp` | The planar path: profiles, polygon and circle pairs |
| `src/shape.cpp`, `src/polyhedron.cpp` | Shape building, mass properties, bounds |
| `src/query.cpp` | Ray casts |

## World and handles

`World` holds a pointer to `World::Impl`, which owns every body, shape and contact in flat
vectors. Bodies and shapes are addressed by `BodyId` and `ShapeId`: an index into the
vector plus a generation. Destroying a body bumps its generation and puts the slot on a
free list, so a stale id fails `isValid` instead of reaching whatever reused the slot.

Copying a `World` copies the `Impl` by value. Nothing inside points into anything else, so
a copy is a full snapshot and steps identically. The viewer's timeline is built on this.

## Bodies

A body is static, kinematic or dynamic. Static bodies never move on their own. Kinematic
bodies move at the velocity they are given and are never pushed. `setKinematicTarget` sets
the velocity that reaches a pose in one step, and the step snaps the body to that pose at
the end so rounding never drifts it.

Mass and inertia come from the shapes and their densities, about the combined center of
mass. `setMassProperties` overrides them until `resetMassProperties`.

### Locks

Locks remove degrees of freedom along the axes of the body's lock frame. Each body keeps
two projection matrices, one for linear and one for angular velocity, built from the axes
that are still free. Velocities are projected after every change, so a locked axis never
picks up motion from gravity, contacts or the user.

The world inverse inertia is inverted only in the free subspace (`invertSubspace`). A body
that can only turn about z gets the inverse of its inertia about z, not the z entry of the
full inverse. The inverse mass matrix is the linear projection over the mass.

`LockPlanar` is translation z plus rotation x and y. A planar body is still a normal 3D body
with locks. See [decision 1](decisions/0001-planar-bodies-are-locked-3d-bodies.md).

## Shapes

| Geometry | Stored as |
| --- | --- |
| `Sphere` | A point with a radius |
| `Box`, `Prism` | A polyhedron |
| `Cylinder` | A 24-sided prism in 3D, an exact circle on the planar path |
| `HalfSpace` | Everything below the local xy plane. Static bodies only |

See [decision 5](decisions/0005-cylinders-are-prisms-in-3d.md) for why cylinders are prisms.

Every polyhedron keeps vertices, face loops with planes, and an edge list with the two
faces on each edge. The edge list is what the SAT edge queries walk.

## The step

`World::step(dt)` runs these in order:

1. **Contacts.** Broadphase, then narrowphase for every candidate pair. New contact points
   inherit impulses from the old ones they continue.
2. **Step state.** Every body records its pose at the start of the step, resets its motion
   since that pose (`dp`, `dq`), and refreshes its world inverse mass.
3. **Active contacts.** Contacts with an awake dynamic body take part in the solve.
4. **Substeps.** `substeps` times (4 by default):
   - integrate velocities
   - warm start
   - solve contacts with bias
   - integrate positions
   - relax: solve contacts again without bias
5. **Restitution.** One pass over points that were approaching faster than the threshold.
6. **Clean up.** Forces and torques clear, and kinematic bodies snap to their targets.
7. **Sleep.** Islands are rebuilt and the ones that have been slow long enough fall asleep.

### Broadphase

Sort and sweep on x, rebuilt from scratch every step. Each shape's bounds grow by the
speculative distance, and a moving body's grow by how far it can travel in the step at its
current velocity. A pair is kept when the bounds overlap, the filters agree, at least one
body is dynamic, and at least one is moving: awake if dynamic, or kinematic with a nonzero
velocity. Contacts between bodies that are both still are carried over untouched, which is
what keeps a sleeping pile cheap.

See [decision 4](decisions/0004-sort-and-sweep-rebuilt-every-step.md).

### Narrowphase in 3D

- **Polyhedron pairs** use SAT: face normals of A, face normals of B, then edge pairs. Edge
  pairs are pruned with the Gauss map test, so only pairs that form a face of the Minkowski
  difference are tried. Edge pairs closer to parallel than a small sine are skipped. Their
  axis is unstable and a face axis gives the same answer.
- **Face contacts** clip the incident face against the side planes of the reference face.
- **Edge contacts** use the closest points of the two edges as a single contact.
- **Spheres** against polyhedra, spheres and half-spaces have their own direct tests. So do
  polyhedra against half-spaces.

Every pair produces candidate points with a separation. Speculative points, up to the
margin apart, are kept so the solver can stop an approach before it turns into overlap.

### Manifold reduction

A manifold keeps at most four points. Each candidate gets a predicted separation: its
current separation plus the relative normal velocity at that point times `dt`. The point
predicted deepest goes in first. Among the points within a small band of it, the next goes
to the one farthest away, then the one giving the largest triangle, then the one farthest
outside that triangle. Free slots take the deepest of the rest.

Ranking on predicted separation rather than current separation matters for tumbling
bodies. The vertex about to land is kept even though it is not touching yet.

### The planar path

A pair tries the 2D path first when every dynamic body in it is planar and their plane
normals are parallel. Static and kinematic bodies fit any plane. `makeProfile` turns
each shape into a 2D profile in the plane: a polygon for a box or prism, a circle for a
cylinder, plus the depth range it covers along the normal. It fails for anything that does
not extrude along the normal, and then the pair uses the 3D path.

Profiles whose depth ranges do not overlap never touch, so layers at different depths pass
each other. Otherwise polygons collide with reference and incident edges as in Box2D, and
the contacts are placed at the middle of the shared depth range. The normal lies in the
plane, so a contact can never push two layers apart in depth.

### Contact matching and warm starting

Contacts are stored sorted by a key made of the two shape indices. A new point stores its
position in body A's local frame. If an old point of the same pair is within the
speculative distance of it, the new point takes that point's normal impulse and the
tangent part of its friction impulse. This is all the warm starting there is, and it is
enough for stacks to hold.

## Solver

The solver is the soft step from Box2D v3. Contacts are springs with a stiffness in hertz
and a damping ratio. That gives a bias velocity, a mass scale and an impulse scale for each
contact. Defaults are 30 Hz and a damping ratio of 10. The stiffness is capped at a quarter
of the substep rate, and contacts with a static, kinematic or sleeping body use twice the
stiffness. Pushout from overlap is capped at `maxPushoutSpeed`, 3 m/s by default. See
[decision 2](decisions/0002-soft-step-contact-solver.md).

Positions are never corrected directly. Each body keeps how far it has moved since the
start of the step (`dp`) and how far it has turned (`dq`). A contact's current separation
is its separation at the start plus how much `dp` and `dq` have closed it. So the solver
sees separation change between substeps without collision running again.

Friction is a two-axis cone, clamped by the point's normal impulse times the mixed
friction. The default mixing is the geometric mean for friction and the maximum for
restitution. Both can be replaced with `setFrictionMix` and `setRestitutionMix`.

## Sleeping and waking

A body's sleep timer runs while its speed, counting rotation at its farthest point, stays
under `sleepThreshold`. Islands are rebuilt every step with union-find over touching
contacts between awake dynamic bodies. An island sleeps when every body in it has been slow
for `timeToSleep`. A moving kinematic body keeps whatever touches it awake.

Waking floods through touching contacts from the bodies that need it:

- A dynamic body wakes its whole island.
- A sleeping body touched by a moving body wakes during the narrowphase.

Removing a shape or body wakes the bodies that touched it.

## Units and tolerances

Every tolerance is in meters and scales with `WorldSettings::lengthUnit`. Linear slop is
5 mm. The speculative distance is four times that, 20 mm. Pushout speed, the restitution
threshold and the sleep threshold scale the same way. A scene in centimeters sets
`lengthUnit` to 100 and behaves exactly like the same scene in meters.

## Determinism

The step is single threaded and never iterates in an order that depends on memory
addresses or hashing. Candidate pairs and contacts are sorted by shape index, and bodies
are visited by index. The tests run the same scene twice and compare every pose bit for
bit, and check that a copied world steps like the original. The viewer's verify command
replays a whole session from a fresh scene and compares every frame.

See [decision 3](decisions/0003-single-threaded-and-deterministic.md).

## Known limits

- No joints, triangle meshes, shape casts, overlap queries or contact callbacks yet.
- The broadphase rebuilds every step. It is fine for a few thousand shapes and will need a
  tree after that.
- A planar body on a half-space goes through the 3D path. It behaves correctly but does
  not get the exact circles of the planar path.
- No continuous collision. A body that speeds up within a step, from a hit or a large
  force, can pass through a thin shape.
