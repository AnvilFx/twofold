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
| `include/twofold/hull.h` | `computeHull` |
| `include/twofold/outline.h` | Alpha tracing, contour simplification, convex decomposition |
| `src/internal.h` | Internal types shared by the source files |
| `src/world.cpp` | Bodies, shapes, the step, broadphase, sleeping, continuous collision |
| `src/solver.cpp` | Contact preparation, warm starting, the contact solve, restitution |
| `src/collide3d.cpp` | SAT for polyhedra, sphere and half-space pairs, manifold reduction |
| `src/collide2d.cpp` | The planar path: profiles, polygon and circle pairs |
| `src/distance.cpp` | GJK distance, sweeps, time of impact |
| `src/hull.cpp` | Quickhull with a vertex limit |
| `src/outline.cpp` | Marching squares, simplification, ear clipping, Hertel-Mehlhorn |
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
| `Hull` | A polyhedron from `computeHull` |
| `HalfSpace` | Everything below the local xy plane. Static bodies only |

See [decision 5](decisions/0005-cylinders-are-prisms-in-3d.md) for why cylinders are prisms.

Every polyhedron keeps vertices, face loops with planes, and an edge list with the two
faces on each edge. The edge list is what the SAT edge queries walk.

`computeHull` is quickhull. It stops adding points once the hull reaches `maxVertices`,
always taking the point farthest outside the current hull, so a limited hull keeps the
corners that matter most. Faces whose normals are within a cosine of 0.9998 merge into one
polygon. A flat, collinear or single-point cloud is thickened into a thin slab first, so
every hull has volume.

## The step

`World::step(dt)` runs these in order:

1. **Contacts.** Broadphase, then narrowphase for every candidate pair. New contact points
   inherit impulses from the old ones they continue.
2. **Step state.** Every body records its pose at the start of the step, resets its motion
   since that pose (`dp`, `dq`), and refreshes its world inverse mass.
3. **Active contacts.** Contacts with an awake dynamic body take part in the solve. Pairs
   whose relative motion over the step exceeds the speculative distance are also listed as
   fast.
4. **Substeps.** `substeps` times (4 by default):
   - halfway through, fast pairs are collided again (see below)
   - integrate velocities
   - warm start
   - solve contacts with bias
   - integrate positions
   - relax: solve contacts again without bias
5. **Restitution.** One pass over points that were approaching faster than the threshold.
6. **Clean up.** Forces and torques clear, and kinematic bodies snap to their targets.
7. **Continuous collision.** Bodies that moved well past what their speculative contacts
   covered are swept and moved back.
8. **Sleep.** Islands are rebuilt and the ones that have been slow long enough fall asleep.

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
  If the faces miss each other along the normal, which happens to small faces in a deep
  overlap, the deepest vertex stands in so the pair never loses its contact.
- **Edge contacts** clip the face of B nearest the axis against the face of A, seen along
  the edge axis. An edge lying across a face rests on two points instead of one. When that
  clip is empty, the closest points of the two edges make a single contact.
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

Two cases in the polygon test need care. Separated polygons whose nearest features are two
vertices get one contact along the line between them, not a clipped edge. An incident edge
that lies entirely past the end of the reference edge cannot be clipped, so its deepest
vertex is used instead.

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

A direction the bodies cannot move in, such as the depth axis of a planar body, has an
effective mass of zero. Rounding can leave such a direction with a tiny nonzero response,
which would turn into a huge effective mass. Any response smaller than a millionth of the
two inverse masses is treated as zero. See
[decision 9](decisions/0009-locked-directions-have-no-effective-mass.md).

### Mid-step refresh

Contact points come from the start of the step. A body spinning fast carries those points
away from where the overlap really is, and the solver pushes on points that no longer mean
anything. Pairs listed as fast at the start of the step are collided again halfway through
the substeps.

Before that, every active contact is rebased. Its anchors become the bodies' current lever
arms, its stored separation absorbs the distance already moved, and `dp`, `dq` and the
base orientation of every body reset. New points then take impulses from old points
within the speculative distance, and the restitution data of the nearest old point.

See [decision 8](decisions/0008-collide-fast-pairs-again-mid-step.md).

## Continuous collision

Speculative contacts cover motion at the velocities a step starts with. A body that speeds
up within the step, from a hit or a large force, can still pass through a thin shape. After
the solve, a shape counts as fast when it moved and turned further than the step expected
plus half its inner extent.

Fast shapes are swept against everything near their path in two stages:

1. A core sphere, at the shape's centroid with half its inner extent as radius, is swept
   with conservative advancement. If the core never reaches the other shape, nothing
   tunneled and the body keeps its pose.
2. If the core hits, the full shape is swept to find where it first touched. The body moves
   back to that time on its sweep.

Dynamic bodies that hit each other move back together, to the earliest time in the group,
so none is left inside another. Time is lost for a moved body. The next step carries on
with the velocity the solve gave it.

For a planar body, the inner extent ignores the faces facing along the plane normal. A
thin layer is wide in the plane, and that width is what keeps it from tunneling.

See [decision 7](decisions/0007-continuous-collision-past-the-speculative-margin.md).

## Sleeping and waking

A body's sleep timer runs while its speed, counting rotation at its farthest point, stays
under `sleepThreshold`. Islands are rebuilt every step with union-find over touching
contacts between awake dynamic bodies. An island sleeps when every body in it has been slow
for `timeToSleep`. A moving kinematic body keeps whatever touches it awake.

Waking floods through touching contacts from the bodies that need it:

- A dynamic body wakes its whole island.
- A static or kinematic body wakes the dynamic bodies touching it, without spreading past
  them, and every sleeping body whose bounds overlap its own. For a half-space the overlap
  test is against its plane.
- A sleeping body touched by a moving body wakes during the narrowphase.

Removing a shape or body wakes the dynamic bodies that touched it.

See [decision 10](decisions/0010-moved-anchors-wake-what-they-overlap.md).

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

## Outlines

`alphaPieces` turns an alpha image into convex polygons ready to use as prisms:

1. **Trace.** Marching squares on the alpha channel, with crossings interpolated between
   pixel centers. Holes wind the opposite way to the outlines around them.
2. **Simplify.** Douglas-Peucker on each closed contour. Corners then move to where the
   fitted straight runs beside them meet, which restores corners that tracing rounded off.
   Contours are simplified together, and any edges that touch or cross are refined until
   none do.
3. **Decompose.** Holes are bridged into their outline with a zero-width cut. The result
   is ear clipped, then Hertel-Mehlhorn removes diagonals, longest first, wherever both ends
   stay convex and the piece stays within `maxVertices`.

See [decision 6](decisions/0006-traced-outlines-become-convex-pieces.md).

## Known limits

- No joints, triangle meshes, shape casts, overlap queries or contact callbacks yet.
- The broadphase rebuilds every step. It is fine for a few thousand shapes and will need a
  tree after that.
- A planar body on a half-space goes through the 3D path. It behaves correctly but does
  not get the exact circles of the planar path.
- Fast spinners that end up deep inside each other under a pile can stay 150 to 220 mm
  deep for 8 to 12 frames before they recover. This also happens without the mid-step
  refresh. Pushing deep pairs out harder fixes it but keeps some hull piles from ever
  sleeping.
