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
| `include/twofold/mesh.h` | `TriangleMesh` |
| `include/twofold/outline.h` | Alpha tracing, contour simplification, convex decomposition |
| `src/internal.h` | Internal types shared by the source files |
| `src/world.cpp` | Bodies, shapes, the step, broadphase, sleeping, continuous collision |
| `src/tree.h`, `src/tree.cpp` | The dynamic AABB tree behind the broadphase and queries |
| `src/solver.cpp` | Contact preparation, warm starting, the contact solve, restitution |
| `src/joint.cpp` | Joint warm starting, the joint solve and the settle pass |
| `src/collide3d.cpp` | SAT for polyhedra, sphere, half-space and mesh pairs, manifold reduction |
| `src/collide2d.cpp` | The planar path: profiles, polygon and circle pairs |
| `src/distance.cpp` | GJK distance, sweeps, time of impact |
| `src/hull.cpp` | Quickhull with a vertex limit |
| `src/mesh.cpp` | Mesh building, edge flags, the triangle tree, mesh rays |
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
| `Mesh` | A shared `TriangleMesh`. Static and kinematic bodies only |

See [decision 5](decisions/0005-cylinders-are-prisms-in-3d.md) for why cylinders are prisms.

A `TriangleMesh` merges vertices at the same position and drops triangles whose corners
are nearly collinear. Its tree splits at the median of the triangle centers along their
widest axis, down to four triangles a leaf. Each triangle keeps its normal and six edge
flags, one per edge and side. An edge has a normal of its own on a side when the neighbor
folds away from that side by more than a sine of 0.087, and on both sides when it is open
or shared by more than two triangles. Meshes have no mass and don't collide with each
other or with half-spaces.

Every polyhedron keeps vertices, face loops, face planes, its bounds, and an edge list with
the two faces on each edge. Each edge also stores its length and its arc on the Gauss map,
the cross of its two face normals. The edge list is what the SAT edge queries walk.

`computeHull` is quickhull. It stops adding points once the hull reaches `maxVertices`,
always taking the point farthest outside the current hull, so a limited hull keeps the
corners that matter most. Faces whose normals are within a cosine of 0.9998 merge into one
polygon. A flat, collinear or single-point cloud is thickened into a thin slab first, so
every hull has volume.

## The step

`World::step(dt)` runs these in order:

1. **Contacts.** Broadphase, then narrowphase for every candidate pair. A pair that has
   barely moved keeps its points. New contact points inherit impulses from the old ones they
   continue.
2. **Step state.** Every body records its pose at the start of the step, resets its motion
   since that pose (`dp`, `dq`), and refreshes its world inverse mass.
3. **Active contacts.** Contacts with an awake dynamic body take part in the solve. Pairs
   whose relative motion over the step exceeds the speculative distance are also listed as
   fast.
4. **Substeps.** `substeps` times (4 by default):
   - halfway through, fast pairs are collided again (see below)
   - integrate velocities
   - warm start joints, then contacts
   - solve joints with bias, then contacts, then again for impact contacts
   - integrate positions
   - relax: solve joints and contacts again without bias, then again for impact contacts
5. **Restitution.** One pass over points that were approaching faster than the threshold.
   Then the settle pass moves jointed bodies to close what joint error is left.
6. **Clean up.** Forces and torques clear, and kinematic bodies snap to their targets.
7. **Continuous collision.** Bodies that moved well past what their speculative contacts
   covered are swept and moved back.
8. **Sleep.** Islands are rebuilt and the ones that have been slow long enough fall asleep.

### Broadphase

Every shape except half-spaces has a leaf in a dynamic AABB tree. A leaf covers the
shape's bounds grown by the speculative distance and a 100 mm margin, and moves in the
tree only when the shape leaves it. Half-spaces are kept in a list of their own.

The world keeps the sorted pairs of shapes whose leaves overlap, on different bodies with
at least one dynamic. Each step, shapes whose leaves moved find their new pairs, and pairs
whose leaves no longer overlap are dropped.

From those pairs and the half-spaces, a candidate goes to the narrowphase when the bounds
overlap, the filters agree, and at least one body is moving: awake if dynamic, or
kinematic with a nonzero velocity. Here each shape's bounds grow by the speculative
distance, and a moving body's grow by how far it can travel in the step at its current
velocity. Contacts between bodies that are both still are carried over untouched, which
is what keeps a sleeping pile cheap.

Ray casts, the sweeps of continuous collision and waking query the same tree.

See [decision 16](decisions/0016-bounds-tree-with-persistent-pairs.md).

### Narrowphase in 3D

- **Polyhedron pairs** use SAT in A's local frame: face normals of A, face normals of B,
  then edge pairs. Each face query starts with the face toward the other shape and tests
  the other shape's deepest vertex first. Only edges reaching into the overlap of the two
  bounds are paired, and when that skips any edge a GJK distance check confirms the shapes
  are within the margin. See
  [decision 14](decisions/0014-edge-queries-pair-edges-near-the-overlap.md). Edge pairs are
  pruned with the Gauss map test, so only pairs that form a face of the Minkowski
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
- **Meshes** collide with each triangle whose bounds reach the other shape, in the mesh's
  frame. See below.

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

### Meshes

A triangle is a polyhedron with two faces, front and back, so the face queries, face
contacts and edge contacts above apply to it unchanged. Its edge query is its own: the arc
of a triangle edge on the Gauss map is the half circle from the front normal through the
edge's outward direction to the back normal. A sphere takes the closest point on the
triangle.

When a contact's normal is off the triangle's normal, the triangle feature closest to its
deepest point is found. If that edge or corner has no normal of its own on the side the
body is on, the contact is made again as a face contact on that side. If any of its points
is deeper than SAT found the shapes to overlap, the contact belongs to a neighbor and is
dropped. A sphere's contact on such a feature is dropped unless the sphere's center lies
over the face. See [decision 19](decisions/0019-meshes-collide-one-triangle-at-a-time.md).

Contacts from triangles whose normals are within a cosine of 0.996 merge under the normal
of the deepest one and are reduced to four points. Each group becomes its own contact.

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
position in body A's local frame and continues the nearest old point of the same pair
within the speculative distance. It takes that point's normal impulse and the tangent part
of its friction impulse. When several new points continue the same old point, they split
its impulse evenly, so a manifold that changes shape doesn't push harder than before. See
[decision 12](decisions/0012-old-contact-points-share-their-impulse.md).

A pair whose bodies moved less than `contactReuseDistance` and turned less than
`contactReuseAngle` relative to each other since it was last collided skips the
narrowphase. It keeps its points and impulses, measures each separation again from where
the point sits on each body, and turns its normal with body A. See
[decision 15](decisions/0015-pairs-that-barely-moved-keep-their-contacts.md).

A mesh pair can have several contacts, ordered by a part index after the key. They are
reused only together, when all of them still fit. A new contact takes its impulses from
the old contact of the pair with the closest normal. See
[decision 20](decisions/0020-a-mesh-pair-has-a-contact-per-facing.md).

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

Each point direction's lever arms, velocity change per unit impulse and effective mass are
built once per step. The angular responses are updated after each position integration,
so a turned body responds through its turned inertia. The effective masses stay as they
were at the start of the step.

Friction is a two-axis cone, clamped by the point's normal impulse times the mixed
friction. The biased pass solves friction points in order, the relax pass in reverse. With
one fixed order small stacks rock forever and never sleep. See
[decision 13](decisions/0013-relax-pass-takes-friction-in-reverse.md).

A body in a pair that closes faster than the restitution threshold has all its contacts
solved twice in each pass. One pass doesn't carry an impact through a stack within a
substep, so a heavy body landing on a light one would crush it into the floor. Resting
bodies pay nothing for this. See
[decision 11](decisions/0011-impact-contacts-get-a-second-pass.md). The default mixing is the geometric mean for friction and the maximum for
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
base orientation of every body reset. New points then take impulses from old points the
same way as after the narrowphase, and the restitution data of the nearest old point. A
mesh pair is collided once, and each of its contacts takes the new group whose normal is
closest to its own.

See [decision 8](decisions/0008-collide-fast-pairs-again-mid-step.md).

## Joints

A joint stores its anchors in each body's local frame, and for a hinge or slider its axis
in both. An unset body A is the world, a static body at the origin. Values such as the
hinge angle and slider translation are body B relative to A, zero in the joined pose.

| Joint | Rows |
| --- | --- |
| Point | A 3D point row |
| Fixed | A 3D turn row, then the point row |
| Hinge | Motor, limits, a 2D turn row across the axis, then the point row |
| Slider | Motor, limits, a 3D turn row, then a 2D point row across the axis |
| Distance | One row along the line between the anchors, rigid, sprung or limited |

The point and turn rows are soft, at `jointHertz` and `jointDampingRatio`, capped at half
the substep rate. Limits only push. Like speculative contacts, they start acting before
the limit is reached, so a fast swing stops at it instead of past it. A motor is clamped
to its maximum torque or force times the substep. Every row warm starts from the impulse it
ended the last step with.

A hinge's motor and limits act about the anchor, not the body's center. Their impulse is
a turn about the axis plus the point impulse that turn needs to keep the anchors
together. A slider's motor and limits likewise carry the turn that keeps the bodies from
rotating. Without this a door bounces off its limit. See
[decision 18](decisions/0018-limits-and-motors-move-the-whole-joint.md).

The relax pass solves every row rigid. After restitution, two settle iterations move the
bodies of each active joint directly, turn first, then any exceeded limit, then the point
or slide, without changing velocities. See
[decision 17](decisions/0017-soft-joints-with-a-settle-pass.md).

A joint is active when one of its bodies is an awake dynamic body. Joints tie their
dynamic bodies into one island, so they sleep and wake together, and a moving kinematic
end keeps the other body awake. Bodies joined directly don't collide unless
`collideConnected` is set. Changing a joint's limits, motor or spring wakes it.

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
so none is left inside another. Against a mesh, the full shape is swept against each
triangle under the path. The triangles are found in the mesh's frame at its poses at both
ends of the step.
Meshes are never fast themselves. Time is lost for a moved body. The next step carries on
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

- No shape casts, overlap queries or contact callbacks yet.
- Meshes can't be dynamic. Use hulls or several convex shapes for a moving concave body.
- A light chain holding a much heavier body stretches. Twelve links holding a ball about
  200 times heavier than each link rest 0.17 m longer in total.
- A planar body on a half-space goes through the 3D path. It behaves correctly but does
  not get the exact circles of the planar path.
- Contacts are springs, so a stack sags under its own weight. A three box tower rests
  with 7 mm of overlap between the bottom two boxes. More substeps or a higher
  `contactHertz` make it stiffer.
- Fast spinners that end up deep inside each other under a pile can stay 150 to 220 mm
  deep for 8 to 12 frames before they recover. This also happens without the mid-step
  refresh. Pushing deep pairs out harder fixes it but keeps some hull piles from ever
  sleeping.
