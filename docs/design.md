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
| `src/world.cpp` | Bodies, shapes, the step, broadphase, sleeping, continuous collision, events, sensors |
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
| `src/query.cpp` | Ray casts, overlap tests and shape casts against one shape |

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
| `Cylinder` | A radius and a half height along local z |
| `Hull` | A polyhedron from `computeHull` |
| `HalfSpace` | Everything below the local xy plane. Static bodies only |
| `Mesh` | A shared `TriangleMesh`. Static and kinematic bodies only |

See [decision 21](decisions/0021-cylinders-are-exact-in-3d.md) for how cylinders collide.

Boxes and prisms are cleaned when the shape is built, and `geometry` returns the cleaned
values. Negative half extents flip, and non-finite values are dropped or set to zero. A
prism polygon that repeats points, has collinear corners, is not convex or winds more than
once is replaced by its convex outline. A box or prism with no thickness in some direction
gets the same slab `computeHull` gives a flat cloud, so every polyhedron has volume and
every face has a normal.

`setLocalTransform` moves a shape on its body. Like `setGeometry`, it updates the bounds
and the mass, keeps the shape's contacts from being reused, and wakes the body.

A `TriangleMesh` merges vertices at the same position, drops vertices that are not finite,
and drops triangles whose widest corner has a sine under 1e-3. Its tree splits at the
median of the triangle centers along their widest axis, down to four triangles a leaf.
Each triangle keeps its normal and six edge flags, one per edge and side. An edge has a
normal of its own on a side when the neighbor folds away from that side by more than a
sine of 0.087, and on both sides when it is open or shared by more than two triangles.
Meshes have no mass and don't collide with each other or with half-spaces.

Every polyhedron keeps vertices, face loops, face planes, its bounds, and an edge list with
the two faces on each edge. Each edge also stores its length and its arc on the Gauss map,
the cross of its two face normals. The edge list is what the SAT edge queries walk.

`computeHull` is quickhull. It stops adding points once the hull reaches `maxVertices`,
always taking the point farthest outside the current hull, so a limited hull keeps the
corners that matter most. Faces whose normals are within a cosine of 0.9998 merge into one
polygon. A flat, collinear or single-point cloud is thickened into a thin slab first, so
every hull has volume. The slab is at least ten times the build tolerance, which grows
with distance from the origin, and a cloud is thickened at most once. Non-finite points are
skipped.

A face next to the visible region is also removed when the new face over their shared edge
would be a sliver or fold back past the tolerance. Points dropped while just above a face
are tested again against the finished hull, and a merged face whose outline is not convex
is grown again one triangle at a time. The result is convex within the tolerance and holds
every input point.

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
   - halfway through, fast pairs are collided again, and fast pairs without a planar body
     before every substep (see below)
   - integrate velocities
   - warm start joints, then contacts
   - solve joints with bias, then contacts, then again for impact contacts
   - integrate positions
   - relax: solve joints and contacts again without bias, then again for impact contacts
5. **Restitution.** One pass over points that were approaching faster than the threshold.
   Then the settle pass moves jointed bodies to close what joint error is left, and pairs
   that started or stopped touching are reported.
6. **Clean up.** Forces and torques clear, and kinematic bodies snap to their targets.
7. **Continuous collision.** Bodies that moved well past what their speculative contacts
   covered are swept and moved back. Then sensors find what overlaps them.
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

Ray casts, shape casts, overlap queries, the sweeps of continuous collision and waking
query the same tree. See [Queries](#queries).

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
- **Spheres** against polyhedra, spheres, cylinders and half-spaces have their own direct
  tests. So do polyhedra and cylinders against half-spaces.
- **Cylinders** against polyhedra, cylinders and mesh triangles test a set of feature axes:
  the polyhedron's face normals, the cylinder axes, and the crosses of edges and axes with
  each other, turned radial where they are parallel. Two GJK directions join them, one
  between the shapes and one between their cores, each cylinder shrunk by a quarter of its
  smaller dimension and rounded by the same amount. A GJK direction is used only when it
  beats every feature axis, which happens where a rim meets an edge or another rim. Against
  a mesh triangle GJK is skipped where it can't change the result: when the triangle is out
  of reach of the capsule around the cylinder's axis, or when a pair of points found by
  projecting back and forth between the two is already as close as the best axis. Each
  side then offers a feature along the normal: a cap as 24 rim points, the line of the
  side, a polyhedron face, or its one or two deepest vertices. A face clips the other
  feature, and two lines touch along their overlap, or at their closest points with the
  depth their ends reach along the normal. A spinning cylinder also offers the side lines
  next to the one touching, so the next one to land is already in the manifold. In an
  overlap deeper than twice the speculative distance the GJK directions are rough or
  missing, so the best axis then climbs by small turns while the gap grows, starting from
  itself and from each cylinder's cap normal and side direction. No direction gives more
  gap than the shapes really have, so the climb never overshoots.
- **Meshes** collide with each triangle whose bounds reach the other shape, in the mesh's
  frame. See below.

Every pair produces candidate points with a separation. Speculative points, up to the
margin apart, are kept so the solver can stop an approach before it turns into overlap.

A point sits at the middle of the gap, so both bodies push on the same spot. A point
farther apart than the speculative distance is anchored on each body's own surface
instead. A fast body only reaches such a point late in the step, and a midpoint anchor
would put its friction half the gap away from where it lands. Planar pairs keep the
midpoint.

### Manifold reduction

A manifold keeps at most four points. Each candidate gets a predicted separation: its
current separation plus the relative normal velocity at that point times `dt`. The point
predicted deepest goes in first. Among the points within a small band of it, the next goes
to the one farthest away, then the one giving the largest triangle, then the one farthest
outside that triangle. Free slots take the deepest of the rest.

Ranking on predicted separation rather than current separation matters for tumbling
bodies. The vertex about to land is kept even though it is not touching yet.

For a fast pair the spread is chosen among every point predicted to reach the surface
within the step, not only those near the deepest. A thin plate landing at speed otherwise
keeps the corners of the edge that hits first, pivots on that edge and flips up.

### Meshes

A triangle is a polyhedron with two faces, front and back, so the face queries, face
contacts and edge contacts above apply to it unchanged. Its edge query is its own: the arc
of a triangle edge on the Gauss map is the half circle from the front normal through the
edge's outward direction to the back normal. A sphere takes the closest point on the
triangle.

When a contact's normal is off the triangle's normal, the triangle feature it touches is
the one reaching farthest along that normal, an edge when its far corner is within slop.
The deepest contact point is a poor guide, since on an edge contact it lies on the edge
and rounding puts it on the face as often as not. If that edge or corner has no normal of
its own, neither on the side the body is on nor on the side the normal leans to, the
contact is made again as a face contact on the side the body is on. If any of its points
is deeper than SAT found the shapes to overlap and lies outside the triangle, the contact
belongs to a neighbor and is dropped. A sphere's contact on such a feature is dropped
unless the sphere's center lies over the face. See
[decision 19](decisions/0019-meshes-collide-one-triangle-at-a-time.md).

Contacts from triangles whose normals are within a cosine of 0.996 merge under the normal
of the deepest one and are reduced to four points. Each group becomes its own contact. In
a group the spread is picked among every point predicted to penetrate, not only those near
the deepest. On a curved patch the deepest points all sit at one end, and a manifold of
only those lets the body seesaw.

### The planar path

A pair tries the 2D path first when every dynamic body in it is planar and their plane
normals are parallel. Static and kinematic bodies fit any plane. `makeProfile` turns
each shape into a 2D profile in the plane: a polygon for a box or prism, a circle for a
cylinder, plus the depth range it covers along the normal. A half-space whose plane stands
on the body's plane, within the same tolerance, becomes a half-plane that covers every
depth. It fails for anything that does not extrude along the normal, and then the pair
uses the 3D path. Two half-planes always do.

A polygon meets a half-plane with the edge facing it most, a circle with one point. The
line is taken at the middle of the other profile's depth range, which also covers a
half-space tilted within the tolerance. See
[decision 23](decisions/0023-half-spaces-are-half-planes-on-the-planar-path.md).

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

Rolling resistance is an angular impulse against the relative spin of the two bodies,
limited by the mixed coefficient times the rolling radius times the total normal impulse.
The rolling radius is the mean height of the dynamic body's center of mass above the
contact points. The impulse is solved in one step along the current spin direction and
the total is clamped, so a planar body only takes impulse about the axes it can turn on.
It resists spin about the normal too. The mix is the maximum, and a coefficient of zero
skips it. See [decision 22](decisions/0022-rolling-resistance-follows-the-spin.md).

A direction the bodies cannot move in, such as the depth axis of a planar body, has an
effective mass of zero. Rounding can leave such a direction with a tiny nonzero response,
which would turn into a huge effective mass. Any response smaller than a millionth of the
two inverse masses is treated as zero. See
[decision 9](decisions/0009-locked-directions-have-no-effective-mass.md).

### Mid-step refresh

Contact points come from the start of the step. A body spinning fast carries those points
away from where the overlap really is, and the solver pushes on points that no longer mean
anything. Pairs listed as fast at the start of the step are collided again halfway through
the substeps. Fast pairs without a planar body are collided again before every substep
after the first, with the margin covering only the motion left in the step.

Before that, every active contact that isn't collided again is rebased. Its anchors become
the bodies' current lever arms, its stored separation absorbs the distance already moved,
and `dp`, `dq` and the base orientation of every body reset. New points then take impulses
from old points the same way as after the narrowphase, and the restitution data of the
nearest old point. A mesh pair is collided once, and each of its contacts takes the new
group whose normal is closest to its own.

See [decision 8](decisions/0008-collide-fast-pairs-again-mid-step.md) and
[decision 24](decisions/0024-fast-solid-pairs-collide-every-substep.md).

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
   tunneled and the body keeps its pose. Each advance uses the gap along the GJK normal
   when it is smaller than the GJK distance, since GJK can stop above the true distance.
2. If the core hits, the full shape is swept to find where it first touched. The body moves
   back to that time on its sweep.

Dynamic bodies that hit each other move back together, to the earliest time in the group,
so none is left inside another. A kinematic body keeps moving after the hit, so whatever
first hit one rides along with it from the time of impact to the end of the step. Against
a mesh, the full shape is swept against each
triangle under the path. The triangles are found in the mesh's frame at its poses at both
ends of the step.
Meshes are never fast themselves. Time is lost for a moved body. The next step carries on
with the velocity the solve gave it.

For a planar body, the inner extent ignores the faces facing along the plane normal. A
thin layer is wide in the plane, and that width is what keeps it from tunneling.

See [decision 7](decisions/0007-continuous-collision-past-the-speculative-margin.md).

## Queries

Queries read the world and never change it. Each takes a category mask, culls with the
bounds tree, tests half-spaces from their own list, and works against every shape kind. A
planar body's shapes are ordinary 3D shapes to a query. Results come in shape index order
unless noted.

- **Rays.** `castRay` returns the closest hit and `castRayAll` one hit per shape, nearest
  first and by shape index at equal fractions. Each shape kind has an exact ray test, and a
  mesh walks its triangle tree.
- **Points.** `overlapPoint` tests each shape analytically: the sphere's radius, the
  cylinder's rim and caps, a polyhedron's face planes, the half-space's plane. A mesh has no
  inside, so a point counts only on one of its triangles.
- **Boxes.** `overlapAabb` compares bounds only. A half-space counts when the box reaches
  past its plane.
- **Convex geometry.** `overlapGeometry` takes a sphere, box, cylinder, prism or hull at a
  pose and runs GJK against each candidate. Against a mesh, GJK runs on the triangles under
  the query's bounds until one touches.
- **Shape casts.** `castGeometry` moves the same geometry along a translation without
  turning and returns the first shape it touches. Each candidate is swept by conservative
  advancement on the same GJK distance, stepping by the larger of the gap along the GJK
  normal over the closing speed and the distance over the speed. A step that lands inside
  is bisected back. Against a mesh, the query's bounds slide down the triangle tree and
  each triangle under the path is swept. The point and normal come from the narrowphase,
  run just before the hit, with the point moved onto the hit shape. A cast that starts
  touching a shape reports fraction 0 and `initialOverlap`, with the deepest contact point
  and a normal pushing the geometry out, or a zero normal when there is none.

Shapes within 0.01 mm, scaled by `lengthUnit`, touch. Overlaps and casts use the same GJK
distance, query geometry first, so a cast's hit is found by an overlap query at the hit
pose. GJK can stop above the true distance on a cylinder by an amount that depends on
argument order, so near touch both orders are tried and the smaller kept. Casts don't use
`timeOfImpact`, which stops early on cylinder caps.

See [decision 25](decisions/0025-queries-share-one-touch-test.md).

## Sleeping and waking

A body's sleep timer runs while its speed, counting rotation at its farthest point, stays
under `sleepThreshold`. Islands are rebuilt every step with union-find over contacts with
points between awake dynamic bodies, speculative points included. The narrowphase wakes a
sleeping body near a moving one, so islands split only at touching contacts would keep
waking each other. An island sleeps when every body in it has been slow for `timeToSleep`.
A moving kinematic body keeps whatever touches it awake.

Waking floods through contacts with points and through joints, from the bodies that need
it. Each body's links are built when first needed after a step or a new joint:

- A dynamic body wakes its whole island.
- A static or kinematic body wakes the dynamic bodies touching it, without spreading past
  them, and every sleeping body whose bounds overlap its own. For a half-space the overlap
  test is against its plane.
- A sleeping body touched by a moving body wakes during the narrowphase.

Removing a shape or body wakes the dynamic bodies that touched it. Its contacts are marked
removed and dropped by the next step, so the cost is what the body touched. A body that
stops being dynamic loses its contacts with static and kinematic bodies.

See [decision 10](decisions/0010-moved-anchors-wake-what-they-overlap.md).

## Events and sensors

Each step collects contact, hit and sensor events, read with `contactEvents`, `hitEvents`
and `sensorEvents` until the next step. Contact and hit events are ordered by the shape
pair, lower index first, and sensor events by the sensor and then the visitor. Events from
calls made between steps come first, in the order of the calls. A copied world reports
exactly what the original does.

### Touches

Every contact keeps whether its pair was touching after the last step. All the contacts of
a mesh pair agree, so a pair begins and ends once however many facings it has.

After the solve, a pair is touching when the solver pushed on any of its points in some
substep, or a point overlaps. A speculative point only takes an impulse when its gap would
close within the substep, so shapes that pass close without meeting report nothing. Once
touching, a pair holds until nothing pushes on it and every point is further apart than the
linear slop. Resting contacts are springs that hover around zero overlap, and without that
margin a box at rest could end and begin again.

A pair also ends when its contacts disappear at the start of the step: its bounds moved
apart, the narrowphase found nothing, or a filter or motion type change dropped it. New
contacts of a pair take the state of its old ones, removed ones included, so a pair that
still touches after its filter changed reports nothing. Destroying a body or shape, and a
joint that stops its bodies colliding, end the touching pairs at once with the ids from
before the call. They are reported with the next step, and a shape slot reused before then
starts with no touches.

Sleeping pairs are not solved and keep the impulses of their last solve, so their state
holds and they report nothing, before or after they wake.

A hit is the fastest point of a pair that the solver pushed on, when it approached faster
than `hitThreshold` at the start of the step, 1 m/s by default. It carries that speed, the
point and the pair's normal from shape A to shape B. A pair gives at most one hit a step,
and only pairs solved in the step count.

Events are always on. Touches and hits come from the pass over contacts that already counts
the touching ones, and ends at the start of the step from the contacts the narrowphase
drops.

### Sensors

A shape with `ShapeDesc::sensor` set never collides. It has no contacts and no mass, doesn't
count toward its body's extent, and is skipped by continuous collision, waking, rays and
queries. Adding a sensor changes nothing about how bodies move.

After continuous collision, each sensor queries the tree with its bounds. A shape it reaches
is a visitor when it is not a sensor, is on another body, and passes the same category and
mask test as a collision. Shapes on static, kinematic, dynamic and sleeping bodies all
count, and a sensor works on a body of any type. Joints don't hide shapes from sensors.
Overlap is a GJK distance of zero or less between the shapes' 3D geometry, also for planar
bodies, and against a mesh, against each triangle the other shape's bounds reach. Two meshes,
a mesh and a half-space, or two half-spaces never overlap, as they never collide.

A pair where neither body is moving keeps its last result. Adding a shape, or moving,
reshaping or refiltering one, makes the next step test every pair again. The pairs found
are kept sorted, and the difference from the last step gives the begin and end events. A
sensor never wakes what it finds, and a sleeping body inside one stays inside.

With no sensors the pass returns at once. Each sensor costs a tree query a step and an
overlap test for each moving shape it reaches. A half-space sensor visits every shape.

See [decision 26](decisions/0026-events-follow-the-solve-and-sensors-test-overlap.md).

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

Strict float flags keep the compiler from fusing or reordering float math. Results still
depend on the platform's `sin`, `cos` and `atan2`, which turns, joints and outlines use, so
two platforms can differ.

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

- Shape casts don't turn, and take only convex geometry. A half-space or mesh as query
  geometry finds nothing.
- A cast that only grazes a shape, closing on it at almost no speed, can report a miss.
- Meshes can't be dynamic. Use hulls or several convex shapes for a moving concave body.
- A light chain holding a much heavier body stretches. Twelve links holding a ball about
  200 times heavier than each link rest 0.17 m longer in total.
- Contacts are springs, so a stack sags under its own weight. A three box tower rests
  with 7 mm of overlap between the bottom two boxes. More substeps or a higher
  `contactHertz` make it stiffer.
- Fast spinners that end up deep inside each other under a pile can stay 150 to 220 mm
  deep for 8 to 12 frames before they recover. This also happens without the mid-step
  refresh. Pushing deep pairs out harder fixes it but keeps some hull piles from ever
  sleeping.
- A thin body hitting the ground at 2000 m/s can still come back up at 5 to 6 m/s with
  restitution 0. Up to 500 m/s it stays under 0.6 m/s.
- Cylinder pairs overlapping by less than 4 cm keep the feature axis estimate, which in
  rare poses is off by up to 30 mm.
- A hull built thousands of units from the origin can have faces off by more than the
  tolerance, since the face planes come from absolute coordinates. Build hulls around the
  body and place them with the body's pose.
- Spheres and cylinders with a radius or height of zero, below zero or not finite are not
  cleaned up.
- `alphaPieces` is quadratic in the number of holes.
