# Devlog

Notes as I go, newest first. The longer reasoning behind some of this is in
[decisions](decisions/).

## 2026-10-30 12:54

Fast pairs are collided again before every substep, and that was a quarter to a third of
a step in Rain, Hulls, Meshes and Sweeper. Most of those collisions find what the last one
found. Pairs of polyhedra now keep the edge pair they last met on and test only that one,
as long as they turned less than 0.1 rad since the step started and the pair's arcs still
cross. The face queries still run every time, so a switch to a face is never missed.

Whole steps with sleep off, ms, before and after: Hulls 0.43 and 0.38, Sweeper 0.32 and
0.30. With sleep on, Hulls 0.32 and 0.27. Rain and Meshes are the same, since their time
goes to cylinders and triangles. Results change where fast polyhedra touch, and over 12
sets of 32 seeds the overlap counts move within noise.

Reusing faces too was faster still, but 8% of the reused features weren't the ones a
search would pick, and Hulls overlapped more. Dropping the turn limit did the same. See
[decision 34](decisions/0034-fast-pairs-keep-their-edge-pair-within-the-step.md).

## 2026-10-30 03:46

Collision is most of what's left between Twofold and Jolt in Rain, Hulls, Sweeper and
Meshes, and nearly all of it is pairs that touch. The edge query's Gauss map tests and the
face queries' vertex scans now take four candidates at a time, and only the pairs that pass
are finished one by one, in the old order. The GJK check after a pruned edge query stops
as soon as its bounds settle the answer. Two loops recomputed the best dot product on
every step, and the edge box test now marks each vertex once. None of it changes a result:
every scene keeps its hash, with and without SIMD.

Whole steps with sleep off, ms, before and after: Hulls 0.529 and 0.470, Rain 0.687 and
0.653, Sweeper 0.328 and 0.308, Meshes 1.414 and 1.384. With sleep on, Hulls goes from
0.437 to 0.326.

Tried and dropped: warm-starting GJK (Rain still averages ten iterations), a looser GJK
tolerance (it moves normals, like the certificate rejected before), and caching separating
axes for pairs with no contact (they only cost 0.03 to 0.06 ms a frame). See
[decision 33](decisions/0033-sat-tests-run-four-candidates-at-a-time.md).

## 2026-10-30 01:47

Contacts are solved four at a time now, with SSE2 or NEON. They are colored first, so no
two contacts of one color move the same body, and each color is cut into batches of four.
Each lane runs the same float operations as the single contact solve, so a build with
`TWOFOLD_SIMD=OFF` gives the same bits. Checked on MSVC with and without AVX2, GCC with
`-march=native` and Clang with FMA enabled. The NEON path builds for aarch64 macOS but
hasn't run yet. See [decision 32](decisions/0032-contacts-are-solved-four-at-a-time-by-color.md).

The first version was slower than the old solve everywhere but Pyramid. Two things cost
the most: copying every contact's axes into its batch after every mid-step refresh, and a
gather that called its source once per component instead of once per lane. Batches now
build their own axes, and warm starting and restitution run on them too.

Whole steps with sleep off, ms, before and after: Pyramid 0.480 and 0.329, Hulls 0.585
and 0.542, Meshes 1.705 and 1.601, Dominoes 0.133 and 0.128. Rain is the same, since its
time goes to collision. Solving in color order changes results: stacks rest about 0.3 mm
deeper, and Dominoes has more frames past 5 mm.

## 2026-10-30 00:02

The repository is public now, and CONTRIBUTING.md says what a pull request needs: a test,
numbers for any claim, passing CI and updated docs.

## 2026-10-29 19:40

The README now says which platforms AnvilFx tests by hand and which only CI checks.

## 2026-10-29 16:12

CI failed on macOS: a degenerate prism in the drop test was still in the air after two
seconds. The macOS build fuses multiply-adds in the tests, so the prism started in a
slightly different pose. The same thing happens on Windows with other poses. A 1 cm by
1 cm by 1 m box rod dropped flat from 1 m in 200 random orientations came back above its
drop height 32 times, the worst at 1.8e8 m.

The world inverse inertia was rebuilt after every position integration and the contact
responses with it, but the effective masses were built once per step. A rod has about
5000 times more inertia across it than along it, so a small turn made the two disagree and
every impulse added energy. The inertia now stays fixed for the step, and the per substep
response update is gone. See
[decision 31](decisions/0031-world-inertia-is-fixed-for-the-step.md).

200 drops of each, above the drop height / fastest spin:

| Body | Before | After |
| --- | --- | --- |
| 1 cm box rod | 32, 2e10 rad/s | 0, 400 rad/s |
| 1 cm cylinder rod | 13, 5.6e15 rad/s | 0, 149 rad/s |
| Degenerate prism | 36, 4.6e9 rad/s | 1, 363 rad/s |

What's left is a rod landing flat at 3 m/s, which can still bounce to 0.83 m with 4
substeps. 8 substeps settle it. Capping the turn per substep at 45 degrees takes the
bounce down to 0.49 m.

The older note that Dominoes needs the response update no longer holds. Over 12 seed sets,
frames past 5 mm:

| Scene | Before | After |
| --- | --- | --- |
| Dominoes | 12969 | 6676 |
| Sweeper | 9972 | 9899 |
| Hulls | 12000 | 11760 |
| Flat layers | 18888 | 18890 |
| Rain, 4 sets | 72944 | 72805 |

Frames past 20 and 50 mm are the same within noise, except Rain, where they fall from
15179 to 14652 and 1769 to 1649. Rain's deepest overlap over the last half second is worse
in 4 of 6 seed sets, 21 mm against 14 averaged over the sets. Tower, Pyramid, Slopes,
Sleeping and Cutouts didn't change.

Tests:

- A thin rod drop test, 40 orientations. It fails before this change.
- The hull well-formed check took face normals from cross products of absolute positions,
  which lose precision far from the origin. They are now relative to a face vertex.
- The far from origin hull test compared the height of a ball that is still rolling after
  two seconds. It now checks that both heights lie between the ball's inner and outer
  radius.
- test_outline.cpp was missing <algorithm>.

## 2026-10-29 04:00

The README now gives examples of what AnvilFx directs: the 2D path, the bounding volume
hierarchy and the substepped soft contact solver.

## 2026-10-29 02:43

Cutouts is layers only, and it took 0.44 ms a step, about three times Jolt with its bodies
locked to the plane. That was the 3D engine with axes locked, which is not the point. Layer
pairs now collide and solve in 2D, and the bounds tree holds bodies instead of shapes.
Mixed pairs, a layer against a free body, are untouched. See
[decision 29](decisions/0029-planar-pairs-are-solved-in-2d.md) and
[decision 30](decisions/0030-bounds-tree-leaves-are-bodies.md).

- Profiles are made once per shape in its body's frame and kept. Each body is placed in the
  plane once per round of collisions, so a pair costs one turn and one shift before the
  polygon test. A box test drops pairs out of reach first, and the separation loop stops
  once an edge separates them.
- Contacts of layers in one plane solve with scalar arms and an inverse inertia about the
  normal. Contacts holding nothing skip warm starting and their friction rows.
- A planar body no longer rebuilds its world inertia every substep.
- One leaf per body. A cutout of several pieces kept a pair for every two pieces whose
  leaves met. Bodies of more than 64 shapes keep a tree of their own, or one static body of
  9,800 boxes went from 0.19 to 0.35 ms a step. With it, that is back to 0.18.
- Smaller: box and prism bounds from their outline, the planar inner extent kept per
  shape, step lists kept between steps, and min and max that select instead of branching.

Pinned minimum over 8 seeds, ms a step:

| Scene | Before | After |
| --- | --- | --- |
| Cutouts | 0.442 | 0.125 |
| Flat layers | 0.067 | 0.044 |
| Rain | 0.606 | 0.582 |
| Sweeper | 0.319 | 0.297 |
| Hulls | 0.432 | 0.415 |
| Meshes | 1.47 | 1.40 |
| Joints | 0.192 | 0.184 |

Tower, Pyramid, Dominoes, Slopes and Sleeping are the same. In the Jolt harness Cutouts is
now 0.13 ms a frame to Jolt's 0.145, with a worst overlap of 17 mm to Jolt's 141.

Only Flat layers, Cutouts, Joints and Meshes change. The pure 3D scenes give the same
results bit for bit, with fewer candidates tested. Against 2771959:

- Cutouts, 48 sets: frame mean 4.333 to 4.321 mm, frames past 10 mm 14051 to 13968, past
  20 mm 415 to 390, past 50 mm 39 to 18, contact overlap 0.2690 to 0.2694.
- Flat layers, 24 sets: frame mean 2.533 to 2.522, past 10 mm 10233 to 10189, past 20 mm
  1585 to 1619, none past 50 mm, contact overlap 0.3307 to 0.3292.
- Joints, 12 sets: the same.
- Meshes, 48 sets: frame mean 4.018 to 3.991, past 10 mm 113657 to 110258, past 20 mm
  12770 to 12890, past 50 mm 300 to 297, contact overlap 0.1993 to 0.1987.

Meshes only changes by a different path through chaos, and single sets swing by 480 frames
past 20 mm either way. Its worst set has two cylinders 20 mm into each other for a while.
Their kept contact reads no overlap while a fresh collision finds 20 mm. That is the
contact reuse in 3D, reached by a different trajectory, not the 2D work.

Dropping speculative points that weren't closing looked like a cheap win and wasn't. Over
twelve Cutouts sets frames past 20 mm went from 107 to 649 and past 50 mm from 1 to 78, and
a test of a near miss failed.

Builds with GCC and Clang now. GCC stopped at a missing `<cstdint>` in the hull builder,
which MSVC's headers pulled in on their own. Clang compiled everything, but four shape cast
tests failed from -O1 up: `Vec3::operator[]` read `(&x)[i]`, which indexes past `x` and is
undefined, and casts against a mesh missed. It picks the member by index now. MSVC gives
the same results bit for bit and the same timings.

On macOS the viewer keeps a display scale of one, since the system already scales the
window. A CI workflow builds and tests on Linux with GCC and Clang, on macOS and on
Windows.

## 2026-10-29 00:50

Meshes was the slowest bench scene and the furthest behind Jolt, about six times its time.
It is down about a fifth, 1.74 to 1.35 ms a step pinned. Nothing else moves.

- A fast mesh pair keeps the triangles it ruled out at the start of the step with their
  distances, and the collisions within the step skip those still out of reach.
- Against a cylinder, where the capsule around its axis is closest to the triangle between
  the caps, the direction between the closest points is exact and GJK isn't needed.
- Near a cap GJK still ran about 75 times a step, and changed the answer in almost all of
  them. Projecting back and forth between triangle and cylinder now runs longer, and once
  the direction between the pair separates the shapes to within 0.05 mm of the pair's
  distance it stands in for GJK. That leaves about 24 calls, which start from that
  direction. A 0.5 mm bound was faster but let the direction drift well off, since the gap
  is flat near the best one, and overlaps grew (frames past 20 mm 1123 to 1812 over six
  sets). Pairs closer than the slop always go to GJK.
- A pruned polyhedron triangle skips the GJK check when its deepest vertex already reaches.
- Smaller: cap points from a table, both signs of an axis in one pass, and the solver no
  longer turns anchors on bodies that haven't turned or pushes anchored bodies.

Only Meshes changes. Single seed sets swing a lot here: starting GJK elsewhere, which moves
results only within its tolerance, changed one set's frames past 20 mm from 150 to 408. So
this was measured over twelve sets of 32. Against 65373dd: frame mean 3.946 to 3.899,
frames past 10 mm 26080 to 25334, past 20 mm 2821 to 3102, past 50 mm 59 both, mean
contact overlap 0.1997 to 0.1972. The same within noise.

## 2026-10-28 20:26

Continuous collision stopped short of flat faces, a few millimeters on a box or a cylinder
cap and 15 mm on 5 m cylinders. Near a face GJK ends on a long thin triangle, and the float
barycentric weights tilted its normal by about a hundredth of a radian. The distance stayed
right but the gap along the normal didn't, and that gap sets each advance. The weights are
now found in double, and swept shapes stop within about a millimeter of the slop.
Rain and Meshes change. Over six seed sets Meshes overlaps less (frames past 10 mm from
15099 to 11529), Rain is the same within noise and Hulls is unchanged. Both scenes run a
little faster.

## 2026-10-28 19:06

Swing-twist joints, for ragdolls. The swing stays in an elliptical cone and the twist in a
range. The first limit rows turned about the limb's axis and the swing axis, and a limb
held at its twist limit tilted out of its cone. The rows follow the exact rates now, the
spin that changes only the twist and the one that changes the room left in the cone. A
ragdoll dropped from 10 m, or six piled up, stays together and holds its limits within
0.002 rad at rest. See [decision 27](decisions/0027-swing-twist-rows-follow-the-exact-rates.md).

Friction on hinges, sliders and swing-twist joints, alongside motors and limits. A hinge
with only a motor at speed zero did the same, but then it couldn't have a motor too.

`jointInfo` gives the force and torque a joint applied in the last step. Reading the last
substep alone was 2% high for a hanging body, so it is the sum over the step. Joints break
past `breakForce` or `breakTorque`, before the settle pass, and `brokenJoints` lists them.

Spheres and cylinders with a radius or height of zero, below zero or not finite are
cleaned like boxes now.

A hull built 10000 units from its body had normals off by 0.012, and at 100000 units
quickhull lost points. Hulls are built about their own center. Bench scenes with hulls and
prisms change hashes, from normals moving by 5e-6 at most. See
[decision 28](decisions/0028-hulls-are-built-about-their-own-center.md).

## 2026-10-28 18:43

Contact and sensor events. Before this the only way to know two shapes met was to read
`contacts()` every step and diff it, which has a contact per facing for mesh pairs and
holds speculative contacts that never touch.

A pair touches when the solver pushed on one of its points in some substep, or a point
overlaps. Having points at all would begin a touch for bodies that only pass close, and
overlap alone flickers on a resting box, whose contacts hover around zero. It holds until
nothing pushes and every point is a slop apart. Each contact keeps the flag and all
contacts of a pair agree, so a mesh pair begins and ends once. Hits are the fastest pushed
point of a pair over `hitThreshold`, one per pair a step.

Destroying a body or shape ends its touches right away, with the old ids, and they come
first in the next step's list. A slot reused before that step can't inherit a touch.
Filter and motion type changes only mark contacts removed, so the step decides. A pair that
still touches under its new filter reports nothing. Sleeping pairs keep their impulses and
stay quiet.

Sensors are a shape flag and stay out of the contact list. Sensor contacts would need
skipping in waking, islands, the mid-step refresh and continuous collision. A pass after
continuous collision queries the tree with each sensor's bounds and runs GJK on what it
finds, and pairs where neither body moved keep their result. Sensors have no mass, wake
nothing, and rays and queries go through them. They see static, kinematic and sleeping
bodies, follow the collision filter, and ignore other sensors and their own body.

With no sensors the hashes are the same as before and the bench scenes time the same
within noise. See [decision 26](decisions/0026-events-follow-the-solve-and-sensors-test-overlap.md).

## 2026-10-28 18:39

Added overlap queries for points, boxes and convex geometry, shape casts, and a ray cast
that returns every hit. All of them take a mask, cull with the bounds tree and return
shapes in index order.

The first version of the casts used `timeOfImpact`. Against cylinder caps it stopped up to
4.4 mm short, because it accepts a hit once the gap along the GJK normal is under its
target, and that normal tilts on a cap. Changing it would change stepping, so the casts
have their own advancement on the plain GJK distance, bisecting back when a step lands
inside.

The random test that checks casts against overlap queries still failed after that. GJK on
a cylinder sometimes runs out of iterations and stops a little above the distance, by a
different amount depending on which shape goes first. The cast and the overlap test now
pass the query geometry first, and near touch with a cylinder both orders run and the
smaller distance wins. Hit normals across a gap of a few micrometers were off by up to 7
degrees on cylinder rims, so the normal now comes from the narrowphase half a millimeter
before the hit.

Scene hashes are unchanged.

## 2026-10-28 04:28

Won back most of what colliding every substep cost, with results bit for bit the same.

Most of Meshes went to GJK between cylinders and terrain triangles. A triangle out of reach
of the capsule around the cylinder's axis is dropped before it now, and GJK is skipped when
two points found by projecting back and forth are already as close as the best feature
axis. That halves the calls. Contacts about to be collided again were rebased first for
nothing, and MSVC turns `std::fmax` and `std::fmin` into library calls, so `maxf` and
`minf` do the same inline. A body's extent is kept on the body, and the rebuilt contact
list keeps its storage between steps.

Over the bench scenes a step took 2.75 ms before the plate fix, 3.61 ms after it, and
3.01 ms now. Every scene but Meshes is at or under where it was before the fix. Meshes is
at 1.44 ms against 1.14 ms.

## 2026-10-28 03:01

A thin plate hitting the ground at 100 m/s or more bounced back up with restitution 0, at
up to 34 m/s. Reduction kept the corners of the edge that hit first, so the plate pivoted
on it, and friction on speculative points meters apart acted half the gap away from where
the plate landed.

Fixing both stopped the bounce but made Cutouts, Flat layers, Meshes and Rain worse.
Planar pairs were most of that, so they keep midpoint anchors now. Fast pairs without a
planar body are also collided again before every substep instead of only halfway. That
costs up to 46% more in Meshes, but overlap went down in Rain, Hulls, Meshes, Sweeper and
Joints. The plate stays under 0.6 m/s up to 500 m/s. See
[decision 24](decisions/0024-fast-solid-pairs-collide-every-substep.md).

## 2026-10-28 02:09

A dynamic body resting on the ground and then made static kept its contact with the ground
forever, two static bodies touching. A body that stops being dynamic now drops its
contacts with static and kinematic bodies, and the step won't carry such a contact either.

Removing contacts rebuilt a list of every contact and joint for each body to wake from, on
every call. Destroying 8000 resting pairs one at a time took 1.1 s. Removed contacts are
now only marked and dropped by the next step, and waking walks links built once per step.
The same takes 2.5 ms.

`setLocalTransform`, to move a shape on its body. It does what `setGeometry` does.

A ball knocked into a fast kinematic paddle rewound to the time of impact while the paddle
kept its end pose, so it could end up 0.3 m inside it. A body whose first hit is with a
kinematic body now rides along with it to the end of the step.

Time of impact advanced by the GJK distance, which is an upper bound when GJK stops early.
Nearly parallel cylinders stop early often, and some long sweeps ended up to 4 mm inside.
It advances by the gap along the GJK normal now when that is smaller.

Ray casts said the lowest shape index wins a tie, but a later shape at the same fraction
replaced the earlier one. Fixed.

Planar bodies against a half-space always took the 3D path. That is most pairs in a 2D
scene. A half-space standing on the body's plane is now a half-plane on the planar path.
See [decision 23](decisions/0023-half-spaces-are-half-planes-on-the-planar-path.md). The
planar path also keeps small polygons inline instead of allocating two vectors per pair.

`alphaPieces` was quadratic on comb-like outlines: 1000 teeth took 1.3 s, 4000 took 26 s.
The farthest point search in Douglas-Peucker walks a box tree now, and ear clipping only
tests vertices in nearby grid cells. 25 ms and 110 ms. Near ties between points on long,
nearly flat runs are now settled exactly, which changes the result on a few images by a
vertex or two.

Thin bodies hitting the ground at 100 m/s and more still come back up with restitution 0.
Reduction drops the speculative points about to land, and contacts anchored at the
midpoint give a fast spinner the wrong lever arms. Keeping those points and anchoring on
each body's surface stops it up to 500 m/s, but leaves resting cylinders a little deeper
and breaks the spinning cylinder test, so it waits.

## 2026-10-28 02:08

Went looking for ways to break things with bad input, and found a few.

A box with a zero half extent, or a prism whose polygon repeats its first point, has faces
with no normal. Contacts against them came out with a zero normal and zero separation, and
a box dropped on one fell straight through. An empty prism read past the end of its
polygon, and a box with one negative half extent was inside out, with negative volume.
Boxes and prisms are now cleaned when the shape is built. Negative half extents flip,
repeated and collinear points go, a polygon that isn't convex becomes its convex outline,
and anything with no thickness gets the same slab as a flat hull cloud.

Quickhull could return hulls that weren't convex. About one random sphere cloud in a
hundred had a vertex up to 8% of the hull size in front of a face, and squashed clouds
were much worse. An eye point just over the tolerance above one face but under it above
the next left that neighbor on the horizon, and the new face over the shared edge folded
back. Such neighbors now go too. Flat clouds also left points outside the hull, because
points dropped as nearly on a face were never looked at again, and merged faces could have
reflex corners. A single point far from the origin recursed forever, since the slab didn't
grow with the tolerance. Apart from hulls far from the origin, the fuzz cases all come out
convex now, with every point inside.

Meshes had ghost contacts on flat interior edges. Whether a contact was on an edge or the
face came from the triangle feature closest to the deepest point, and for an edge contact
that point lies right on the edge, so rounding often called it the face and the tilted
normal stayed. On a flat grid, up to two thirds of box and cylinder poses got a tilted
manifold at 10 mm overlap. The feature now comes from the triangle's support along the
normal, and there are none.

That made landing on hilly terrain sink twice as deep. The ghosts had been propping boxes
up in the valleys, and the face contact that replaced them only kept points at the deep
end, so the box rocked. A mesh group now picks its spread from every point predicted to
penetrate. On the test terrain the deepest point per frame while things land averaged
10 mm, and now 3.5 mm.

Cylinder side lines tilted against the normal reported the gap where the two segments pass
closest, which is about zero when they cross, even 86 mm into a box. They report how far
their ends reach now. Deep cylinder overlaps went wrong the other way, up to 0.44 m too
deep, since past the core shrink there are no GJK directions. In deep overlaps the best
axis now climbs while the gap grows.

Meshes drop vertices that aren't finite, which made the sorts undefined, and slivers whose
widest corner has a sine under 1e-3.

## 2026-10-27 22:13

Exact cylinders. The 24-sided prisms bobbed and lost speed on every side, which shows on
anything that rolls. Cylinder pairs now test feature axes from faces, cylinder axes and
edge and axis crosses, with GJK directions where a rim meets an edge. A cap contact is 24
rim points.

The first version had contacts 80 to 240 mm deep between cylinders that were really 2 to
7 mm into each other. GJK was calling a nearly flat tetrahedron an overlap. Curved supports
make those flat simplices common. Faces too flat to tell their inside from rounding now
count as facing the origin.

Rolling resistance, so the exact shapes can stop. One angular impulse along the relative
spin, clamped as a vector, so a planar disc only takes impulse about the axis it turns on.

## 2026-10-27 19:19

Triangle meshes, for static and kinematic bodies. `TriangleMesh` takes vertices and
indices, merges vertices at the same position, drops degenerate triangles and builds a
median split tree over the rest. Shapes share it through a `shared_ptr`, so one terrain
can back any number of bodies. Both sides of every triangle collide, since imported models
are rarely closed or consistently wound.

Each triangle goes through the existing convex code as a two-faced polyhedron. Face
queries and clipping needed nothing new. The edge query did: a triangle edge's arc on the
Gauss map is a half circle from the front normal round to the back, not the short arc
between two faces, so it has its own version.

Seams are the usual problem. Left alone, a box sliding across a flat grid hits the edges
of the next triangle and hops. Each edge carries a flag per side, set when the
neighbor folds away from that side, and a contact on a flagged-off edge or corner is made
again against the face. Boxes slide and balls roll straight across.

A box in a V crease needs two normals at once, so a pair can now hold more than one
contact. Only mesh pairs ever do, and every older scene steps bit for bit as before.

The Meshes scene found a bad one. The deepest penetration plot sat at 300 to 800 mm with
nothing visibly wrong. A body resting near a ridge is behind the plane of the triangle on
the far side, that triangle's ridge edge is flagged off from behind, and the face contact
that replaced it pushed straight down through the body. The side also came from the sign
of a nearly sideways normal, which can point either way. The side now comes from the face query, and a
face contact deeper than SAT says the shapes overlap is dropped, since the neighbor
already covers it. Deepest contact on a hilly mesh is now a few centimeters while things
land.

Rays use Möller-Trumbore down the tree. Continuous collision sweeps a fast body against
each triangle near its path.

Scene: 131 bodies on a 4608 triangle terrain and a rocking kinematic dish, about 1 ms a
step averaged over 20 seconds.

## 2026-10-27 18:41

Joints: point, fixed, hinge, slider and distance. Hinges and sliders take limits and a
motor, and a motor at speed zero works as friction. A distance joint is a rod, a rope with
only a maximum, or a damped spring. Leaving body A unset joins B to the world, like the
ground body in Box2D but without having to make one. Joined bodies skip collision with
each other, joints tie islands together for sleeping, and destroying a body takes its
joints with it.

They solve like the contacts, soft in the biased pass and rigid in the relax pass, at
60 Hz. A 12 link chain with a ball about 200 times heavier than a link stretched more than
I liked, so two settle iterations after restitution move the bodies to close what's left.
Dropping to 30 Hz made the chain four times worse. Put the same chain in Jolt at its
default 10 velocity and 2 position steps to see where that stands:

| Case | Total stretch | Jolt |
| --- | --- | --- |
| Swinging | 0.20 m | 2.22 m |
| Hanging still | 0.17 m | 0.30 m |
| Hanging still, ball 21 times a link | 0.023 m | 0.046 m |

Good enough to leave without a tree solver.

The test scene found a bounce. A door swinging into its limit came back off it at 14% of
the speed it hit with, a slider with its anchor off center at 38%. The limit row only
turned the door about its own center, the point row then had to pull the center back
round the hinge, and between them the warm started point impulse came back as spin the
wrong way. Taking warm starting off the limits or turning off the settle pass did
nothing. Jolt's order, limits last, was worse. Softer joints helped but brought the
stretch back.

What fixed it was making the hinge's limit and motor rows carry the point impulse that a
turn about the axis needs to keep the anchors together, so the row's mass is that of
turning about the anchor. Same for the slider's rows with the turn that keeps it from
rotating. Both now stop dead at the limit.

Added a Joints scene to the viewer with every type: chain, weld, door, pendulums with and
without friction, a rail, a motor lift, rod, rope, spring, a planar chain on a peg, a
plank bridge and a paddle wheel. Joints draw in their own color with their axes, and a
red line wherever two anchors that should meet have come apart. Every older scene steps
bit for bit as before.

## 2026-10-27 16:16

Broadphase. Sort and sweep sorted every shape every step, sleeping or not, and Cutouts spent
a quarter of its frame on it. Replaced it with a dynamic AABB tree like Box2D's: fat leaves
with a 100 mm margin, surface area cost on insert, rotations to keep it balanced.
Half-spaces stay out of the tree since their bounds are infinite.

First version had every awake shape query the tree every step. Faster on most scenes but
Cutouts got slower, 0.126 to 0.16 ms. A traced cutout is dozens of small pieces, and each
one queried even when nothing near it changed. Now the world keeps the pairs whose leaves
overlap and only shapes whose leaf moved look for new ones.

Kept the candidate test exactly as it was, so every scene and every ray hit comes out bit
for bit the same as before. Rays sort what the tree finds by shape index so ties still go
to the lowest. One thing the old code got for free: pairs where neither body is dynamic
aren't kept, so a static body made dynamic has to look for its pairs again. Without that
it sank into what it sat on until its leaf moved.

| Scene | Broadphase before ms | After ms |
| --- | --- | --- |
| Pyramid | 0.010 | 0.002 |
| Rain | 0.034 | 0.024 |
| Hulls | 0.016 | 0.007 |
| Cutouts | 0.127 | 0.097 |
| 9,800 shapes, 40 awake | 2.46 | 0.13 |

Continuous collision in Cutouts also dropped from 0.041 to 0.015 ms, since the sweeps only
look at what's near the path now. The big scene went from 5.2 to 2.7 ms a step overall.

## 2026-10-27 15:58

Speed. Rain was taking 1.2 ms a frame against Jolt's 0.14. Profiled it first: about half in
the narrowphase, half in the solver, the broadphase barely shows.

Narrowphase:

- Collision runs in A's local frame now. Only B gets transformed, and with a matrix
  instead of a quaternion per vertex.
- Face planes and each edge's arc on the Gauss map are stored once per shape instead of
  rebuilt for every query.
- The face query starts with the face pointing at the other shape and tests B's deepest
  vertex first, so most faces are ruled out with one dot product.
- Edge queries only pair edges that reach into the overlap of the two bounds. Two rocks
  went from 96 x 96 pairs to a handful.

That last one was wrong at first. Ran old and new collision side by side on the same
random pairs and about 1 in 1600 came back touching when the old code said apart, some
27 mm deep. When the shapes really are apart, the edge pair that proves it can sit outside
the overlap of the bounds. A GJK distance check now runs whenever edges were skipped, and
the pair is dropped when the shapes are further apart than the margin. 400,000 pairs agree
after that, except 6 at the margin that the old code kept as contacts and the new code
doesn't, since they're actually further apart.

Solver: the lever arms, the velocity change per unit impulse and the effective mass of
every point direction are built once per step. Turned bodies need their responses updated
after each position integration. Without that, Dominoes frames past 5 mm doubled.

Then contact reuse, like Jolt's body pair cache. A pair that moved less than 0.5 mm and
turned less than 2 degrees relative to itself since it was last collided keeps its points
and only measures the separations again. 1 mm made the dominoes worse, 926 frames past
5 mm against 500. 0.5 mm leaves them as they were.

| Scene | Before ms | After ms | Jolt ms | Frames past 50 mm before / after |
| --- | --- | --- | --- | --- |
| Rain | 1.76 | 0.65 | 0.18 | 441 / 503 |
| Hulls | 1.56 | 0.39 | 0.09 | 35 / 20 |
| Pyramid | 0.115 | 0.070 | 0.014 | 0 / 0 |
| Dominoes | 0.143 | 0.101 | 0.041 | 0 / 0 |
| Cutouts | 0.74 | 0.65 | 0.14 | 0 / 0 |

Rain past 100 mm went from 29 to 24, so the tail is about where it was. The rest of the
gap to Jolt is the edge query on big shapes and the solver. Cutouts barely moved since
its time is in the broadphase and the planar path, which this didn't touch.

## 2026-10-27 14:48

A heavy box landing on a light sphere pushed the sphere 70 mm into the floor. The box hits
the sphere, but the sphere's contact with the floor only gets its turn after the sphere has
already been shoved down. One pass per substep doesn't carry an impact through a stack. Any
body in a pair closing faster than the restitution threshold now gets all its contacts
solved twice per pass. Down to 47 mm, and nothing at rest pays for it. Rain is 16% slower.

The Tower preset sits with 7 mm of overlap between the bottom two boxes, even asleep.
That's just the springs. 30 Hz holding up the weight above sags 7.2 mm. Box2D does the
same. 8 substeps or a higher `contactHertz` stiffen it. Leaving the defaults.

While watching the towers: three box stacks never fall asleep. They rock and twist at about
9 Hz forever. Two things feeding it.

- Warm starting gave a new point the full impulse of the first old point in reach. When the
  manifold changes shape, two new points can take the same old point and that side gets
  pushed twice as hard. Now each new point takes the nearest old point and they split it.
- Friction points are always solved in the same order, so the first absorbs the error and
  the last overshoots, every substep. The relax pass now goes through them backwards.

Tried reversing the whole relax pass first. Stops the rocking too, but the overlap tail in
Rain doubled. Reversing only friction keeps the tail where it was. 108 seeds:

| | Rain, frames past 50 / 100 / 150 mm | Hulls awake at 900 frames |
| --- | --- | --- |
| Before | 519 / 29 / 1 | 38 |
| Nearest point, split | 527 / 32 / 5 | 16 |
| Plus whole relax reversed | 531 / 62 / 18 | 21 |
| Plus relax friction reversed | 441 / 29 / 0 | 27 |

Towers, Cutouts and Hulls now all end asleep. Flat layers goes from 16 bodies awake after a
minute to 5.

Also ran the same scenes through Jolt to see where things stand. Overlap is lower here in
every scene but Pyramid, where Jolt with linear casting and 2 collision steps rests at
1 mm against 2.2. Jolt is still 4 to 9 times faster per frame.

## 2026-10-27 12:36

Found a hole in waking. A box sleeping on the ground, then a static box teleported with
`setTransform` so the two overlap by 0.9. After 60 steps the box was still asleep and
still overlapping. Waking follows existing contacts, and these two had never touched. The
broadphase also skips any pair where neither body is moving, so the overlap was never
even looked at.

Two ways to fix it. Flag moved anchors so the next broadphase pairs them with sleeping
bodies, or wake every sleeping body whose bounds overlap the anchor's new bounds right
away, which is what Box2D does. Went with the second. The cost lands on the call that
moved the body instead of the step, and it covers `setGeometry` and adding shapes for
free. Half-spaces have infinite bounds, so they get a plane test instead. Otherwise
nudging the ground would wake the whole world.

While in there, noticed that removing a body woke everything else on its floor. The
floor was in the list of touched bodies, and waking from a static body spreads to all of
its contacts since the change on the 26th. Removal now only wakes dynamic bodies.

## 2026-10-27 12:17

Fast spinners. A cylinder spinning at 25 rad/s can sink 50 mm into another before the
contact catches up, because its contact points are computed once at the start of the step
and the spin carries the real overlap away from them. Pairs that move more than the
speculative distance in a step are now collided again halfway through the substeps.

First version made rolling discs slow down for no reason. The new anchors were measured
in the wrong frame and acted like rolling resistance. Fixed by rebasing every contact to
the bodies' current lever arms before re-colliding. The second problem was bounce:
restitution needs the approach speed from the start of the step, and the new points
didn't have it, so a full-restitution ball lost about 15% of its height. New points now
copy it from the nearest old point.

36 seeds each:

| Scene | Frames past 50 / 100 mm, before | After |
| --- | --- | --- |
| Rain | 907 / 110 | 589 / 53 |
| Hulls | 298 / 20 | 129 / 0 |

Costs 13% on Rain and 20% on Hulls. Pyramid, Tower and Dominoes don't change at all.

Also tried pushing deep pairs out harder. It clears the worst overlaps faster but some
hull piles then never go to sleep, so that's out. There is still a case where two spinning
cylinders get trapped under a pile, 150 to 220 mm deep for 8 to 12 frames. It happens
without the refresh too. Not fixed.

Chasing Cutouts pieces that leave the bowl and glide 100+ m across the ground.
A planar box sliding at 3 m/s on a half-space went 27 m. On a box floor, 0.84 m. Planar
bodies on a half-space go through the 3D path, and the half-space normal has rounding in
z from the rotation that made it. The friction tangent then has a tiny free component,
so its effective mass is enormous, and that one direction eats the whole friction cone.
Directions with a response under a millionth of the inverse masses now get no effective
mass. Same box now stops at 0.83 m.

One test broke from that: a traced ring with a bar through it no longer settled in 480
frames. With real friction it rolls back and forth between the bar ends, same as it does
on a box floor, and sleeps by frame 671. Gave it 900 frames. In Cutouts, bodies still
awake after a minute went from 34 (gliding) to 6 (swinging on a peg).

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
