# Twofold

Rigid body physics for scenes that mix 2D and 3D bodies in one world.

A body locked to a plane is still a 3D body, so it collides with everything else in the
world. Two planar bodies that share a plane collide and solve through a dedicated 2D
path: their profiles are tested as polygons and circles, their contacts are solved in the
plane, and contacts never push them apart along the plane normal.

## Who does what

The code, tests and documentation are written by Claude, Anthropic's AI model. AnvilFx
does not write code. AnvilFx directs the approaches and techniques the engine uses and
makes sure they are modern and fast, rather than leaving those choices to Claude.

Every claim Claude makes about a change, such as that it works or that it is faster,
is tested by AnvilFx before it ships.

## Features

- Sphere, box, cylinder, extruded convex polygon and half-space shapes
- Convex hulls of point clouds, with a vertex limit
- Outlines traced from alpha images and cut into convex pieces, for exact 2D shapes with
  holes and concave edges
- Two-sided triangle meshes for static and kinematic bodies, with no catching on the
  edges between triangles
- Bodies made of several shapes, with mass and inertia computed from density
- Static, kinematic and dynamic bodies
- Per-axis translation and rotation locks in any frame, including planar bodies
- Substepped soft contact solver with friction, rolling resistance, restitution and warm
  starting
- Speculative contacts for fast bodies
- Continuous collision for bodies knocked or pushed through thin shapes within a single
  step, against static, kinematic and dynamic bodies
- Sleeping islands, and bodies created asleep until touched
- Point, fixed, hinge, slider, distance and swing-twist joints, with limits, motors,
  springs and friction, a force and torque readout, and breaking
- Collision filtering by category and mask
- Ray casts against every shape type, for the closest hit or every hit
- Shape casts and overlap queries for points, boxes and convex geometry, against every
  shape type
- Begin, end and hit events for touching shapes, and sensor shapes that report what
  overlaps them without colliding
- Deterministic single-threaded stepping, and world copies that step identically
- Per-step timings and counts, contact details and sleep islands for debugging tools

Not yet implemented: dynamic meshes.

## Building

Requires CMake 3.21 and a C++17 compiler. Builds with MSVC on Windows, GCC and Clang on
Linux, and Apple Clang on macOS.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release
```

Tests and the viewer build only when Twofold is the top-level project.

On Linux the viewer needs the X11 and OpenGL development packages. On Debian or Ubuntu:

```
sudo apt install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl1-mesa-dev
```

To use Twofold from another CMake project:

```cmake
add_subdirectory(twofold)
target_link_libraries(your_target PRIVATE twofold::twofold)
```

## Viewer

`twofold_viewer` is an interactive debugger for the demo scenes. It uses
[raylib](https://www.raylib.com), [Dear ImGui](https://github.com/ocornut/imgui) and
[rlImGui](https://github.com/raylib-extras/rlImGui), which CMake downloads the first time
you configure. Pass `-DTWOFOLD_BUILD_VIEWER=OFF` to skip it.

- **Timeline.** Every frame of the last minute is recorded. Step forward and back, or drag
  the slider to any recorded frame. Edits, shots and drags are logged, so rewinding and
  playing again reproduces the run exactly. If stepping through recorded frames ever gives
  a different state, the viewer reports the first frame that differs. Verify determinism
  replays the whole run from a fresh scene and compares every frame.
- **Inspector.** Select a body to see and edit its type, pose, velocities, mass, locks,
  damping and materials, its contacts with separations and forces, and plots of its
  height and speed.
- **Bodies.** A sortable table of every body, for finding the fastest, the highest or the
  ones that never fall asleep.
- **Step.** Time spent in each phase of the step, pair and contact counts, and plots of step
  time, kinetic energy, deepest penetration and awake bodies.
- **Draw.** Color by material, motion state, sleep island or speed. Overlays show contact
  points, speculative contacts, normal and friction forces, bounds, centers of mass, body
  axes, velocities, lock planes, joints and body numbers, for all bodies or only the
  selected one.
- **World.** Every world setting can be changed while the scene runs.

| Input | Action |
| --- | --- |
| Left click | Select a body. Drag to pull it |
| B | Shoot a ball at the cursor |
| Space | Play or pause |
| Right, Left | Step forward or back |
| Home, End | First or last recorded frame |
| R | Reset the scene |
| 1 to 9, 0 | Change scene |
| F | Focus the selected body |
| C | Show contacts |
| Delete | Delete the selected body |
| Escape | Clear the selection |
| Right drag | Orbit |
| Middle drag | Pan |
| Wheel | Zoom |

## Example

```cpp
#include <twofold/twofold.h>

using namespace twofold;

World world;

BodyDesc groundDesc;
groundDesc.type = MotionType::Static;
groundDesc.rotation = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });
BodyId ground = world.createBody(groundDesc);
ShapeDesc groundShape;
groundShape.geometry = HalfSpace{};
world.addShape(ground, groundShape);

BodyDesc boxDesc;
boxDesc.position = { 0, 4, 0 };
BodyId box = world.createBody(boxDesc);
ShapeDesc boxShape;
boxShape.geometry = Box{ { 0.5f, 0.5f, 0.5f } };
world.addShape(box, boxShape);

for (int i = 0; i < 120; ++i) {
    world.step(1.0f / 60.0f);
}
Vec3 p = world.position(box);
```

## Planar bodies

Set `BodyDesc::locks` to `LockPlanar` to keep a body in the xy plane of its
`lockFrame`, turning only about that frame's z axis. The default frame is the world
frame. Boxes, prisms and cylinders that extrude along the plane normal use the 2D
path against each other. Any other pair, such as a planar box against a free sphere,
uses the 3D path.

## Meshes

Build a `TriangleMesh` once from vertices and three indices per triangle, then share it
between as many shapes as you like. Winding doesn't matter, since both sides collide.

```cpp
auto terrain = std::make_shared<TriangleMesh>(vertices, indices);

BodyDesc floorDesc;
floorDesc.type = MotionType::Static;
BodyId floor = world.createBody(floorDesc);
ShapeDesc floorShape;
floorShape.geometry = Mesh{ terrain };
world.addShape(floor, floorShape);
```

Meshes go on static and kinematic bodies only. A kinematic mesh moved with
`setKinematicTarget` carries what rests on it.

## Joints

Anchors and axes are given in world space, at the poses the bodies have when the joint is
made. Leave `bodyA` unset to join body B to the world.

```cpp
JointDesc door;
door.bodyB = panel;
HingeJoint hinge{ { 0, 1, 0 }, { 0, 1, 0 } };
hinge.limited = true;
hinge.lowerAngle = -1.6f;
hinge.upperAngle = 1.6f;
door.kind = hinge;
JointId id = world.createJoint(door);
```

A swing-twist joint suits a shoulder or hip. Its swing stays inside an elliptical cone and
its twist about the axis inside a range.

```cpp
JointDesc shoulder;
shoulder.bodyA = torso;
shoulder.bodyB = upperArm;
SwingTwistJoint st;
st.anchor = { 0.2f, 1.5f, 0 };
st.axis = { 1, 0, 0 };
st.normal = { 0, 0, 1 };
st.swingY = 0.8f;
st.swingZ = 0.5f;
st.lowerTwist = -0.5f;
st.upperTwist = 0.5f;
st.frictionTorque = 0.05f;
shoulder.kind = st;
shoulder.breakForce = 500.0f;
JointId arm = world.createJoint(shoulder);
```

Hinges and sliders take `frictionTorque` and `frictionForce`, alongside their motors and
limits. `jointInfo` gives the force and torque a joint applied in the last step. A joint
breaks when either goes over `breakForce` or `breakTorque`. A broken joint stops acting,
shows up in `brokenJoints` for that step, and keeps its handle until destroyed.

Joined bodies don't collide with each other unless `collideConnected` is set, even after
the joint breaks. Destroying a body destroys its joints.

## Queries

Queries never change the world. Each takes a category mask, skips sensors, and results come
in a fixed order, so the same world gives the same answer.

```cpp
std::vector<ShapeId> found;
world.overlapGeometry(Sphere{ 0.5f }, { { 0, 1, 0 } }, found);

CastHit hit;
if (world.castGeometry(Box{ { 0.5f, 0.5f, 0.5f } }, { { 0, 5, 0 } }, { 0, -10, 0 }, hit)) {
    Vec3 rest = Vec3{ 0, 5, 0 } + Vec3{ 0, -10, 0 } * hit.fraction;
}
```

A cast that starts inside a shape reports it at fraction 0 with `initialOverlap` set.
Queries take spheres, boxes, cylinders, prisms and hulls. A mesh has no inside, so only
its triangles count.

## Events

After each step, `contactEvents` lists the shape pairs that began or stopped touching,
`hitEvents` the pairs that struck faster than `WorldSettings::hitThreshold`, and
`sensorEvents` the shapes that entered or left a sensor.

```cpp
ShapeDesc zone;
zone.geometry = Box{ { 2, 1, 2 } };
zone.sensor = true;
ShapeId trigger = world.addShape(ground, zone);

world.step(1.0f / 60.0f);
for (const SensorEvent& e : world.sensorEvents()) {
    if (e.sensor == trigger && e.begin) {
        // e.visitor entered the zone
    }
}
```

A sensor has no mass and never collides or wakes anything. It sees static, kinematic and
sleeping bodies, but not other sensors or the shapes on its own body.

## Units

Twofold works in any length unit. Set `WorldSettings::lengthUnit` to the number of
world units per meter so that tolerances and speed limits scale with the scene.

## License

MIT. See [LICENSE](LICENSE).
