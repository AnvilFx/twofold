# Twofold

Rigid body physics for scenes that mix 2D and 3D bodies in one world.

A body locked to a plane is still a 3D body, so it collides with everything else in the
world. Two planar bodies that share a plane collide through a dedicated 2D path: their
profiles are tested as polygons and circles, and contacts never push them apart along
the plane normal.

## Who does what

The code, tests and documentation are written by Claude, Anthropic's AI model. AnvilFx
does not write code. AnvilFx directs the approaches and techniques the engine uses and
makes sure they are modern and fast, rather than leaving those choices to Claude.

Every claim Claude makes about a change, such as that it works or that it is faster,
is tested by AnvilFx before it ships.

## Features

- Sphere, box, cylinder, extruded convex polygon and half-space shapes
- Bodies made of several shapes, with mass and inertia computed from density
- Static, kinematic and dynamic bodies
- Per-axis translation and rotation locks in any frame, including planar bodies
- Substepped soft contact solver with friction, restitution and warm starting
- Speculative contacts for fast bodies
- Sleeping islands, and bodies created asleep until touched
- Collision filtering by category and mask
- Ray casts against every shape type
- Deterministic single-threaded stepping, and world copies that step identically
- Per-step timings and counts, contact details and sleep islands for debugging tools

Not yet implemented: joints, convex hulls, triangle meshes, shape casts and overlap
queries, and contact callbacks. Cylinders collide as 24-sided prisms, except against a
planar body in the same plane, where they are exact circles.

## Building

Requires CMake 3.21 and a C++17 compiler.

```
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Tests build only when Twofold is the top-level project.

To use Twofold from another CMake project:

```cmake
add_subdirectory(twofold)
target_link_libraries(your_target PRIVATE twofold::twofold)
```

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

## Units

Twofold works in any length unit. Set `WorldSettings::lengthUnit` to the number of
world units per meter so that tolerances and speed limits scale with the scene.

## License

MIT. See [LICENSE](LICENSE).
