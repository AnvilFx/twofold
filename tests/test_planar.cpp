#include <cmath>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

// A floor through the origin, leaning toward +z by tilt radians.
void addFloor(World& w, bool halfSpace, float tilt) {
    const Quat lean = Quat::axisAngle({ 1, 0, 0 }, tilt);
    BodyDesc bd;
    bd.type = MotionType::Static;
    ShapeDesc sd;
    sd.friction = 0.6f;
    if (halfSpace) {
        bd.rotation = lean * kUp;
        sd.geometry = HalfSpace{};
    } else {
        bd.rotation = lean;
        bd.position = lean.rotate({ 0, -0.5f, 0 });
        sd.geometry = Box{ { 50, 0.5f, 2 } };
    }
    w.addShape(w.createBody(bd), sd);
}

BodyId addLayer(World& w, MotionType type, const Vec3& pos, const Geometry& g) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    w.addShape(b, sd);
    return b;
}

struct Landing {
    Vec3 position;
    Vec3 velocity;
    float spin = 0.0f;
    // Largest component of a contact normal along the plane normal.
    float offPlane = 0.0f;
};

Landing land(bool halfSpace, float tilt, const Geometry& g, const BodyDesc& start, float rollingResistance = 0.0f) {
    World w;
    addFloor(w, halfSpace, tilt);
    BodyDesc bd = start;
    bd.locks = LockPlanar;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    sd.friction = 0.6f;
    sd.rollingResistance = rollingResistance;
    w.addShape(b, sd);
    Landing out;
    for (int i = 0; i < 600; ++i) {
        w.step(kDt);
        for (const ContactInfo& c : w.contacts()) { out.offPlane = std::fmax(out.offPlane, std::fabs(c.normal.z)); }
    }
    out.position = w.position(b);
    out.velocity = w.linearVelocity(b);
    out.spin = w.angularVelocity(b).z;
    return out;
}

}

TEST(planar_box_lands_on_tilted_half_space_as_on_a_box_floor) {
    const Box box{ { 0.5f, 0.3f, 0.1f } };
    BodyDesc tumble;
    tumble.position = { 0, 1, 0 };
    tumble.linearVelocity = { 4, 0, 0 };
    tumble.angularVelocity = { 0, 0, 5 };
    BodyDesc slide;
    slide.position = { 0, 0.3f, 0 };
    slide.linearVelocity = { 3, 0, 0 };
    for (float tilt : { 0.0f, 0.005f }) {
        for (const BodyDesc& start : { tumble, slide }) {
            const Landing floor = land(false, tilt, box, start);
            const Landing half = land(true, tilt, box, start);
            CHECK_NEAR(half.position.x, floor.position.x, 1e-3f);
            CHECK_NEAR(half.position.y, floor.position.y, 1e-4f);
            CHECK(half.offPlane == 0.0f);
            CHECK(length(half.velocity) < 1e-3f);
        }
    }
}

TEST(planar_disc_rolls_on_tilted_half_space_as_on_a_box_floor) {
    const Cylinder disc{ 0.5f, 0.1f };
    BodyDesc roll;
    roll.position = { 0, 0.5f, 0 };
    roll.linearVelocity = { 2, 0, 0 };
    roll.angularVelocity = { 0, 0, -4 };
    BodyDesc spin;
    spin.position = { 0, 0.5f, 0 };
    spin.angularVelocity = { 0, 0, -6 };
    for (float tilt : { 0.0f, 0.005f }) {
        const Landing floor = land(false, tilt, disc, roll, 0.1f);
        const Landing half = land(true, tilt, disc, roll, 0.1f);
        CHECK_NEAR(half.position.x, floor.position.x, 1e-3f);
        CHECK_NEAR(half.position.y, 0.5f, 1e-4f);
        CHECK(half.offPlane == 0.0f);

        const Landing turning = land(true, tilt, disc, spin);
        CHECK_NEAR(turning.spin, -2.0f, 0.02f);
        CHECK_NEAR(turning.velocity.x, 1.0f, 0.01f);
        CHECK_NEAR(turning.position.y, 0.5f, 1e-4f);
    }
}

TEST(ray_hits_the_first_of_coincident_layers) {
    World w;
    const Prism layer{ { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } }, -0.05f, 0.05f };
    BodyId first;
    for (int i = 0; i < 3; ++i) {
        const BodyId b = addLayer(w, MotionType::Static, { 0, 1, 0 }, layer);
        if (i == 0) { first = b; }
    }
    RayHit hit;
    CHECK(w.castRay({ 0.3f, 1.2f, 5 }, { 0, 0, -10 }, hit));
    CHECK(hit.body == first);
    CHECK(w.castRay({ 5, 1.3f, 0.01f }, { -10, 0, 0 }, hit));
    CHECK(hit.body == first);

    const BodyId front = addLayer(w, MotionType::Static, { 0, 1, 0.2f }, layer);
    CHECK(w.castRay({ 0.3f, 1.2f, 5 }, { 0, 0, -10 }, hit));
    CHECK(hit.body == front);
}
