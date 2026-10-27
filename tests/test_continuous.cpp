#include <cmath>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;

BodyId addBody(World& w, MotionType type, const Vec3& pos, const Geometry& g, const Quat& q = {}, uint8_t locks = 0) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    bd.rotation = q;
    bd.locks = locks;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    w.addShape(b, sd);
    return b;
}

World weightless() {
    World w;
    w.setGravity({});
    return w;
}

void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) { w.step(kDt); }
}

}

TEST(knocked_box_stays_behind_thin_wall) {
    for (int trial = 0; trial < 20; ++trial) {
        World w = weightless();
        addBody(w, MotionType::Static, {}, Box{ { 0.005f, 2, 2 } });
        const BodyId target = addBody(w, MotionType::Dynamic, { -0.45f - 0.01f * trial, 0, 0 }, Box{ { 0.1f, 0.1f, 0.1f } },
                                      Quat::axisAngle({ 0.3f, 1, 0.2f }, 0.3f * (float)trial));
        const BodyId hammer = addBody(w, MotionType::Dynamic, { -2.0f, 0.03f * trial - 0.3f, 0.01f * trial },
                                      Box{ { 0.2f, 0.2f, 0.2f } });
        w.setLinearVelocity(hammer, { 60, 0, 0 });
        run(w, 60);
        CHECK(w.position(target).x < 0.0f);
    }
}

TEST(knocked_box_stays_in_front_of_thin_dynamic_plate) {
    for (int trial = 0; trial < 20; ++trial) {
        World w = weightless();
        const BodyId plate = addBody(w, MotionType::Dynamic, {}, Box{ { 0.005f, 1, 1 } });
        const BodyId target = addBody(w, MotionType::Dynamic, { -0.45f - 0.01f * trial, 0, 0 }, Box{ { 0.1f, 0.1f, 0.1f } });
        const BodyId hammer = addBody(w, MotionType::Dynamic, { -2.0f, 0.02f * trial, 0 }, Box{ { 0.2f, 0.2f, 0.2f } });
        w.setLinearVelocity(hammer, { 60, 0, 0 });
        run(w, 2);
        CHECK(w.transform(plate).applyInverse(w.position(target)).x < 0.0f);
    }
}

TEST(knocked_planar_layer_stays_behind_thin_band) {
    const auto square = [](float h) { return Prism{ { { -h, -h }, { h, -h }, { h, h }, { -h, h } }, -0.05f, 0.05f }; };
    for (int trial = 0; trial < 20; ++trial) {
        World w = weightless();
        const Prism band{ { { -0.005f, -2 }, { 0.005f, -2 }, { 0.005f, 2 }, { -0.005f, 2 } }, -0.05f, 0.05f };
        addBody(w, MotionType::Static, {}, band, {}, LockPlanar);
        const BodyId target = addBody(w, MotionType::Dynamic, { -0.45f - 0.01f * trial, 0, 0 }, square(0.1f), {}, LockPlanar);
        const BodyId hammer = addBody(w, MotionType::Dynamic, { -2.0f, 0.02f * trial, 0 }, square(0.2f), {}, LockPlanar);
        w.setLinearVelocity(hammer, { 60, 0, 0 });
        run(w, 60);
        CHECK(w.position(target).x < 0.0f);
    }
}

TEST(strong_force_cannot_push_box_through_thin_wall) {
    for (int trial = 0; trial < 20; ++trial) {
        World w = weightless();
        addBody(w, MotionType::Static, {}, Box{ { 0.005f, 2, 2 } });
        const BodyId b = addBody(w, MotionType::Dynamic, { -0.5f - 0.03f * trial, 0, 0 }, Box{ { 0.05f, 0.05f, 0.05f } });
        const float m = w.massProperties(b).mass;
        for (int i = 0; i < 60; ++i) {
            w.addForce(b, { m * 6000.0f, 0, 0 });
            w.step(kDt);
        }
        CHECK(w.position(b).x < 0.0f);
    }
}

TEST(heavy_gravity_cannot_pull_sphere_through_thin_floor) {
    for (int trial = 0; trial < 20; ++trial) {
        World w;
        w.setGravity({ 0, -3000, 0 });
        addBody(w, MotionType::Static, {}, Box{ { 3, 0.005f, 3 } });
        const BodyId b = addBody(w, MotionType::Dynamic, { 0, 0.3f + 0.037f * trial, 0 }, Sphere{ 0.05f });
        run(w, 30);
        CHECK(w.position(b).y > 0.0f);
    }
}

TEST(slow_bodies_are_not_swept) {
    World w;
    addBody(w, MotionType::Static, { 0, -0.5f, 0 }, Box{ { 5, 0.5f, 5 } });
    for (int i = 0; i < 5; ++i) { addBody(w, MotionType::Dynamic, { 0, 0.6f + 1.1f * (float)i, 0 }, Box{}); }
    for (int i = 0; i < 120; ++i) {
        w.step(kDt);
        CHECK(w.profile().sweptPairs == 0);
    }
}
