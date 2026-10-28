#include <chrono>
#include <cmath>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kPi = 3.14159265358979f;

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

BodyId addGround(World& w, float friction = 0.6f) {
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    sd.friction = friction;
    w.addShape(g, sd);
    return g;
}

BodyId addBox(World& w, const Vec3& pos, const Vec3& half, MotionType type = MotionType::Dynamic) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ half };
    sd.friction = 0.6f;
    w.addShape(b, sd);
    return b;
}

BodyId addSphere(World& w, const Vec3& pos, float r, float restitution = 0.0f) {
    BodyDesc bd;
    bd.position = pos;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Sphere{ r };
    sd.restitution = restitution;
    w.addShape(b, sd);
    return b;
}

void run(World& w, float seconds) {
    const int steps = (int)std::lround(seconds / kDt);
    for (int i = 0; i < steps; ++i) { w.step(kDt); }
}

float upAxisY(const Quat& q) { return std::fabs(q.rotate({ 0, 1, 0 }).y); }

float deepestContact(const World& w) {
    float d = 0.0f;
    for (const ContactInfo& c : w.contacts()) {
        for (int i = 0; i < c.pointCount; ++i) { d = std::fmax(d, -c.separations[i]); }
    }
    return d;
}

}

TEST(mass_of_box_prism_and_cylinder) {
    World w;
    BodyDesc bd;
    const BodyId box = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 1.0f, 1.5f } };
    sd.density = 2.0f;
    w.addShape(box, sd);
    const MassProperties mb = w.massProperties(box);
    CHECK_NEAR(mb.mass, 12.0f, 1e-4f);
    CHECK_NEAR(mb.inertia.r[0].x, 12.0f * (4.0f + 9.0f) / 12.0f, 1e-3f);
    CHECK_NEAR(mb.inertia.r[2].z, 12.0f * (1.0f + 4.0f) / 12.0f, 1e-3f);
    CHECK_NEAR(mb.inertia.r[0].y, 0.0f, 1e-4f);

    const BodyId prism = w.createBody(bd);
    sd.geometry = Prism{ { { -0.5f, -1.0f }, { -0.5f, 1.0f }, { 0.5f, 1.0f }, { 0.5f, -1.0f } }, -1.5f, 1.5f };
    w.addShape(prism, sd);
    const MassProperties mp = w.massProperties(prism);
    CHECK_NEAR(mp.mass, mb.mass, 1e-4f);
    CHECK_NEAR(mp.inertia.r[1].y, mb.inertia.r[1].y, 1e-3f);

    const BodyId cyl = w.createBody(bd);
    sd.geometry = Cylinder{ 1.0f, 0.5f };
    sd.density = 1.0f;
    w.addShape(cyl, sd);
    CHECK_NEAR(w.massProperties(cyl).mass, kPi, 1e-4f);
}

TEST(offset_shape_moves_center_of_mass) {
    World w;
    const BodyId b = w.createBody({});
    ShapeDesc sd;
    sd.geometry = Sphere{ 0.5f };
    sd.position = { 1, 0, 0 };
    w.addShape(b, sd);
    sd.position = { -1, 0, 0 };
    sd.density = 3.0f;
    w.addShape(b, sd);
    CHECK_NEAR(w.massProperties(b).center.x, -0.5f, 1e-4f);
}

TEST(box_rests_on_half_space) {
    World w;
    addGround(w);
    const BodyId b = addBox(w, { 0, 2, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 3.0f);
    CHECK_NEAR(w.position(b).y, 0.5f, 0.02f);
    CHECK(upAxisY(w.rotation(b)) > 0.999f);
}

TEST(box_rests_on_static_box) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 10, 0.5f, 10 }, MotionType::Static);
    const BodyId b = addBox(w, { 0, 2, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 3.0f);
    CHECK_NEAR(w.position(b).y, 0.5f, 0.02f);
    CHECK_NEAR(w.position(b).x, 0.0f, 0.01f);
}

TEST(sphere_rests_on_static_box) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 10, 0.5f, 10 }, MotionType::Static);
    const BodyId s = addSphere(w, { 0, 2, 0 }, 0.5f);
    run(w, 3.0f);
    CHECK_NEAR(w.position(s).y, 0.5f, 0.02f);
}

TEST(tilted_box_settles_flat) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 3, 0 };
    bd.rotation = Quat::axisAngle({ 1, 0.3f, 0.7f }, 0.9f);
    bd.angularVelocity = { 1, 2, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    run(w, 5.0f);
    CHECK_NEAR(w.position(b).y, 0.5f, 0.03f);
    const Quat q = w.rotation(b);
    const float best = std::fmax(std::fabs(q.rotate({ 1, 0, 0 }).y),
                                 std::fmax(std::fabs(q.rotate({ 0, 1, 0 }).y), std::fabs(q.rotate({ 0, 0, 1 }).y)));
    CHECK(best > 0.995f);
}

TEST(crossed_edges_do_not_pass_through) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.25f * kPi);
    const BodyId base = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(base, sd);

    BodyDesc fall;
    fall.position = { 0, 2, 0 };
    fall.rotation = Quat::axisAngle({ 1, 0, 0 }, 0.25f * kPi);
    fall.linearVelocity = { 0, -4, 0 };
    fall.locks = LockRotationX | LockRotationY | LockRotationZ;
    const BodyId b = w.createBody(fall);
    w.addShape(b, sd);
    run(w, 1.0f);
    const std::vector<ContactInfo> cs = w.contacts();
    CHECK(cs.size() == 1);
    CHECK(cs[0].pointCount == 1);
    CHECK(cs[0].separations[0] > -0.02f);
    CHECK_NEAR(w.position(b).y, std::sqrt(2.0f), 0.02f);
}

TEST(box_across_slab_edge_does_not_sink) {
    for (int trial = 0; trial < 40; ++trial) {
        World w;
        BodyDesc sd;
        sd.type = MotionType::Static;
        sd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.4f);
        const BodyId slab = w.createBody(sd);
        ShapeDesc ss;
        ss.geometry = Box{ { 3, 0.3f, 2 } };
        w.addShape(slab, ss);
        const float t = (float)trial;
        BodyDesc bd;
        bd.position = { 0.07f * t - 1.4f, 2.0f, 1.9f + 0.009f * t };
        bd.rotation = Quat::axisAngle(normalize(Vec3{ std::sin(t), std::cos(1.7f * t), std::sin(2.3f * t) }), 0.37f * t);
        const BodyId b = w.createBody(bd);
        ShapeDesc bs;
        const Vec3 half{ 0.2f + 0.3f * std::fabs(std::sin(0.9f * t)), 0.2f + 0.3f * std::fabs(std::cos(1.3f * t)), 0.35f };
        bs.geometry = Box{ half };
        w.addShape(b, bs);
        float worst = 0.0f;
        for (int i = 0; i < 240; ++i) {
            w.step(kDt);
            if (i >= 30) { worst = std::fmax(worst, deepestContact(w)); }
        }
        CHECK(worst < 0.02f);
    }
}

TEST(spinning_cylinders_meet_without_sinking) {
    float worst = 0.0f;
    for (int trial = 0; trial < 60; ++trial) {
        World w;
        addBox(w, { 0, -0.5f, 0 }, { 4, 0.5f, 4 }, MotionType::Static);
        const float t = (float)trial;
        const Quat q = Quat::axisAngle(normalize(Vec3{ std::sin(t), std::cos(1.7f * t), std::sin(2.3f * t) }), 0.7f * t);
        const Vec3 spin = normalize(Vec3{ std::cos(1.3f * t), 0.3f * std::sin(t), std::sin(1.3f * t) }) * 25.0f;
        for (int k = 0; k < 2; ++k) {
            const float side = k ? 1.0f : -1.0f;
            BodyDesc bd;
            bd.position = { 0.7f * side, 0.6f + (float)k * 0.1f * std::sin(t), 0 };
            bd.rotation = k ? q.conjugate() : q;
            const BodyId b = w.createBody(bd);
            ShapeDesc sd;
            sd.geometry = Cylinder{ 0.3f, 0.2f };
            w.addShape(b, sd);
            w.setAngularVelocity(b, spin * -side);
            w.setLinearVelocity(b, { -3.0f * side, 0, 0 });
        }
        for (int i = 0; i < 90; ++i) {
            w.step(kDt);
            worst = std::fmax(worst, deepestContact(w));
        }
    }
    CHECK(worst < 0.02f);
}

TEST(box_stack_stays_upright) {
    World w;
    addGround(w);
    BodyId top;
    for (int i = 0; i < 10; ++i) { top = addBox(w, { 0, 0.5f + (float)i, 0 }, { 0.5f, 0.5f, 0.5f }); }
    run(w, 6.0f);
    CHECK_NEAR(w.position(top).y, 9.5f, 0.1f);
    CHECK_NEAR(w.position(top).x, 0.0f, 0.05f);
    CHECK_NEAR(w.position(top).z, 0.0f, 0.05f);
}

TEST(bounce_with_full_restitution) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    sd.restitution = 1.0f;
    w.addShape(g, sd);
    const BodyId s = addSphere(w, { 0, 3, 0 }, 0.5f, 1.0f);
    float peak = 0.0f;
    bool landed = false;
    for (int i = 0; i < 240; ++i) {
        w.step(kDt);
        const float y = w.position(s).y;
        if (y < 0.6f) { landed = true; }
        if (landed) { peak = std::fmax(peak, y); }
    }
    CHECK(landed);
    CHECK(peak > 2.6f);
}

TEST(friction_holds_on_slope) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.35f) * kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    sd.friction = 0.8f;
    w.addShape(g, sd);

    BodyDesc boxDesc;
    boxDesc.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.35f);
    boxDesc.position = boxDesc.rotation.rotate({ 0, 0.5f, 0 });
    const BodyId b = w.createBody(boxDesc);
    ShapeDesc bs;
    bs.geometry = Box{};
    bs.friction = 0.8f;
    w.addShape(b, bs);
    const Vec3 start = w.position(b);
    run(w, 2.0f);
    CHECK(length(w.position(b) - start) < 0.05f);
}

TEST(low_friction_slides_down_slope) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.35f) * kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    sd.friction = 0.05f;
    w.addShape(g, sd);

    BodyDesc boxDesc;
    boxDesc.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.35f);
    boxDesc.position = boxDesc.rotation.rotate({ 0, 0.5f, 0 });
    const BodyId b = w.createBody(boxDesc);
    ShapeDesc bs;
    bs.geometry = Box{};
    bs.friction = 0.05f;
    w.addShape(b, bs);
    const Vec3 start = w.position(b);
    run(w, 1.0f);
    CHECK(w.position(b).x < start.x - 0.5f);
}

TEST(body_falls_asleep) {
    World w;
    addGround(w);
    const BodyId b = addBox(w, { 0, 1, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 4.0f);
    CHECK(!w.isAwake(b));
}

TEST(near_neighbors_do_not_toggle_sleep) {
    World w;
    addGround(w);
    const BodyId a = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    const BodyId b = addBox(w, { 1.01f, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    w.setCanSleep(a, false);
    run(w, 3.0f);
    int toggles = 0;
    bool last = w.isAwake(b);
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        if (w.isAwake(b) != last) { ++toggles; }
        last = w.isAwake(b);
    }
    CHECK(toggles == 0);
}

TEST(sleeping_body_wakes_when_hit) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.5f, 0 };
    bd.awake = false;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    CHECK(!w.isAwake(b));
    run(w, 0.5f);
    CHECK(!w.isAwake(b));

    BodyDesc ball;
    ball.position = { -3, 0.5f, 0 };
    ball.linearVelocity = { 6, 0, 0 };
    const BodyId s = w.createBody(ball);
    sd.geometry = Sphere{ 0.5f };
    w.addShape(s, sd);
    run(w, 1.0f);
    CHECK(w.position(b).x > 0.1f);
}

TEST(dormant_body_floats_until_touched) {
    World w;
    BodyDesc bd;
    bd.position = { 0, 5, 0 };
    bd.awake = false;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    w.addForce(b, { 0, 100, 0 }, false);
    run(w, 1.0f);
    CHECK_NEAR(w.position(b).y, 5.0f, 1e-6f);
}

TEST(kinematic_body_pushes_dynamic) {
    World w;
    addGround(w, 0.3f);
    const BodyId k = addBox(w, { -2, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f }, MotionType::Kinematic);
    const BodyId d = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    for (int i = 0; i < 120; ++i) {
        const float x = -2.0f + 3.0f * (float)(i + 1) / 120.0f;
        w.setKinematicTarget(k, { x, 0.5f, 0 }, {}, kDt);
        w.step(kDt);
    }
    CHECK_NEAR(w.position(k).x, 1.0f, 1e-4f);
    CHECK(w.position(d).x > 1.9f);
}

TEST(rotation_locks_hold) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.5f, 0 };
    bd.locks = LockRotationX | LockRotationY | LockRotationZ;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    w.addImpulseAtPoint(b, { 5, 0, 0 }, { 0, 1, 0.5f });
    run(w, 1.0f);
    CHECK(upAxisY(w.rotation(b)) > 0.99999f);
    CHECK(w.position(b).x > 0.1f);
}

TEST(planar_bodies_stay_in_plane) {
    World w;
    addGround(w);
    BodyId top;
    for (int i = 0; i < 6; ++i) {
        BodyDesc bd;
        bd.position = { 0.05f * (float)i, 0.5f + (float)i, 0.004f * (float)(i % 2) };
        bd.locks = LockPlanar;
        top = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Box{ { 0.5f, 0.5f, 0.01f } };
        w.addShape(top, sd);
    }
    run(w, 4.0f);
    CHECK_NEAR(w.position(top).z, 0.004f, 1e-6f);
    CHECK_NEAR(w.position(top).y, 5.5f, 0.1f);
    const Vec3 n = w.rotation(top).rotate({ 0, 0, 1 });
    CHECK(n.z > 0.99999f);
}

TEST(planar_thin_layers_do_not_separate_in_depth) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 0.5f, 0 };
    const BodyId a = w.createBody(bd);
    bd.position = { 0.6f, 0.5f, 0.001f };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 0.5f, 0.002f } };
    w.addShape(a, sd);
    w.addShape(b, sd);
    run(w, 1.0f);
    CHECK(w.position(b).x - w.position(a).x > 0.99f);
    CHECK_NEAR(w.position(b).z, 0.001f, 1e-6f);
}

TEST(planar_box_slides_to_rest_on_half_space) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 0.5f, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 0.5f, 0.1f } };
    w.addShape(b, sd);
    w.setLinearVelocity(b, { 3, 0, 0 });
    run(w, 2.0f);
    CHECK(w.position(b).x < 1.0f);
    CHECK(std::fabs(w.linearVelocity(b).x) < 0.01f);
}

// A spinning disc set down settles into rolling at a third of its spin.
TEST(planar_disc_rolls_on_box) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 10, 0.5f, 1 }, MotionType::Static);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 0.5f, 0 };
    const BodyId wheel = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.5f, 0.1f };
    w.addShape(wheel, sd);
    w.setAngularVelocity(wheel, { 0, 0, -6 });
    run(w, 1.0f);
    CHECK_NEAR(w.angularVelocity(wheel).z, -2.0f, 0.02f);
    CHECK_NEAR(w.linearVelocity(wheel).x, 1.0f, 0.01f);
    CHECK_NEAR(w.position(wheel).y, 0.5f, 0.01f);
}

TEST(planar_body_meets_3d_body) {
    World w;
    addGround(w, 0.2f);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 0.5f, 0 };
    const BodyId flat = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 0.5f, 0.05f } };
    w.addShape(flat, sd);

    BodyDesc ball;
    ball.position = { -3, 0.4f, 0 };
    ball.linearVelocity = { 8, 0, 0 };
    const BodyId s = w.createBody(ball);
    sd.geometry = Sphere{ 0.4f };
    sd.density = 4.0f;
    w.addShape(s, sd);

    run(w, 1.0f);
    CHECK(w.position(flat).x > 0.3f);
    CHECK_NEAR(w.position(flat).z, 0.0f, 1e-6f);
}

TEST(rotated_lock_frame_defines_the_plane) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.lockFrame = Quat::axisAngle({ 0, 1, 0 }, 0.5f * kPi);
    bd.position = { 0, 2, 0 };
    bd.linearVelocity = { 3, 0, 2 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    run(w, 1.0f);
    CHECK_NEAR(w.position(b).x, 0.0f, 1e-5f);
    CHECK(w.position(b).z > 1.0f);
}

TEST(identical_runs_match_exactly) {
    auto simulate = []() {
        World w;
        addGround(w);
        std::vector<BodyId> ids;
        for (int i = 0; i < 20; ++i) {
            BodyDesc bd;
            bd.position = { 0.3f * (float)(i % 4), 1.0f + 1.1f * (float)i, 0.2f * (float)(i % 3) };
            bd.rotation = Quat::axisAngle({ 1, 1, 0 }, 0.1f * (float)i);
            const BodyId b = w.createBody(bd);
            ShapeDesc sd;
            if (i % 2) {
                sd.geometry = Sphere{ 0.4f };
            } else {
                sd.geometry = Box{ { 0.4f, 0.3f, 0.5f } };
            }
            w.addShape(b, sd);
            ids.push_back(b);
        }
        run(w, 3.0f);
        std::vector<float> out;
        for (BodyId b : ids) {
            const Transform t = w.transform(b);
            out.insert(out.end(), { t.p.x, t.p.y, t.p.z, t.q.x, t.q.y, t.q.z, t.q.w });
        }
        return out;
    };
    CHECK(simulate() == simulate());
}

TEST(destroyed_ids_become_invalid) {
    World w;
    const BodyId a = addBox(w, { 0, 0, 0 }, { 0.5f, 0.5f, 0.5f });
    w.destroyBody(a);
    CHECK(!w.isValid(a));
    const BodyId b = addBox(w, { 0, 0, 0 }, { 0.5f, 0.5f, 0.5f });
    CHECK(b.index == a.index);
    CHECK(w.isValid(b));
    CHECK(!w.isValid(a));
}

TEST(filters_skip_collision) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    w.addShape(g, sd);

    const BodyId b = w.createBody({});
    ShapeDesc bs;
    bs.geometry = Box{};
    bs.mask = 0;
    w.addShape(b, bs);
    w.setTransform(b, { 0, 1, 0 }, {});
    run(w, 1.0f);
    CHECK(w.position(b).y < -1.0f);
}

TEST(contacts_are_reported) {
    World w;
    const BodyId g = addGround(w);
    const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 0.5f);
    const std::vector<ContactInfo> cs = w.contacts();
    CHECK(cs.size() == 1);
    CHECK((cs[0].bodyA == g && cs[0].bodyB == b) || (cs[0].bodyA == b && cs[0].bodyB == g));
    CHECK(cs[0].pointCount == 4);
}

TEST(copied_world_steps_identically) {
    World w;
    addGround(w);
    std::vector<BodyId> ids;
    for (int i = 0; i < 10; ++i) {
        ids.push_back(addBox(w, { 0.2f * (float)(i % 3), 0.6f + 1.05f * (float)i, 0 }, { 0.5f, 0.5f, 0.5f }));
    }
    run(w, 0.5f);
    World copy = w;
    run(w, 1.0f);
    run(copy, 1.0f);
    for (BodyId b : ids) {
        const Transform a = w.transform(b), c = copy.transform(b);
        CHECK(a.p.x == c.p.x && a.p.y == c.p.y && a.p.z == c.p.z);
        CHECK(a.q.x == c.q.x && a.q.y == c.q.y && a.q.z == c.q.z && a.q.w == c.q.w);
    }
}

TEST(ray_hits_closest_shape) {
    World w;
    const BodyId box = addBox(w, { 5, 0, 0 }, { 0.5f, 0.5f, 0.5f }, MotionType::Static);
    addSphere(w, { 10, 0, 0 }, 1.0f);
    RayHit hit;
    CHECK(w.castRay({ 0, 0, 0 }, { 20, 0, 0 }, hit));
    CHECK(hit.body == box);
    CHECK_NEAR(hit.fraction, 4.5f / 20.0f, 1e-5f);
    CHECK_NEAR(hit.point.x, 4.5f, 1e-4f);
    CHECK_NEAR(hit.normal.x, -1.0f, 1e-5f);

    CHECK(w.castRay({ 10, 5, 0 }, { 0, -10, 0 }, hit));
    CHECK_NEAR(hit.point.y, 1.0f, 1e-4f);
    CHECK_NEAR(hit.normal.y, 1.0f, 1e-4f);
    CHECK(!w.castRay({ 0, 0, 0 }, { 4, 0, 0 }, hit));
    CHECK(!w.castRay({ 0, 0, 0 }, { 20, 0, 0 }, hit, 0));
}

TEST(ray_hits_ground_and_ignores_inner_origin) {
    World w;
    const BodyId g = addGround(w);
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.position = { 0, 2, 0 };
    bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.25f * kPi);
    const BodyId tilted = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(tilted, sd);

    RayHit hit;
    CHECK(w.castRay({ 3, 5, 0 }, { 0, -10, 0 }, hit));
    CHECK(hit.body == g);
    CHECK_NEAR(hit.point.y, 0.0f, 1e-4f);
    CHECK(w.castRay({ 0, 5, 0 }, { 0, -10, 0 }, hit));
    CHECK(hit.body == tilted);
    CHECK_NEAR(hit.point.y, 2.0f + std::sqrt(0.5f), 1e-4f);
    CHECK(w.castRay({ 0, 2, 0 }, { 0, -10, 0 }, hit));
    CHECK(hit.body == g);
}

TEST(islands_follow_contacts) {
    World w;
    const BodyId g = addGround(w);
    const BodyId lower = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    const BodyId upper = addBox(w, { 0, 1.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    const BodyId apart = addBox(w, { 5, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 0.1f);
    CHECK(w.island(g) == UINT32_MAX);
    CHECK(w.island(lower) == w.island(upper));
    CHECK(w.island(lower) != w.island(apart));
    CHECK(w.profile().awakeBodies == 3);
    CHECK(w.profile().awakeIslands == 2);
    CHECK(w.profile().touchingContacts == 3);
    CHECK(w.profile().total >= w.profile().solver);
}

TEST(speculative_contacts_are_reported_on_request) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.51f, 0 };
    bd.gravityScale = 0.0f;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{};
    w.addShape(b, sd);
    w.step(kDt);
    CHECK(w.contacts().empty());
    const std::vector<ContactInfo> cs = w.contacts(true);
    CHECK(cs.size() == 1);
    CHECK(cs[0].separations[0] > 0.0f);
}

TEST(planar_near_miss_reports_the_gap) {
    World w;
    w.setGravity({});
    BodyDesc ad;
    ad.type = MotionType::Static;
    ad.locks = LockPlanar;
    const BodyId a = w.createBody(ad);
    ShapeDesc sa;
    sa.geometry = Prism{ { { 0.296f, 0.338f }, { 0.315f, 0.376f }, { 0.152f, 0.466f } }, -0.05f, 0.05f };
    w.addShape(a, sa);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.linearVelocity = { 3.4f, -2.0f, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sb;
    sb.geometry = Prism{ { { 0.125f, 0.570f }, { 0.122f, 0.587f }, { 0.062f, 0.533f }, { 0.096f, 0.528f } }, -0.05f, 0.05f };
    w.addShape(b, sb);
    w.step(kDt);
    const std::vector<ContactInfo> cs = w.contacts(true);
    CHECK(cs.size() == 1);
    for (const ContactInfo& c : cs) {
        for (int i = 0; i < c.pointCount; ++i) { CHECK(c.separations[i] > 0.08f); }
    }
}

TEST(cylinder_on_its_side_comes_to_rest) {
    World w;
    addGround(w, 0.5f);
    BodyDesc bd;
    bd.position = { 0, 0.6f, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.3f, 0.3f };
    w.addShape(b, sd);
    run(w, 1.0f);
    float deepest = 0.0f;
    for (int f = 0; f < 120; ++f) {
        w.step(kDt);
        deepest = std::fmax(deepest, deepestContact(w));
    }
    CHECK(!w.isAwake(b));
    CHECK(deepest < 0.005f);
}

TEST(tumbling_cylinders_do_not_sink_into_ground) {
    float deepest = 0.0f;
    for (int i = 0; i < 24; ++i) {
        World w;
        addGround(w);
        const float t = (float)i;
        BodyDesc bd;
        bd.position = { 0, 2.0f + 0.25f * t, 0 };
        bd.rotation = Quat::axisAngle({ std::sin(1.3f * t), std::cos(2.1f * t), 0.5f }, 0.4f * t);
        bd.angularVelocity = Vec3{ std::sin(0.7f * t), std::cos(1.9f * t), std::sin(2.9f * t) } * 8.0f;
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Cylinder{ 0.3f, 0.3f };
        w.addShape(b, sd);
        for (int f = 0; f < 240; ++f) {
            w.step(kDt);
            deepest = std::fmax(deepest, deepestContact(w));
        }
    }
    CHECK(deepest < 0.02f);
}

TEST(cylinder_rolls_on_its_side_at_a_constant_height) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.3f, 0 };
    bd.linearVelocity = { 2, 0, 0 };
    bd.angularVelocity = { 0, 0, -2.0f / 0.3f };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.3f, 0.3f };
    w.addShape(b, sd);
    float lowest = 1.0f, highest = 0.0f;
    for (int f = 0; f < 120; ++f) {
        w.step(kDt);
        lowest = std::fmin(lowest, w.position(b).y);
        highest = std::fmax(highest, w.position(b).y);
    }
    CHECK(lowest > 0.298f && highest < 0.302f);
    CHECK_NEAR(w.linearVelocity(b).x, 2.0f, 0.02f);
}

TEST(cylinders_rest_on_a_box_lying_and_standing) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 3, 0.5f, 3 }, MotionType::Static);
    BodyDesc lying;
    lying.position = { -1, 0.35f, 0 };
    const BodyId a = w.createBody(lying);
    BodyDesc standing;
    standing.position = { 1, 0.45f, 0 };
    standing.rotation = kUp;
    const BodyId b = w.createBody(standing);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.3f, 0.4f };
    w.addShape(a, sd);
    w.addShape(b, sd);
    run(w, 2.0f);
    CHECK(!w.isAwake(a) && !w.isAwake(b));
    CHECK(deepestContact(w) < 0.005f);
    CHECK_NEAR(w.position(a).y, 0.3f, 0.005f);
    CHECK_NEAR(w.position(b).y, 0.4f, 0.005f);
}

TEST(standing_cylinders_stack) {
    World w;
    addGround(w);
    std::vector<BodyId> ids;
    for (int i = 0; i < 4; ++i) {
        BodyDesc bd;
        bd.position = { 0.03f * (float)(i % 2), 0.3f + 0.61f * (float)i, 0 };
        bd.rotation = kUp;
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Cylinder{ 0.4f, 0.3f };
        w.addShape(b, sd);
        ids.push_back(b);
    }
    run(w, 3.0f);
    for (BodyId b : ids) { CHECK(!w.isAwake(b)); }
    CHECK(deepestContact(w) < 0.005f);
    CHECK_NEAR(w.position(ids.back()).y, 2.1f, 0.01f);
}

// Tilted cylinders a few millimeters into each other. Contacts must not report an overlap much deeper than
// the true one.
TEST(crossed_cylinder_contacts_match_their_overlap) {
    float deepest = 0.0f;
    for (int k = 0; k < 2000; ++k) {
        World w;
        w.setGravity({});
        auto jitter = [k](int j) { return 0.004f * std::sin(12.9898f * (float)(k * 7 + j) + 78.233f); };
        BodyDesc a;
        a.position = { 1.1586f, 0.3697f, 0.5008f };
        a.rotation = normalize(Quat{ 0.880915f, 0.291875f, -0.359438f, -0.097993f });
        BodyDesc b;
        b.position = Vec3{ 1.3181f, 0.7676f, 1.1147f } + Vec3{ jitter(1), jitter(2), jitter(3) };
        b.rotation = normalize(Quat{ -0.879479f + jitter(4), -0.195779f + jitter(5), -0.218103f + jitter(6), 0.374991f });
        const BodyId ia = w.createBody(k % 2 ? b : a);
        const BodyId ib = w.createBody(k % 2 ? a : b);
        ShapeDesc sa;
        sa.geometry = Cylinder{ 0.361f, 0.237f };
        ShapeDesc sb;
        sb.geometry = Cylinder{ 0.343f, 0.298f };
        w.addShape(ia, k % 2 ? sb : sa);
        w.addShape(ib, k % 2 ? sa : sb);
        w.step(kDt);
        deepest = std::fmax(deepest, deepestContact(w));
    }
    CHECK(deepest < 0.02f);
}

TEST(ray_hits_cylinder_side_and_cap) {
    World w;
    BodyDesc bd;
    bd.type = MotionType::Static;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.5f, 1.0f };
    w.addShape(b, sd);
    RayHit hit;
    CHECK(w.castRay({ -5, 0, 0 }, { 10, 0, 0 }, hit));
    CHECK_NEAR(hit.fraction, 0.45f, 1e-5f);
    CHECK_NEAR(hit.normal.x, -1.0f, 1e-5f);
    CHECK(w.castRay({ -5, 0.499f, 0.3f }, { 10, 0, 0 }, hit));
    CHECK_NEAR(hit.point.x, -std::sqrt(0.25f - 0.499f * 0.499f), 1e-4f);
    CHECK(!w.castRay({ -5, 0.501f, 0.3f }, { 10, 0, 0 }, hit));
    CHECK(w.castRay({ 0.2f, 0.1f, 5 }, { 0, 0, -10 }, hit));
    CHECK_NEAR(hit.point.z, 1.0f, 1e-4f);
    CHECK_NEAR(hit.normal.z, 1.0f, 1e-5f);
    CHECK(!w.castRay({ 0, 0, 0 }, { 10, 0, 0 }, hit));
}

TEST(boxes_slide_on_slope_without_sinking) {
    float deepest = 0.0f;
    for (float friction : { 0.02f, 0.25f, 0.9f }) {
        World w;
        addGround(w);
        BodyDesc bd;
        bd.type = MotionType::Static;
        bd.position = { 0, 2, 0 };
        bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.35f);
        const BodyId slope = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Box{ { 5, 0.25f, 1.2f } };
        sd.friction = friction;
        w.addShape(slope, sd);

        BodyDesc boxDesc;
        boxDesc.rotation = bd.rotation;
        boxDesc.position = bd.position + bd.rotation.rotate({ 3.5f, 0.65f, 0 });
        const BodyId b = w.createBody(boxDesc);
        ShapeDesc bs;
        bs.geometry = Box{ { 0.4f, 0.4f, 0.4f } };
        bs.friction = friction;
        w.addShape(b, bs);
        for (int f = 0; f < 600; ++f) {
            w.step(kDt);
            deepest = std::fmax(deepest, deepestContact(w));
        }
    }
    CHECK(deepest < 0.005f);
}

TEST(sleeping_body_falls_when_its_floor_moves) {
    for (int mode = 0; mode < 3; ++mode) {
        World w;
        const MotionType type = mode == 0 ? MotionType::Static : MotionType::Kinematic;
        const BodyId floor = addBox(w, { 0, -0.5f, 0 }, { 3, 0.5f, 3 }, type);
        const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
        run(w, 3.0f);
        CHECK(!w.isAwake(b));
        if (mode < 2) {
            w.setTransform(floor, { 0, -2.5f, 0 }, {});
        } else {
            for (int i = 1; i <= 30; ++i) {
                w.setKinematicTarget(floor, { 0, -0.5f - 2.0f * (float)i / 30.0f, 0 }, {}, kDt);
                w.step(kDt);
            }
            w.setLinearVelocity(floor, {});
        }
        run(w, 1.5f);
        CHECK_NEAR(w.position(b).y, -1.5f, 0.02f);
    }
}

TEST(sleeping_body_wakes_when_a_body_is_moved_into_it) {
    for (int mode = 0; mode < 3; ++mode) {
        World w;
        addBox(w, { 0, -0.5f, 0 }, { 1, 0.5f, 1 }, MotionType::Static);
        addBox(w, { 10, 9.5f, 0 }, { 1, 0.5f, 1 }, MotionType::Static);
        const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
        const BodyId far = addBox(w, { 10, 10.5f, 0 }, { 0.5f, 0.5f, 0.5f });
        const BodyId mover = mode < 2 ? addBox(w, { 5, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f }, mode ? MotionType::Kinematic : MotionType::Static)
                                      : addGround(w);
        if (mode == 2) { w.setTransform(mover, { 0, -5, 0 }, kUp); }
        run(w, 3.0f);
        CHECK(!w.isAwake(b));
        CHECK(!w.isAwake(far));
        if (mode < 2) {
            w.setTransform(mover, { 0.1f, 0.5f, 0 }, {});
        } else {
            w.setTransform(mover, { 0, 0.4f, 0 }, kUp);
        }
        CHECK(w.isAwake(b));
        CHECK(!w.isAwake(far));
        run(w, 1.0f);
        CHECK(deepestContact(w) < 0.01f);
    }
}

TEST(resized_shape_settles_on_its_new_size) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.5f, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 0.5f, 0.5f } };
    const ShapeId s = w.addShape(b, sd);
    run(w, 1.0f);
    w.setGeometry(s, Box{ { 0.5f, 0.25f, 0.5f } });
    run(w, 1.0f);
    CHECK_NEAR(w.position(b).y, 0.25f, 0.02f);
}

TEST(static_body_made_dynamic_rests_on_the_static_body_below) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 4, 0.5f, 4 }, MotionType::Static);
    const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f }, MotionType::Static);
    run(w, 0.5f);
    w.setMotionType(b, MotionType::Dynamic);
    float lowest = 1.0f;
    for (int i = 0; i < 60; ++i) {
        w.step(kDt);
        lowest = std::fmin(lowest, w.position(b).y);
    }
    CHECK(lowest > 0.49f);
}

TEST(removing_a_body_leaves_others_on_its_floor_asleep) {
    World w;
    addGround(w);
    const BodyId a = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    const BodyId b = addBox(w, { 5, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 3.0f);
    CHECK(!w.isAwake(b));
    w.destroyBody(a);
    CHECK(!w.isAwake(b));
}

TEST(heavy_body_landing_on_a_light_one_does_not_crush_it_into_the_floor) {
    float worst = 0.0f;
    for (int trial = 0; trial < 8; ++trial) {
        World w;
        addBox(w, { 0, -0.5f, 0 }, { 4, 0.5f, 4 }, MotionType::Static);
        const BodyId light = addSphere(w, { 0, 0.235f, 0 }, 0.235f);
        const BodyId heavy = addBox(w, { 0.05f * (float)trial, 1.0f, 0.03f * (float)trial }, { 0.4f, 0.4f, 0.4f });
        w.setLinearVelocity(heavy, { 0, -6.0f - 1.0f * (float)trial, 0 });
        for (int i = 0; i < 30; ++i) {
            w.step(kDt);
            worst = std::fmax(worst, deepestContact(w));
        }
        CHECK(w.isValid(light));
    }
    CHECK(worst < 0.06f);
}

TEST(ball_keeps_rolling_without_rolling_resistance) {
    World w;
    addGround(w);
    const BodyId ball = addSphere(w, { 0, 0.5f, 0 }, 0.5f);
    w.setLinearVelocity(ball, { 2, 0, 0 });
    w.setAngularVelocity(ball, { 0, 0, -4 });
    run(w, 3.0f);
    CHECK_NEAR(w.linearVelocity(ball).x, 2.0f, 0.02f);
}

// Rolling resistance c takes c * g * r of work per unit of rolled distance, so a ball rolling at v stops
// after 1.4 v^2 / (2 c g) and a disc after 1.5 v^2 / (2 c g).
TEST(rolling_ball_stops_after_the_expected_distance) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.5f, 0 };
    bd.linearVelocity = { 2, 0, 0 };
    bd.angularVelocity = { 0, 0, -4 };
    const BodyId ball = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Sphere{ 0.5f };
    sd.rollingResistance = 0.1f;
    w.addShape(ball, sd);
    run(w, 4.0f);
    CHECK_NEAR(w.position(ball).x, 1.4f * 4.0f / (2.0f * 0.1f * 9.81f), 0.15f);
    CHECK(length(w.linearVelocity(ball)) < 0.01f);
    CHECK_NEAR(w.position(ball).y, 0.5f, 0.01f);
}

TEST(rolling_planar_disc_stops_after_the_expected_distance) {
    World w;
    addBox(w, { 0, -0.5f, 0 }, { 10, 0.5f, 1 }, MotionType::Static);
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 0.5f, 0 };
    bd.linearVelocity = { 2, 0, 0 };
    bd.angularVelocity = { 0, 0, -4 };
    const BodyId disc = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Cylinder{ 0.5f, 0.1f };
    sd.rollingResistance = 0.1f;
    w.addShape(disc, sd);
    run(w, 4.0f);
    CHECK_NEAR(w.position(disc).x, 1.5f * 4.0f / (2.0f * 0.1f * 9.81f), 0.15f);
    CHECK(length(w.linearVelocity(disc)) < 0.01f);
}

TEST(rolling_resistance_holds_a_ball_on_a_gentle_slope) {
    for (float slope : { 0.05f, 0.2f }) {
        World w;
        BodyDesc gd;
        gd.type = MotionType::Static;
        gd.rotation = Quat::fromTo({ 0, 0, 1 }, normalize(Vec3{ -slope, 1, 0 }));
        const BodyId g = w.createBody(gd);
        ShapeDesc ground;
        ground.geometry = HalfSpace{};
        w.addShape(g, ground);
        BodyDesc bd;
        bd.position = { 0, 0.5f, 0 };
        const BodyId ball = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Sphere{ 0.5f };
        sd.rollingResistance = 0.1f;
        w.addShape(ball, sd);
        run(w, 3.0f);
        if (slope < 0.1f) {
            CHECK(std::fabs(w.position(ball).x) < 0.05f);
        } else {
            CHECK(w.position(ball).x < -1.0f);
        }
    }
}

TEST(body_made_static_or_kinematic_drops_its_contact_with_the_ground) {
    for (MotionType type : { MotionType::Static, MotionType::Kinematic }) {
        World w;
        addGround(w);
        const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
        run(w, 0.5f);
        CHECK(w.contacts().size() == 1);
        w.setMotionType(b, type);
        CHECK(w.contacts(true).empty());
        run(w, 0.5f);
        CHECK(w.contacts(true).empty());
    }
}

TEST(destroying_sleeping_bodies_costs_what_they_touch) {
    // Pairs of stacked boxes, each its own island, destroyed one body at a time.
    auto seconds = [](int pairs) {
        double best = 1e9;
        for (int trial = 0; trial < 3; ++trial) {
            World w;
            addBox(w, { 0, -0.5f, 0 }, { 200, 0.5f, 200 }, MotionType::Static);
            std::vector<BodyId> ids;
            const int side = (int)std::ceil(std::sqrt((float)pairs));
            for (int i = 0; i < pairs; ++i) {
                const float x = 1.5f * (float)(i % side), z = 1.5f * (float)(i / side);
                ids.push_back(addBox(w, { x, 0.5f, z }, { 0.5f, 0.5f, 0.5f }));
                ids.push_back(addBox(w, { x, 1.5f, z }, { 0.5f, 0.5f, 0.5f }));
            }
            run(w, 1.0f);
            const auto start = std::chrono::steady_clock::now();
            for (BodyId b : ids) { w.destroyBody(b); }
            best = std::fmin(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        }
        return best;
    };
    CHECK(seconds(2000) < 8.0 * seconds(500));
}

TEST(moved_shape_settles_at_its_new_offset) {
    World w;
    addGround(w);
    const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 3.0f);
    CHECK(!w.isAwake(b));
    const ShapeId s = w.shapes(b)[0];
    w.setLocalTransform(s, { 0, -0.2f, 0 }, Quat::axisAngle({ 0, 1, 0 }, 0.5f));
    CHECK(w.isAwake(b));
    CHECK_NEAR(w.localTransform(s).p.y, -0.2f, 1e-6f);
    CHECK_NEAR(w.bounds(s).lower.y, w.position(b).y - 0.7f, 1e-4f);
    CHECK_NEAR(w.centerOfMass(b).y, w.position(b).y - 0.2f, 1e-4f);
    run(w, 1.0f);
    CHECK_NEAR(w.position(b).y, 0.7f, 0.01f);
    CHECK(deepestContact(w) < 0.01f);
}

TEST(moving_a_floor_shape_away_wakes_what_rests_on_it) {
    World w;
    const BodyId floor = addBox(w, { 0, -0.5f, 0 }, { 3, 0.5f, 3 }, MotionType::Static);
    const BodyId b = addBox(w, { 0, 0.5f, 0 }, { 0.5f, 0.5f, 0.5f });
    run(w, 3.0f);
    CHECK(!w.isAwake(b));
    w.setLocalTransform(w.shapes(floor)[0], { 0, -2, 0 }, {});
    CHECK(w.isAwake(b));
    run(w, 1.5f);
    CHECK_NEAR(w.position(b).y, -1.5f, 0.02f);
}

TEST(shape_turned_in_place_lands_on_its_new_side) {
    World w;
    addGround(w);
    const BodyId b = addBox(w, { 0, 0.8f, 0 }, { 0.5f, 0.8f, 0.5f });
    run(w, 3.0f);
    CHECK(!w.isAwake(b));
    w.setLocalTransform(w.shapes(b)[0], {}, Quat::axisAngle({ 1, 0, 0 }, 0.5f * kPi));
    run(w, 1.0f);
    CHECK_NEAR(w.position(b).y, 0.5f, 0.01f);
}

TEST(body_knocked_into_moving_paddle_ends_outside_it) {
    WorldSettings ws;
    ws.gravity = { 0, 0, 0 };
    World w(ws);
    const BodyId paddle = addBox(w, { 1.7f, 0, 0 }, { 0.1f, 1, 1 }, MotionType::Kinematic);
    w.setLinearVelocity(paddle, { -30, 0, 0 });
    const BodyId ball = addSphere(w, { 0, 0, 0 }, 0.25f, 1.0f);
    const BodyId striker = addSphere(w, { -0.51f, 0, 0 }, 0.25f, 1.0f);
    w.setLinearVelocity(striker, { 200, 0, 0 });
    for (int i = 0; i < 20; ++i) {
        w.step(kDt);
        CHECK(std::fabs(w.position(ball).x - w.position(paddle).x) > 0.34f);
    }
    CHECK(w.linearVelocity(ball).x < -100.0f);
}
