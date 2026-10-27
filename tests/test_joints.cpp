#include <cmath>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kPi = 3.14159265358979f;

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

void addGround(World& w) {
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = kUp;
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    w.addShape(g, sd);
}

BodyId addBody(World& w, const Vec3& pos, const Geometry& g, MotionType type = MotionType::Dynamic, uint8_t locks = 0) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    bd.locks = locks;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    w.addShape(b, sd);
    return b;
}

JointId join(World& w, BodyId a, BodyId b, const JointKind& kind, bool collide = false) {
    JointDesc jd;
    jd.bodyA = a;
    jd.bodyB = b;
    jd.kind = kind;
    jd.collideConnected = collide;
    return w.createJoint(jd);
}

float valueOf(const World& w, JointId j) {
    for (const JointInfo& info : w.joints()) {
        if (info.joint == j) { return info.value; }
    }
    return NAN;
}

Transform relative(const Transform& a, const Transform& b) { return { a.applyInverse(b.p), a.q.conjugate() * b.q }; }

float angleBetween(const Quat& a, const Quat& b) {
    return 2.0f * std::acos(std::fmin(std::fabs(dot(a, b)), 1.0f));
}

}

TEST(invalid_joints_are_refused) {
    World w;
    const BodyId a = addBody(w, {}, Box{});
    CHECK(!w.isValid(join(w, a, BodyId{}, PointJoint{})));
    CHECK(!w.isValid(join(w, a, a, PointJoint{})));
    const BodyId gone = addBody(w, { 3, 0, 0 }, Box{});
    w.destroyBody(gone);
    CHECK(!w.isValid(join(w, gone, a, PointJoint{})));
    CHECK(w.isValid(join(w, BodyId{}, a, PointJoint{})));
}

TEST(pendulum_swings_at_its_period) {
    World w;
    const float len = 1.0f;
    const float start = 0.1f;
    const BodyId bob = addBody(w, { len * std::sin(start), 5.0f - len * std::cos(start), 0 }, Sphere{ 0.02f });
    join(w, BodyId{}, bob, PointJoint{ { 0, 5, 0 } });
    float last = 0.0f, firstCross = -1.0f, period = 0.0f;
    int crossings = 0;
    float worstStretch = 0.0f;
    for (int i = 1; i <= 600; ++i) {
        w.step(kDt);
        const Vec3 p = w.position(bob);
        worstStretch = std::fmax(worstStretch, std::fabs(length(p - Vec3{ 0, 5, 0 }) - len));
        if (last > 0.0f && p.x <= 0.0f) {
            const float t = (float)i * kDt - kDt * p.x / (p.x - last);
            if (firstCross < 0.0f) {
                firstCross = t;
            } else {
                period = (t - firstCross) / (float)crossings;
            }
            ++crossings;
        }
        last = p.x;
    }
    const float expected = 2.0f * kPi * std::sqrt(len / 9.81f);
    CHECK(crossings >= 4);
    CHECK_NEAR(period, expected, 0.02f * expected);
    CHECK(worstStretch < 0.002f);
}

TEST(hanging_chain_holds_its_length) {
    for (int swing = 0; swing < 2; ++swing) {
        World w;
        const int links = 12;
        BodyId prev;
        BodyId last;
        for (int i = 0; i < links; ++i) {
            const float y = 10.0f - 0.5f * (float)i - 0.25f;
            const bool end = i == links - 1;
            last = addBody(w, { 0, y, 0 }, end ? Geometry{ Sphere{ 0.25f } } : Geometry{ Box{ { 0.05f, 0.25f, 0.05f } } });
            join(w, prev, last, PointJoint{ { 0, y + 0.25f, 0 } });
            prev = last;
        }
        if (swing == 1) { w.setLinearVelocity(last, { 4, 0, 0 }); }
        float worstGap = 0.0f;
        for (int i = 0; i < 600; ++i) {
            w.step(kDt);
            for (const JointInfo& info : w.joints()) { worstGap = std::fmax(worstGap, length(info.anchorB - info.anchorA)); }
        }
        CHECK(worstGap < 0.01f);
        if (swing == 0) { CHECK_NEAR(10.0f - w.position(last).y, 0.5f * (float)links - 0.25f, 0.03f); }
    }
}

TEST(fixed_joint_keeps_bodies_in_their_joined_pose) {
    World w;
    addGround(w);
    const BodyId a = addBody(w, { 0, 3, 0 }, Box{ { 0.5f, 0.2f, 0.2f } });
    const BodyId b = addBody(w, { 0.7f, 3.3f, 0 }, Box{ { 0.2f, 0.5f, 0.2f } });
    join(w, a, b, FixedJoint{ { 0.6f, 3.0f, 0 } });
    w.setAngularVelocity(a, { 2, 5, 1 });
    const Transform rest = relative(w.transform(a), w.transform(b));
    float worstMove = 0.0f, worstTurn = 0.0f;
    for (int i = 0; i < 240; ++i) {
        w.step(kDt);
        const Transform rel = relative(w.transform(a), w.transform(b));
        worstMove = std::fmax(worstMove, length(rel.p - rest.p));
        worstTurn = std::fmax(worstTurn, angleBetween(rel.q, rest.q));
    }
    CHECK(worstMove < 0.01f);
    CHECK(worstTurn < 0.02f);
    CHECK(w.position(a).y < 1.0f);
}

TEST(joined_bodies_do_not_collide_unless_asked) {
    for (int collide = 0; collide < 2; ++collide) {
        World w;
        const BodyId a = addBody(w, { 0, 3, 0 }, Box{});
        const BodyId b = addBody(w, { 0.5f, 3, 0 }, Box{});
        w.setGravity({});
        join(w, a, b, PointJoint{ { 0.25f, 3, 0 } }, collide == 1);
        for (int i = 0; i < 30; ++i) { w.step(kDt); }
        CHECK((w.contacts().size() > 0) == (collide == 1));
    }
}

TEST(hinge_turns_only_about_its_axis) {
    World w;
    w.setGravity({});
    const BodyId door = addBody(w, { 0.5f, 1, 0 }, Box{ { 0.5f, 1.0f, 0.05f } });
    join(w, BodyId{}, door, HingeJoint{ { 0, 1, 0 }, { 0, 1, 0 } });
    w.setAngularVelocity(door, { 3, 2, -4 });
    float worst = 0.0f;
    for (int i = 0; i < 120; ++i) {
        w.step(kDt);
        const Vec3 axis = w.rotation(door).rotate({ 0, 1, 0 });
        worst = std::fmax(worst, length(cross(axis, { 0, 1, 0 })));
    }
    const Vec3 spin = w.angularVelocity(door);
    CHECK(worst < 0.01f);
    CHECK(std::fabs(spin.y) > 0.3f);
    CHECK(length(cross(spin, { 0, 1, 0 })) < 0.01f);
    CHECK(length(w.centerOfMass(door) - Vec3{ 0, 1, 0 }) < 0.51f);
}

TEST(hinge_limits_hold_the_swing) {
    World w;
    const BodyId bar = addBody(w, { 1, 5, 0 }, Box{ { 1.0f, 0.1f, 0.1f } });
    HingeJoint hj{ { 0, 5, 0 }, { 0, 0, 1 } };
    hj.limited = true;
    hj.lowerAngle = -0.5f;
    hj.upperAngle = 0.3f;
    const JointId j = join(w, BodyId{}, bar, hj);
    w.setLinearVelocity(bar, { 0, 6, 0 });
    w.setAngularVelocity(bar, { 0, 0, 6 });
    float lowest = 0.0f, highest = 0.0f;
    for (int i = 0; i < 240; ++i) {
        w.step(kDt);
        lowest = std::fmin(lowest, valueOf(w, j));
        highest = std::fmax(highest, valueOf(w, j));
    }
    CHECK(highest > 0.25f && highest < 0.33f);
    CHECK(lowest > -0.53f);
    CHECK_NEAR(valueOf(w, j), -0.5f, 0.02f);
}

TEST(hinge_limit_stops_the_swing_without_a_bounce) {
    World w;
    w.setGravity({});
    BodyDesc bd;
    bd.position = { 0.6f, 0, 0 };
    bd.angularVelocity = { 0, 3, 0 };
    bd.linearVelocity = { 0, 0, -1.8f };
    const BodyId door = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.6f, 1.0f, 0.06f } };
    w.addShape(door, sd);
    HingeJoint hj{ {}, { 0, 1, 0 } };
    hj.limited = true;
    hj.lowerAngle = -1.0f;
    hj.upperAngle = 0.5f;
    const JointId j = join(w, BodyId{}, door, hj);
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    CHECK_NEAR(valueOf(w, j), 0.5f, 0.002f);
    CHECK(length(w.angularVelocity(door)) < 0.01f);
}

TEST(slider_limit_stops_an_offset_body_without_a_bounce) {
    World w;
    w.setGravity({});
    BodyDesc bd;
    bd.linearVelocity = { 3, 0, 0 };
    const BodyId block = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.3f, 0.3f, 0.3f } };
    w.addShape(block, sd);
    SliderJoint sj{ { 0, 2, 0 }, { 1, 0, 0 } };
    sj.limited = true;
    sj.lowerTranslation = -1.0f;
    sj.upperTranslation = 0.5f;
    const JointId j = join(w, BodyId{}, block, sj);
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    CHECK_NEAR(valueOf(w, j), 0.5f, 0.002f);
    CHECK(length(w.linearVelocity(block)) < 0.01f);
}

TEST(hinge_motor_reaches_its_speed) {
    World w;
    const BodyId wheel = addBody(w, { 0, 2, 0 }, Cylinder{ 0.5f, 0.1f });
    HingeJoint hj{ { 0, 2, 0 }, { 0, 0, 1 } };
    hj.motorSpeed = 3.0f;
    hj.maxMotorTorque = 50.0f;
    join(w, BodyId{}, wheel, hj);
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    CHECK_NEAR(w.angularVelocity(wheel).z, 3.0f, 0.01f);
}

TEST(hinge_friction_runs_a_pendulum_down) {
    float swing[2] = {};
    for (int k = 0; k < 2; ++k) {
        World w;
        const BodyId bar = addBody(w, { 1, 5, 0 }, Box{ { 1.0f, 0.1f, 0.1f } });
        HingeJoint hj{ { 0, 5, 0 }, { 0, 0, 1 } };
        hj.maxMotorTorque = k == 1 ? 3.0f : 0.0f;
        join(w, BodyId{}, bar, hj);
        for (int i = 0; i < 300; ++i) { w.step(kDt); }
        for (int i = 0; i < 120; ++i) {
            w.step(kDt);
            swing[k] = std::fmax(swing[k], std::fabs(w.angularVelocity(bar).z));
        }
    }
    CHECK(swing[0] > 2.0f);
    CHECK(swing[1] < 0.05f);
}

TEST(slider_moves_only_along_its_axis_and_stops_at_its_limit) {
    World w;
    const Vec3 axis = normalize(Vec3{ 2, -1, 0.5f });
    const BodyId block = addBody(w, { 0, 5, 0 }, Box{ { 0.3f, 0.2f, 0.4f } });
    SliderJoint sj{ { 0, 5, 0 }, axis };
    sj.limited = true;
    sj.lowerTranslation = -0.5f;
    sj.upperTranslation = 2.0f;
    const JointId j = join(w, BodyId{}, block, sj);
    w.setAngularVelocity(block, { 1, 2, 3 });
    float worstSide = 0.0f, worstTurn = 0.0f, furthest = 0.0f;
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        const Vec3 d = w.position(block) - Vec3{ 0, 5, 0 };
        worstSide = std::fmax(worstSide, length(d - axis * dot(d, axis)));
        worstTurn = std::fmax(worstTurn, angleBetween(w.rotation(block), Quat{}));
        furthest = std::fmax(furthest, valueOf(w, j));
    }
    CHECK(worstSide < 0.005f);
    CHECK(worstTurn < 0.02f);
    CHECK(furthest < 2.03f);
    CHECK_NEAR(valueOf(w, j), 2.0f, 0.01f);
}

TEST(slider_motor_lifts_against_gravity) {
    World w;
    const BodyId lift = addBody(w, { 0, 1, 0 }, Box{ { 0.5f, 0.1f, 0.5f } });
    SliderJoint sj{ { 0, 1, 0 }, { 0, 1, 0 } };
    sj.motorSpeed = 0.5f;
    sj.maxMotorForce = 200.0f;
    const JointId j = join(w, BodyId{}, lift, sj);
    for (int i = 0; i < 120; ++i) { w.step(kDt); }
    CHECK_NEAR(w.linearVelocity(lift).y, 0.5f, 0.01f);
    CHECK_NEAR(valueOf(w, j), 1.0f, 0.03f);
}

TEST(distance_joint_holds_its_length_or_range) {
    World w;
    const BodyId rod = addBody(w, { 2, 5, 0 }, Sphere{ 0.1f });
    const JointId held = join(w, BodyId{}, rod, DistanceJoint{ { 0, 5, 0 }, { 2, 5, 0 } });
    const BodyId rope = addBody(w, { 2, 5, 3 }, Sphere{ 0.1f });
    DistanceJoint dj{ { 0, 5, 3 }, { 2, 5, 3 } };
    dj.limited = true;
    dj.minDistance = 0.0f;
    dj.maxDistance = 3.0f;
    const JointId slack = join(w, BodyId{}, rope, dj);
    float worstHeld = 0.0f, longest = 0.0f;
    for (int i = 0; i < 240; ++i) {
        w.step(kDt);
        worstHeld = std::fmax(worstHeld, std::fabs(valueOf(w, held) - 2.0f));
        longest = std::fmax(longest, valueOf(w, slack));
    }
    CHECK(worstHeld < 0.005f);
    CHECK(longest > 2.95f && longest < 3.03f);
}

TEST(distance_spring_oscillates_at_its_frequency) {
    World w;
    w.setGravity({});
    const BodyId bob = addBody(w, { 1.5f, 0, 0 }, Sphere{ 0.1f });
    DistanceJoint dj{ {}, { 1.5f, 0, 0 } };
    dj.length = 1.0f;
    dj.hertz = 1.0f;
    const JointId j = join(w, BodyId{}, bob, dj);
    float last = valueOf(w, j) - 1.0f, first = -1.0f, period = 0.0f;
    int crossings = 0;
    for (int i = 1; i <= 360; ++i) {
        w.step(kDt);
        const float x = valueOf(w, j) - 1.0f;
        if (last > 0.0f && x <= 0.0f) {
            const float t = (float)i * kDt;
            if (first < 0.0f) {
                first = t;
            } else {
                period = (t - first) / (float)crossings;
            }
            ++crossings;
        }
        last = x;
    }
    CHECK(crossings >= 4);
    CHECK_NEAR(period, 1.0f, 0.05f);
}

TEST(planar_chain_stays_in_its_plane) {
    World w;
    addGround(w);
    BodyId prev;
    BodyId last;
    for (int i = 0; i < 8; ++i) {
        const float x = 0.5f * (float)i + 0.25f;
        last = addBody(w, { x, 6, 0 }, Box{ { 0.25f, 0.05f, 0.05f } }, MotionType::Dynamic, LockPlanar);
        join(w, prev, last, HingeJoint{ { x - 0.25f, 6, 0 }, { 0, 0, 1 } });
        prev = last;
    }
    float worstZ = 0.0f, worstGap = 0.0f;
    for (int i = 0; i < 300; ++i) {
        w.step(kDt);
        for (const JointInfo& info : w.joints()) {
            worstGap = std::fmax(worstGap, length(info.anchorB - info.anchorA));
        }
        worstZ = std::fmax(worstZ, std::fabs(w.position(last).z));
    }
    CHECK(worstZ < 1e-4f);
    CHECK(worstGap < 0.01f);
}

TEST(kinematic_body_carries_what_hangs_from_it) {
    World w;
    const BodyId crane = addBody(w, { 0, 5, 0 }, Box{ { 0.2f, 0.2f, 0.2f } }, MotionType::Kinematic);
    const BodyId load = addBody(w, { 0, 3, 0 }, Sphere{ 0.3f });
    join(w, load, crane, DistanceJoint{ { 0, 3, 0 }, { 0, 5, 0 } });
    for (int i = 0; i < 240; ++i) { w.step(kDt); }
    CHECK(!w.isAwake(load));
    w.setLinearVelocity(crane, { 1, 0, 0 });
    for (int i = 0; i < 120; ++i) { w.step(kDt); }
    CHECK(w.isAwake(load));
    CHECK(w.position(load).x > 1.0f);
}

TEST(joined_bodies_sleep_and_wake_together) {
    World w;
    addGround(w);
    const BodyId a = addBody(w, { 0, 0.5f, 0 }, Box{});
    const BodyId b = addBody(w, { 3, 0.5f, 0 }, Box{});
    join(w, a, b, DistanceJoint{ { 0, 0.5f, 0 }, { 3, 0.5f, 0 } });
    for (int i = 0; i < 120; ++i) { w.step(kDt); }
    CHECK(!w.isAwake(a) && !w.isAwake(b));
    CHECK(w.island(a) == w.island(b));
    w.wake(a);
    CHECK(w.isAwake(b));
}

TEST(destroying_a_body_destroys_its_joints) {
    World w;
    const BodyId a = addBody(w, { 0, 5, 0 }, Box{});
    const BodyId b = addBody(w, { 2, 5, 0 }, Box{});
    const JointId j = join(w, a, b, PointJoint{ { 1, 5, 0 } });
    w.destroyBody(b);
    CHECK(!w.isValid(j));
    CHECK(w.joints().empty());
}

TEST(copied_world_with_joints_steps_identically) {
    World w;
    addGround(w);
    BodyId prev;
    for (int i = 0; i < 6; ++i) {
        const BodyId link = addBody(w, { 0.6f * (float)i, 4, 0 }, Box{ { 0.3f, 0.1f, 0.1f } });
        join(w, prev, link, HingeJoint{ { 0.6f * (float)i - 0.3f, 4, 0 }, { 0, 1, 1 } });
        prev = link;
    }
    for (int i = 0; i < 30; ++i) { w.step(kDt); }
    World copy = w;
    for (int i = 0; i < 60; ++i) {
        w.step(kDt);
        copy.step(kDt);
    }
    const std::vector<BodyId> ids = w.bodies();
    for (BodyId id : ids) {
        const Transform x = w.transform(id), y = copy.transform(id);
        CHECK(x.p == y.p && x.q.x == y.q.x && x.q.y == y.q.y && x.q.z == y.q.z && x.q.w == y.q.w);
    }
}
