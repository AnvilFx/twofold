#include <cmath>
#include <memory>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

ShapeId addGround(World& w) {
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = kUp;
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    return w.addShape(w.createBody(bd), sd);
}

ShapeId addShape(World& w, const BodyDesc& bd, const Geometry& g, bool sensor = false) {
    ShapeDesc sd;
    sd.geometry = g;
    sd.sensor = sensor;
    return w.addShape(w.createBody(bd), sd);
}

ShapeId addShape(World& w, MotionType type, const Vec3& pos, const Geometry& g, bool sensor = false) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    return addShape(w, bd, g, sensor);
}

struct Tally {
    int begins = 0;
    int ends = 0;
};

void count(const World& w, ShapeId a, ShapeId b, Tally& t) {
    for (const ContactEvent& e : w.contactEvents()) {
        if (!((e.shapeA == a && e.shapeB == b) || (e.shapeA == b && e.shapeB == a))) { continue; }
        if (e.begin) {
            ++t.begins;
        } else {
            ++t.ends;
        }
    }
}

void countSensor(const World& w, ShapeId sensor, ShapeId visitor, Tally& t) {
    for (const SensorEvent& e : w.sensorEvents()) {
        if (e.sensor != sensor || e.visitor != visitor) { continue; }
        if (e.begin) {
            ++t.begins;
        } else {
            ++t.ends;
        }
    }
}

bool same(const ContactEvent& a, const ContactEvent& b) {
    return a.bodyA == b.bodyA && a.bodyB == b.bodyB && a.shapeA == b.shapeA && a.shapeB == b.shapeB && a.begin == b.begin;
}

bool same(const HitEvent& a, const HitEvent& b) {
    return a.shapeA == b.shapeA && a.shapeB == b.shapeB && a.point == b.point && a.normal == b.normal &&
           a.speed == b.speed;
}

bool same(const SensorEvent& a, const SensorEvent& b) {
    return a.sensor == b.sensor && a.visitor == b.visitor && a.sensorBody == b.sensorBody &&
           a.visitorBody == b.visitorBody && a.begin == b.begin;
}

template <class T>
bool sameList(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) { return false; }
    for (size_t i = 0; i < a.size(); ++i) {
        if (!same(a[i], b[i])) { return false; }
    }
    return true;
}

}

TEST(box_landing_and_lifting_off_begins_and_ends_once) {
    World w;
    const ShapeId ground = addGround(w);
    const ShapeId box = addShape(w, MotionType::Dynamic, { 0, 2, 0 }, Box{});
    const BodyId body = w.shapeBody(box);
    Tally t;
    for (int i = 0; i < 120; ++i) {
        w.step(kDt);
        for (const ContactEvent& e : w.contactEvents()) { CHECK(e.shapeA == ground && e.shapeB == box && e.bodyB == body); }
        count(w, ground, box, t);
    }
    CHECK(t.begins == 1 && t.ends == 0);
    w.setLinearVelocity(body, { 0, 4, 0 });
    for (int i = 0; i < 30; ++i) {
        w.step(kDt);
        count(w, ground, box, t);
    }
    CHECK(t.begins == 1 && t.ends == 1);
}

TEST(hit_reports_the_landing_speed) {
    World w;
    const ShapeId ground = addGround(w);
    const ShapeId ball = addShape(w, MotionType::Dynamic, { 0, 5.5f, 0 }, Sphere{ 0.5f });
    int hits = 0;
    HitEvent first;
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        for (const HitEvent& h : w.hitEvents()) {
            if (hits++ == 0) { first = h; }
        }
    }
    CHECK(hits == 1);
    CHECK(first.shapeA == ground && first.shapeB == ball);
    CHECK_NEAR(first.speed, std::sqrt(2.0f * 9.81f * 5.0f), 0.3f);
    CHECK(first.normal.y > 0.99f);
    CHECK_NEAR(first.point.y, 0.0f, 0.1f);
}

TEST(speculative_contact_without_touch_reports_nothing) {
    World w;
    w.setGravity({});
    addShape(w, MotionType::Static, {}, Box{ { 5, 1, 0.5f } });
    BodyDesc bd;
    bd.position = { -4, 0, 1.01f };
    bd.linearVelocity = { 2, 0, 0 };
    addShape(w, bd, Box{});
    bool speculative = false;
    for (int i = 0; i < 240; ++i) {
        w.step(kDt);
        speculative = speculative || !w.contacts(true).empty();
        CHECK(w.contactEvents().empty());
        CHECK(w.hitEvents().empty());
    }
    CHECK(speculative);
}

TEST(destroying_a_body_ends_its_touches_and_overlaps) {
    World w;
    const ShapeId ground = addGround(w);
    const ShapeId zone = addShape(w, MotionType::Static, { 0, 1, 0 }, Box{ { 2, 2, 2 } }, true);
    const ShapeId box = addShape(w, MotionType::Dynamic, { 0, 0.5f, 0 }, Box{});
    const BodyId body = w.shapeBody(box);
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    w.destroyBody(body);
    const ShapeId again = addShape(w, MotionType::Dynamic, { 0, 0.5f, 0 }, Box{});
    CHECK(again.index == box.index && again != box);
    w.step(kDt);

    const std::vector<ContactEvent>& e = w.contactEvents();
    CHECK(e.size() == 2);
    CHECK(!e[0].begin && e[0].shapeA == ground && e[0].shapeB == box && e[0].bodyB == body);
    CHECK(e[1].begin && e[1].shapeA == ground && e[1].shapeB == again);

    const std::vector<SensorEvent>& s = w.sensorEvents();
    CHECK(s.size() == 2);
    CHECK(!s[0].begin && s[0].sensor == zone && s[0].visitor == box && s[0].visitorBody == body);
    CHECK(s[1].begin && s[1].sensor == zone && s[1].visitor == again);
}

TEST(filter_change_ends_only_the_touches_it_separates) {
    World w;
    const ShapeId ground = addGround(w);
    const ShapeId box = addShape(w, MotionType::Dynamic, { 0, 0.5f, 0 }, Box{});
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    w.wake(w.shapeBody(box));
    w.setFilter(box, 1, 0xffffffffu);
    w.step(kDt);
    CHECK(w.contactEvents().empty());

    w.setFilter(box, 1, 0);
    Tally t;
    for (int i = 0; i < 30; ++i) {
        w.step(kDt);
        count(w, ground, box, t);
    }
    CHECK(t.begins == 0 && t.ends == 1);
}

TEST(contacts_dropped_by_a_motion_type_change_end) {
    World w;
    const ShapeId ground = addGround(w);
    const ShapeId box = addShape(w, MotionType::Dynamic, { 0, 0.5f, 0 }, Box{});
    for (int i = 0; i < 60; ++i) { w.step(kDt); }
    w.setMotionType(w.shapeBody(box), MotionType::Kinematic);
    Tally t;
    for (int i = 0; i < 10; ++i) {
        w.step(kDt);
        count(w, ground, box, t);
    }
    CHECK(t.begins == 0 && t.ends == 1);
}

TEST(sleeping_pairs_report_nothing) {
    World w;
    addGround(w);
    BodyId top;
    for (int i = 0; i < 3; ++i) { top = w.shapeBody(addShape(w, MotionType::Dynamic, { 0, 0.5f + (float)i, 0 }, Box{})); }
    for (int i = 0; i < 180; ++i) { w.step(kDt); }
    for (BodyId b : w.bodies()) { CHECK(!w.isAwake(b) || w.motionType(b) != MotionType::Dynamic); }
    for (int i = 0; i < 120; ++i) {
        w.step(kDt);
        CHECK(w.contactEvents().empty() && w.hitEvents().empty());
    }
    w.wake(top);
    for (int i = 0; i < 30; ++i) {
        w.step(kDt);
        CHECK(w.contactEvents().empty() && w.hitEvents().empty());
    }
}

TEST(mesh_pair_reports_once_per_shape_pair) {
    World w;
    const float s = std::sin(0.4f), c = std::cos(0.4f);
    const std::vector<Vec3> v = { { -3 * c, 3 * s, -2 }, { -3 * c, 3 * s, 2 }, { 0, 0, -2 },
                                  { 0, 0, 2 },           { 3 * c, 3 * s, -2 }, { 3 * c, 3 * s, 2 } };
    const auto mesh = std::make_shared<TriangleMesh>(v, std::vector<uint32_t>{ 0, 1, 2, 2, 1, 3, 2, 3, 4, 4, 3, 5 });
    const ShapeId valley = addShape(w, MotionType::Static, {}, Mesh{ mesh });
    const ShapeId box = addShape(w, MotionType::Dynamic, { 0, 2, 0 }, Box{ { 0.3f, 0.3f, 0.3f } });
    Tally t;
    for (int i = 0; i < 300; ++i) {
        w.step(kDt);
        count(w, valley, box, t);
    }
    int parts = 0;
    for (const ContactInfo& ci : w.contacts()) {
        if (ci.shapeA == valley && ci.shapeB == box) { ++parts; }
    }
    CHECK(parts == 2);
    CHECK(t.begins == 1 && t.ends == 0);
}

TEST(planar_pairs_report_touches) {
    World w;
    const ShapeId floor = addShape(w, MotionType::Static, { 0, -0.5f, 0 }, Box{ { 5, 0.5f, 1 } });
    BodyDesc bd;
    bd.locks = LockPlanar;
    bd.position = { 0, 1, 0 };
    const ShapeId low = addShape(w, bd, Box{});
    bd.position = { 0.2f, 3, 0 };
    const ShapeId high = addShape(w, bd, Box{});
    Tally lowFloor, lowHigh;
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        count(w, floor, low, lowFloor);
        count(w, low, high, lowHigh);
    }
    CHECK(lowFloor.begins == 1 && lowFloor.ends == 0);
    CHECK(lowHigh.begins == 1 && lowHigh.ends == 0);
}

TEST(static_sensor_reports_a_body_passing_through) {
    World w;
    w.setGravity({});
    const ShapeId zone = addShape(w, MotionType::Static, {}, Box{ { 1, 1, 1 } }, true);
    BodyDesc bd;
    bd.position = { -5, 0, 0 };
    bd.linearVelocity = { 5, 0, 0 };
    const ShapeId ball = addShape(w, bd, Sphere{ 0.5f });
    const BodyId body = w.shapeBody(ball);
    int entered = -1, left = -1;
    Tally t;
    for (int i = 0; i < 120; ++i) {
        w.step(kDt);
        for (const SensorEvent& e : w.sensorEvents()) {
            CHECK(e.sensor == zone && e.visitor == ball && e.sensorBody == w.shapeBody(zone) && e.visitorBody == body);
            (e.begin ? entered : left) = i + 1;
        }
        countSensor(w, zone, ball, t);
        CHECK(w.contactEvents().empty() && w.contacts(true).empty());
    }
    CHECK(t.begins == 1 && t.ends == 1);
    CHECK(entered >= 41 && entered <= 44);
    CHECK(left >= 77 && left <= 80);
    CHECK(w.linearVelocity(body) == Vec3(5, 0, 0));
    CHECK(w.isSensor(zone) && !w.isSensor(ball));
}

TEST(kinematic_sensor_sees_sleeping_bodies_without_waking_them) {
    World w;
    const ShapeId ground = addGround(w);
    BodyDesc bd;
    bd.position = { 0, 0.5f, 0 };
    bd.awake = false;
    const ShapeId box = addShape(w, bd, Box{});
    const BodyId body = w.shapeBody(box);
    BodyDesc kd;
    kd.type = MotionType::Kinematic;
    kd.position = { -5, 0.5f, 0 };
    kd.linearVelocity = { 4, 0, 0 };
    const ShapeId sweeper = addShape(w, kd, Box{ { 0.5f, 1, 0.5f } }, true);
    Tally seen, floor;
    for (int i = 0; i < 150; ++i) {
        w.step(kDt);
        countSensor(w, sweeper, box, seen);
        countSensor(w, sweeper, ground, floor);
    }
    CHECK(seen.begins == 1 && seen.ends == 1);
    CHECK(floor.begins == 1 && floor.ends == 0);
    CHECK(!w.isAwake(body));
    CHECK(w.position(body) == Vec3(0, 0.5f, 0));
}

TEST(sensors_follow_filters_and_ignore_other_sensors) {
    World w;
    w.setGravity({});
    BodyDesc sd;
    sd.type = MotionType::Static;
    ShapeDesc zoneDesc;
    zoneDesc.geometry = Box{ { 1, 1, 1 } };
    zoneDesc.sensor = true;
    zoneDesc.mask = 1;
    const ShapeId zone = w.addShape(w.createBody(sd), zoneDesc);
    const ShapeId other = addShape(w, MotionType::Static, {}, Box{ { 1, 1, 1 } }, true);
    BodyDesc bd;
    bd.position = { -0.5f, 0, 0 };
    const ShapeId seen = addShape(w, bd, Sphere{ 0.2f });
    bd.position = { 0.5f, 0, 0 };
    ShapeDesc hidden;
    hidden.geometry = Sphere{ 0.2f };
    hidden.category = 2;
    const ShapeId unseen = w.addShape(w.createBody(bd), hidden);
    w.step(kDt);
    Tally a, b, c;
    countSensor(w, zone, seen, a);
    countSensor(w, zone, unseen, b);
    countSensor(w, zone, other, c);
    CHECK(a.begins == 1 && b.begins == 0 && c.begins == 0);
    for (const SensorEvent& e : w.sensorEvents()) { CHECK(e.visitor != zone && e.visitor != other); }

    w.setFilter(seen, 1, 2);
    w.step(kDt);
    a = {};
    countSensor(w, zone, seen, a);
    CHECK(a.ends == 1);
}

TEST(queries_skip_sensors) {
    World w;
    w.setGravity({});
    addShape(w, MotionType::Static, {}, Box{ { 1, 1, 1 } }, true);
    const ShapeId box = addShape(w, MotionType::Static, { 3, 0, 0 }, Box{ { 0.5f, 0.5f, 0.5f } });
    RayHit hit;
    CHECK(w.castRay({ -5, 0, 0 }, { 10, 0, 0 }, hit) && hit.shape == box);
    std::vector<RayHit> hits;
    CHECK(w.castRayAll({ -5, 0, 0 }, { 10, 0, 0 }, hits) && hits.size() == 1 && hits[0].shape == box);
    CastHit cast;
    CHECK(w.castGeometry(Sphere{ 0.2f }, Transform{ { -5, 0, 0 }, {} }, { 10, 0, 0 }, cast) && cast.shape == box);
    std::vector<ShapeId> found;
    CHECK(!w.overlapGeometry(Sphere{ 0.2f }, Transform{ {}, {} }, found) && found.empty());
    CHECK(!w.overlapAabb(Aabb{ { -0.5f, -0.5f, -0.5f }, { 0.5f, 0.5f, 0.5f } }, found) && found.empty());
    CHECK(!w.overlapPoint({}, found) && found.empty());
}

TEST(sensors_do_not_change_how_bodies_move) {
    World plain, sensed;
    for (World* w : { &plain, &sensed }) {
        addGround(*w);
        for (int i = 0; i < 4; ++i) {
            BodyDesc bd;
            bd.position = { 0.3f * (float)i, 1.0f + 1.2f * (float)i, 0.1f * (float)i };
            bd.rotation = Quat::axisAngle({ 1, 0, 1 }, 0.3f * (float)i);
            bd.angularVelocity = { 0, 2, 0 };
            addShape(*w, bd, i % 2 ? Geometry{ Sphere{ 0.4f } } : Geometry{ Box{ { 0.4f, 0.3f, 0.5f } } });
        }
    }
    addShape(sensed, MotionType::Static, { 0, 2, 0 }, Box{ { 3, 1, 3 } }, true);
    ShapeDesc ring;
    ring.geometry = Sphere{ 3.0f };
    ring.sensor = true;
    ring.density = 5.0f;
    const BodyId carrier = sensed.bodies()[1];
    sensed.addShape(carrier, ring);
    CHECK(sensed.massProperties(carrier).mass == plain.massProperties(carrier).mass);

    int events = 0;
    for (int i = 0; i < 240; ++i) {
        plain.step(kDt);
        sensed.step(kDt);
        events += (int)sensed.sensorEvents().size();
    }
    CHECK(events > 0);
    for (BodyId b : plain.bodies()) {
        const Transform x = plain.transform(b), y = sensed.transform(b);
        CHECK(x.p.x == y.p.x && x.p.y == y.p.y && x.p.z == y.p.z);
        CHECK(x.q.x == y.q.x && x.q.y == y.q.y && x.q.z == y.q.z && x.q.w == y.q.w);
        CHECK(plain.isAwake(b) == sensed.isAwake(b));
    }

    RayHit hit;
    CHECK(sensed.castRay({ 2.5f, 5, 2.5f }, { 0, -10, 0 }, hit));
    CHECK(!sensed.isSensor(hit.shape));
    CHECK_NEAR(hit.point.y, 0.0f, 1e-4f);
}

TEST(events_are_ordered_and_copy_identically) {
    World w;
    addGround(w);
    addShape(w, MotionType::Static, { 0, 1.5f, 0 }, Box{ { 4, 0.5f, 4 } }, true);
    for (int i = 0; i < 6; ++i) {
        BodyDesc bd;
        bd.position = { 0.9f * (float)(i % 3) - 0.9f, 2.0f + 0.7f * (float)i, 0.5f * (float)(i / 3) };
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = i % 2 ? Geometry{ Sphere{ 0.3f } } : Geometry{ Box{ { 0.3f, 0.3f, 0.3f } } };
        sd.restitution = 0.6f;
        w.addShape(b, sd);
    }
    for (int i = 0; i < 20; ++i) { w.step(kDt); }
    World copy = w;
    int begins = 0, ends = 0, hits = 0, overlaps = 0;
    for (int i = 0; i < 180; ++i) {
        w.step(kDt);
        copy.step(kDt);
        CHECK(sameList(w.contactEvents(), copy.contactEvents()));
        CHECK(sameList(w.hitEvents(), copy.hitEvents()));
        CHECK(sameList(w.sensorEvents(), copy.sensorEvents()));
        const std::vector<ContactEvent>& e = w.contactEvents();
        for (size_t k = 0; k < e.size(); ++k) {
            CHECK(e[k].shapeA.index < e[k].shapeB.index);
            if (k > 0) {
                const bool after = e[k - 1].shapeA.index != e[k].shapeA.index ? e[k - 1].shapeA.index < e[k].shapeA.index
                                                                             : e[k - 1].shapeB.index < e[k].shapeB.index;
                CHECK(after);
            }
            if (e[k].begin) {
                ++begins;
            } else {
                ++ends;
            }
        }
        const std::vector<SensorEvent>& s = w.sensorEvents();
        for (size_t k = 1; k < s.size(); ++k) {
            CHECK(s[k - 1].sensor.index < s[k].sensor.index ||
                  (s[k - 1].sensor.index == s[k].sensor.index && s[k - 1].visitor.index < s[k].visitor.index));
        }
        hits += (int)w.hitEvents().size();
        overlaps += (int)s.size();
    }
    CHECK(begins >= 6 && ends > 0 && hits >= 6 && overlaps > 6);
}
