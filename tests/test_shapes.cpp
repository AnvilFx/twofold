#include <cmath>
#include <cstdint>
#include <limits>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kNan = std::numeric_limits<float>::quiet_NaN();
constexpr float kInf = std::numeric_limits<float>::infinity();

struct Drop {
    float height = 0.0f;
    float deepest = 0.0f;
    bool touching = false;
    bool unitNormals = true;
    bool finite = true;
};

// Drops the shape onto a static floor and watches its contacts. The floor's z axis points up.
Drop drop(const Geometry& g, const Geometry& floor, int frames = 240) {
    World w;
    BodyDesc fd;
    fd.type = MotionType::Static;
    fd.rotation = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });
    const BodyId f = w.createBody(fd);
    ShapeDesc slab;
    slab.geometry = floor;
    w.addShape(f, slab);

    BodyDesc bd;
    bd.position = { 0.1f, 1, 0.05f };
    bd.rotation = Quat::axisAngle(normalize(Vec3{ 1, 0.4f, 0.7f }), 1.1f);
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    w.addShape(b, sd);

    Drop d;
    for (int i = 0; i < frames; ++i) {
        w.step(kDt);
        for (const ContactInfo& c : w.contacts()) {
            if (!(std::fabs(length(c.normal) - 1.0f) < 1e-3f)) { d.unitNormals = false; }
            for (int k = 0; k < c.pointCount; ++k) { d.deepest = std::fmax(d.deepest, -c.separations[k]); }
        }
    }
    const Vec3 p = w.position(b);
    d.finite = std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
    d.height = p.y;
    d.touching = !w.contacts().empty();
    return d;
}

Drop drop(const Geometry& g, int frames = 240) {
    return drop(g, Prism{ { { -4, -4 }, { 4, -4 }, { 4, 4 }, { -4, 4 } }, -1, 0 }, frames);
}

float massOf(const Geometry& g) {
    World w;
    BodyDesc bd;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    w.addShape(b, sd);
    return w.massProperties(b).mass;
}

}

TEST(flat_and_messy_plates_hold_a_falling_box) {
    const std::vector<Geometry> plates = {
        Box{ { 2, 2, 0 } },
        Prism{ { { -2, -2 }, { 2, -2 }, { 2, 2 }, { -2, 2 } }, 0, 0 },
        Prism{ { { -2, -2 }, { -2, -2 }, { 2, -2 }, { 2, 0 }, { 2, 2 }, { -2, 2 }, { -2, -2 } }, -0.25f, 0.25f },
    };
    for (const Geometry& plate : plates) {
        const Drop d = drop(Box{ { 0.3f, 0.3f, 0.3f } }, plate);
        CHECK(d.unitNormals);
        CHECK(d.height > 0.29f && d.height < 0.56f);
    }
}

TEST(mirrored_box_weighs_the_same_as_the_box) {
    const float mass = massOf(Box{ { 0.5f, 1.0f, 1.5f } });
    CHECK_NEAR(massOf(Box{ { -0.5f, 1.0f, -1.5f } }), mass, 1e-5f);
    CHECK_NEAR(massOf(Box{ { 0.5f, -1.0f, 1.5f } }), mass, 1e-5f);
    const Drop d = drop(Box{ { -0.3f, 0.3f, -0.3f } });
    CHECK(d.unitNormals);
    CHECK_NEAR(d.height, 0.3f, 0.02f);
}

TEST(non_finite_box_sizes_become_thin_sides) {
    for (float bad : { kNan, kInf, -kInf }) {
        const float mass = massOf(Box{ { 0.5f, bad, 0.5f } });
        CHECK(mass > 0.0f && mass < 0.1f);
    }
}

TEST(degenerate_prisms_have_volume) {
    const std::vector<Prism> prisms = {
        Prism{},
        Prism{ { { 0.2f, 0.1f } }, -0.5f, 0.5f },
        Prism{ { { -0.5f, 0.0f }, { 0.5f, 0.0f } }, -0.5f, 0.5f },
        Prism{ { { -0.5f, 0.0f }, { 0.0f, 0.0f }, { 0.5f, 0.0f } }, -0.5f, 0.5f },
        Prism{ { { -0.5f, -0.5f }, { 0.5f, -0.5f }, { 0.5f, 0.5f }, { -0.5f, 0.5f } }, 0.25f, 0.25f },
        Prism{ { { -0.5f, -0.5f }, { kNan, 0.0f }, { 0.5f, -0.5f }, { 0.5f, 0.5f } }, kInf, 0.5f },
    };
    for (const Prism& p : prisms) {
        const float mass = massOf(p);
        CHECK(mass > 0.0f && mass < 1.0f);
        const Drop d = drop(p, 120);
        CHECK(d.unitNormals);
        CHECK(d.finite);
        CHECK(d.touching);
        CHECK(d.deepest < 0.05f);
    }
}

TEST(thin_rods_land_without_bouncing_up) {
    uint32_t s = 0x2545f491u;
    auto next = [&s]() {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return (float)(s & 0xffffff) / (float)0xffffff;
    };
    for (int k = 0; k < 40; ++k) {
        const Vec3 axis = normalize(Vec3{ next() * 2 - 1, next() * 2 - 1, next() * 2 - 1 });
        const Quat rotation = Quat::axisAngle(axis, next() * 6.2831853f);
        World w;
        BodyDesc fd;
        fd.type = MotionType::Static;
        fd.rotation = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });
        const BodyId f = w.createBody(fd);
        ShapeDesc slab;
        slab.geometry = Prism{ { { -4, -4 }, { 4, -4 }, { 4, 4 }, { -4, 4 } }, -1, 0 };
        w.addShape(f, slab);
        BodyDesc bd;
        bd.position = { 0.1f, 1, 0.05f };
        bd.rotation = rotation;
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Box{ { 0.005f, 0.005f, 0.5f } };
        w.addShape(b, sd);
        float highest = 0.0f, fastestSpin = 0.0f;
        for (int i = 0; i < 240; ++i) {
            w.step(kDt);
            if (i > 40) { highest = std::fmax(highest, w.position(b).y); }
            fastestSpin = std::fmax(fastestSpin, length(w.angularVelocity(b)));
        }
        CHECK(highest < 1.0f);
        CHECK(fastestSpin < 1000.0f);
    }
}

TEST(prism_geometry_reads_back_cleaned) {
    World w;
    BodyDesc bd;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Prism{ { { 0.5f, 0.5f }, { 0.5f, -0.5f }, { 0.5f, -0.5f }, { -0.5f, -0.5f }, { -0.5f, 0.5f } }, 0.5f, -0.5f };
    const ShapeId s = w.addShape(b, sd);
    const Prism p = std::get<Prism>(w.geometry(s));
    CHECK(p.polygon.size() == 4);
    CHECK(p.zMin < p.zMax);
}

TEST(hull_without_finite_points_is_a_tiny_solid) {
    const float mass = massOf(Hull{ { { kNan, 0, 0 }, { 0, kInf, 0 } } });
    CHECK(mass > 0.0f && mass < 1e-6f);
    CHECK(massOf(Hull{}) > 0.0f);
}

TEST(degenerate_spheres_and_cylinders_stay_finite) {
    const std::vector<Geometry> shapes = {
        Sphere{ 0.0f },
        Sphere{ kNan },
        Sphere{ -kInf },
        Cylinder{ 0.0f, 0.0f },
        Cylinder{ 0.0f, 0.3f },
        Cylinder{ 0.3f, 0.0f },
        Cylinder{ kNan, 0.3f },
        Cylinder{ 0.3f, kInf },
        Cylinder{ -0.3f, -0.2f },
    };
    for (const Geometry& g : shapes) {
        const float mass = massOf(g);
        CHECK(mass > 0.0f && mass < 1.0f);
        const Drop d = drop(g, 120);
        CHECK(d.unitNormals);
        CHECK(d.finite);
        CHECK(d.touching);
        CHECK(d.deepest < 0.05f);
    }
    const Drop d = drop(Sphere{ -0.3f });
    CHECK_NEAR(d.height, 0.3f, 0.01f);
}

TEST(sphere_and_cylinder_geometry_reads_back_cleaned) {
    World w;
    BodyDesc bd;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Sphere{ kNan };
    const ShapeId sphere = w.addShape(b, sd);
    sd.geometry = Cylinder{ -0.5f, 0.0f };
    const ShapeId cylinder = w.addShape(b, sd);
    const float r = std::get<Sphere>(w.geometry(sphere)).radius;
    CHECK(r > 0.0f && r < 1e-3f);
    const Cylinder c = std::get<Cylinder>(w.geometry(cylinder));
    CHECK(c.radius == 0.5f);
    CHECK(c.halfHeight > 0.0f && c.halfHeight < 0.01f);
}
