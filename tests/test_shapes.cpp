#include <cmath>
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
