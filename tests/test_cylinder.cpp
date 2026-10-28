#include <cmath>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;

struct Piece {
    bool cylinder;
    Vec3 half;
    float radius, halfHeight;
    Vec3 position;
    Quat rotation;
};

Piece box(const Vec3& half, const Vec3& p, const Quat& q) { return { false, half, 0.0f, 0.0f, p, normalize(q) }; }

Piece cylinder(float r, float h, const Vec3& p, const Quat& q) { return { true, {}, r, h, p, normalize(q) }; }

float extent(const Piece& s, const Vec3& d) {
    const Vec3 l = s.rotation.unrotate(d);
    if (s.cylinder) {
        const float side = std::sqrt(std::fmax(1.0f - l.z * l.z, 0.0f));
        return dot(s.position, d) + s.halfHeight * std::fabs(l.z) + s.radius * side;
    }
    return dot(s.position, d) + s.half.x * std::fabs(l.x) + s.half.y * std::fabs(l.y) + s.half.z * std::fabs(l.z);
}

struct Overlap {
    float deepest = 0.0f;
    float alongNormal = 0.0f;
};

// Steps the pair once without gravity. alongNormal is how far the pieces overlap along the contact normal.
Overlap collidePieces(const Piece& pa, const Piece& pb) {
    World w;
    w.setGravity({});
    BodyId ids[2];
    const Piece* pieces[2] = { &pa, &pb };
    for (int k = 0; k < 2; ++k) {
        BodyDesc bd;
        bd.position = pieces[k]->position;
        bd.rotation = pieces[k]->rotation;
        ids[k] = w.createBody(bd);
        ShapeDesc sd;
        if (pieces[k]->cylinder) {
            sd.geometry = Cylinder{ pieces[k]->radius, pieces[k]->halfHeight };
        } else {
            sd.geometry = Box{ pieces[k]->half };
        }
        w.addShape(ids[k], sd);
    }
    w.step(kDt);
    Overlap out;
    for (const ContactInfo& c : w.contacts()) {
        const Piece& a = c.bodyA == ids[0] ? pa : pb;
        const Piece& b = c.bodyA == ids[0] ? pb : pa;
        out.alongNormal = std::fmax(out.alongNormal, extent(a, c.normal) + extent(b, -c.normal));
        for (int i = 0; i < c.pointCount; ++i) { out.deepest = std::fmax(out.deepest, -c.separations[i]); }
    }
    return out;
}

}

TEST(tilted_line_contacts_reach_the_overlap) {
    const Vec3 half{ 0.5f, 0.3f, 0.8f };
    const Overlap cases[] = {
        collidePieces(box(half, {}, Quat{ 0.343289f, -0.148065f, 0.350591f, -0.858671f }),
                      cylinder(0.4f, 0.6f, { 1.320397f, -0.171163f, 0.628485f },
                               Quat{ 0.272782f, -0.802671f, 0.164336f, -0.504285f })),
        collidePieces(box(half, {}, Quat{ 0.588604f, -0.533177f, -0.147332f, -0.589543f }),
                      cylinder(0.4f, 0.6f, { 0.427221f, -0.140790f, 0.999956f },
                               Quat{ 0.121856f, 0.782841f, -0.366244f, -0.488034f })),
    };
    for (const Overlap& o : cases) {
        CHECK(o.alongNormal > 0.05f);
        CHECK(o.deepest > o.alongNormal - 0.005f);
    }
}

// The expected depths are the least overlap over all directions, found by search.
TEST(deep_cylinder_overlaps_report_their_depth) {
    const Overlap side = collidePieces(
        cylinder(0.4f, 0.6f, {}, Quat{ -0.673142f, 0.625965f, 0.370122f, 0.134378f }),
        cylinder(0.7f, 0.1f, { 0.174414f, -0.823746f, -0.172703f }, Quat{ 0.488223f, -0.870501f, -0.042355f, 0.045535f }));
    CHECK_NEAR(side.deepest, 0.2707f, 0.005f);
    const Overlap rim = collidePieces(
        box({ 0.5f, 0.3f, 0.8f }, {}, Quat{ 0.503138f, 0.425952f, 0.467092f, 0.589272f }),
        cylinder(0.4f, 0.6f, { -0.570532f, -0.499719f, 0.375650f }, Quat{ -0.416410f, -0.115449f, 0.322337f, -0.842243f }));
    CHECK_NEAR(rim.deepest, 0.4963f, 0.005f);
}
