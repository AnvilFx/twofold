#include <cmath>
#include <memory>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

constexpr int kProbes = 5;
constexpr int kTargets = 7;

bool close(const Vec3& a, const Vec3& b, float tol) { return length(a - b) <= tol; }

ShapeId addShape(World& w, const Geometry& g, const Vec3& position, const Quat& rotation = {}, uint32_t category = 1,
                 MotionType type = MotionType::Static) {
    BodyDesc bd;
    bd.type = type;
    bd.position = position;
    bd.rotation = rotation;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    sd.category = category;
    return w.addShape(b, sd);
}

// A square of side 2 * half in the local yz plane.
std::shared_ptr<TriangleMesh> wall(float half) {
    const std::vector<Vec3> v{ { 0, -half, -half }, { 0, half, -half }, { 0, half, half }, { 0, -half, half } };
    return std::make_shared<TriangleMesh>(v, std::vector<uint32_t>{ 0, 1, 2, 0, 2, 3 });
}

// A flat grid in the xz plane, cells of size one from -half to half.
std::shared_ptr<TriangleMesh> grid(int half) {
    std::vector<Vec3> v;
    std::vector<uint32_t> idx;
    const int n = 2 * half + 1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { v.push_back({ (float)(i - half), 0.0f, (float)(j - half) }); }
    }
    for (int i = 0; i + 1 < n; ++i) {
        for (int j = 0; j + 1 < n; ++j) {
            const uint32_t a = i * n + j, b = a + 1, c = a + n, d = c + 1;
            idx.insert(idx.end(), { a, c, d, a, d, b });
        }
    }
    return std::make_shared<TriangleMesh>(v, idx);
}

// Query geometry reaching 0.5 from its origin along x.
Geometry probe(int kind) {
    switch (kind) {
        case 0: return Sphere{ 0.5f };
        case 1: return Box{ { 0.5f, 0.5f, 0.5f } };
        case 2: return Cylinder{ 0.5f, 0.5f };
        case 3: return Prism{ { { -0.5f, -0.5f }, { 0.5f, -0.5f }, { 0.5f, 0.5f }, { -0.5f, 0.5f } }, -0.5f, 0.5f };
        default: {
            std::vector<Vec3> cube;
            for (int i = 0; i < 8; ++i) { cube.push_back({ i & 1 ? 0.5f : -0.5f, i & 2 ? 0.5f : -0.5f, i & 4 ? 0.5f : -0.5f }); }
            return Hull{ cube, 64 };
        }
    }
}

// A shape of each kind with a surface facing +x at x = 0.5, solid behind it except for the mesh wall.
ShapeId addTarget(World& w, int kind, uint32_t category = 1) {
    if (kind < kProbes) { return addShape(w, probe(kind), {}, {}, category); }
    if (kind == kProbes) { return addShape(w, HalfSpace{}, { 0.5f, 0, 0 }, Quat::fromTo({ 0, 0, 1 }, { 1, 0, 0 }), category); }
    return addShape(w, Mesh{ wall(5.0f) }, { 0.5f, 0, 0 }, {}, category);
}

Transform at(const Vec3& p, const Quat& q = {}) { return { p, q }; }

}

TEST(overlap_point_against_each_shape_kind) {
    for (int t = 0; t < kTargets; ++t) {
        World w;
        const ShapeId id = addTarget(w, t);
        std::vector<ShapeId> found;
        if (t == kTargets - 1) {
            CHECK(w.overlapPoint({ 0.5f, 1, 1 }, found));
            CHECK(!w.overlapPoint({ 0.49f, 1, 1 }, found));
            CHECK(found.empty());
            continue;
        }
        CHECK(w.overlapPoint({ 0.49f, 0, 0 }, found));
        CHECK(found.size() == 1 && found[0] == id);
        CHECK(!w.overlapPoint({ 0.51f, 0, 0 }, found));
        CHECK(found.empty());
        CHECK(w.overlapPoint({ 0, 0.2f, 0.3f }, found));
    }

    // The corner of a sphere's bounds is outside it.
    World w;
    addShape(w, Sphere{ 0.5f }, {});
    std::vector<ShapeId> found;
    CHECK(!w.overlapPoint({ 0.45f, 0.45f, 0.45f }, found));
}

TEST(overlap_geometry_against_each_shape_kind) {
    for (int t = 0; t < kTargets; ++t) {
        for (int g = 0; g < kProbes; ++g) {
            World w;
            const ShapeId id = addTarget(w, t);
            std::vector<ShapeId> found;
            CHECK(w.overlapGeometry(probe(g), at({ 0.99f, 0, 0 }), found));
            CHECK(found.size() == 1 && found[0] == id);
            CHECK(!w.overlapGeometry(probe(g), at({ 1.01f, 0, 0 }), found));
            CHECK(found.empty());
            // The mesh has no inside.
            CHECK(w.overlapGeometry(probe(g), at({ -0.2f, 0, 0 }), found) == (t != kTargets - 1));
        }
    }
}

TEST(overlap_geometry_is_exact_where_bounds_overlap) {
    World w;
    addShape(w, Box{}, {});
    std::vector<ShapeId> found;
    CHECK(!w.overlapGeometry(Sphere{ 0.5f }, at({ 0.9f, 0.9f, 0.9f }), found));
    CHECK(w.overlapGeometry(Sphere{ 0.5f }, at({ 0.75f, 0.75f, 0.75f }), found));

    // A box turned 45 degrees about z reaches 0.707 along x.
    const Quat turn = Quat::axisAngle({ 0, 0, 1 }, 0.78539816f);
    CHECK(w.overlapGeometry(Box{}, at({ 1.197f, 0, 0 }, turn), found));
    CHECK(!w.overlapGeometry(Box{}, at({ 1.217f, 0, 0 }, turn), found));
    CHECK(!w.overlapGeometry(Box{}, at({ 1.197f, 0, 1.01f }, turn), found));

    CHECK(!w.overlapGeometry(HalfSpace{}, at({}), found));
    CHECK(!w.overlapGeometry(Mesh{ wall(1.0f) }, at({}), found));
}

TEST(overlap_aabb_uses_bounds) {
    World w;
    const ShapeId sphere = addShape(w, Sphere{ 0.5f }, {});
    const ShapeId ground = addShape(w, HalfSpace{}, { 0, -1, 0 }, kUp);
    const ShapeId terrain = addShape(w, Mesh{ grid(2) }, { 20, 0, 0 });
    std::vector<ShapeId> found;
    CHECK(w.overlapAabb({ { 0.4f, 0.4f, 0.4f }, { 1, 1, 1 } }, found));
    CHECK(found.size() == 1 && found[0] == sphere);
    CHECK(!w.overlapAabb({ { 5, 0, 5 }, { 6, 1, 6 } }, found));
    CHECK(w.overlapAabb({ { 5, -1.1f, 5 }, { 6, 1, 6 } }, found));
    CHECK(found.size() == 1 && found[0] == ground);
    CHECK(w.overlapAabb({ { 21.5f, -0.1f, 1.5f }, { 25, 0.1f, 5 } }, found));
    CHECK(found.size() == 1 && found[0] == terrain);
    CHECK(w.overlapAabb({ { -30, -30, -30 }, { 30, 30, 30 } }, found));
    CHECK(found.size() == 3 && found[0] == sphere && found[1] == ground && found[2] == terrain);
}

TEST(cast_geometry_hits_each_shape_kind) {
    for (int t = 0; t < kTargets; ++t) {
        for (int g = 0; g < kProbes; ++g) {
            World w;
            const ShapeId id = addTarget(w, t);
            CastHit hit;
            CHECK(w.castGeometry(probe(g), at({ 3, 0, 0 }), { -4, 0, 0 }, hit));
            CHECK(hit.shape == id && hit.body == w.shapeBody(id));
            CHECK(!hit.initialOverlap);
            CHECK_NEAR(hit.fraction, 0.5f, 1e-4f);
            CHECK(close(hit.normal, { 1, 0, 0 }, 1e-3f));
            CHECK_NEAR(hit.point.x, 0.5f, 1e-3f);

            // Stopping just short misses.
            CHECK(!w.castGeometry(probe(g), at({ 3, 0, 0 }), { -1.99f, 0, 0 }, hit));
            CHECK(w.castGeometry(probe(g), at({ 3, 0, 0 }), { -2.01f, 0, 0 }, hit));
            CHECK_NEAR(hit.fraction, 2.0f / 2.01f, 1e-4f);
        }
    }
}

TEST(cast_geometry_that_just_misses) {
    for (int t : { 0, 1, kTargets - 1 }) {
        World w;
        addTarget(w, t);
        CastHit hit;
        // Passing the target's +y side.
        const float side = t == kTargets - 1 ? 5.0f : 0.5f;
        CHECK(!w.castGeometry(Sphere{ 0.5f }, at({ 3, side + 0.51f, 0 }), { -6, 0, 0 }, hit));
        if (t == 1) {
            CHECK(w.castGeometry(Sphere{ 0.5f }, at({ 3, side + 0.49f, 0 }), { -6, 0, 0 }, hit));
            CHECK(hit.normal.y > 0.5f);
        }
        // Moving away from a shape just out of reach.
        CHECK(!w.castGeometry(Box{}, at({ 1.001f, 0, 0 }), { 1, 0, 0 }, hit));
        CHECK(!w.castGeometry(Box{}, at({ 1.001f, 0, 0 }), { 0, 1, 0 }, hit));
    }
}

TEST(cast_geometry_reports_initial_overlap) {
    for (int t = 0; t < kTargets; ++t) {
        World w;
        const ShapeId id = addTarget(w, t);
        CastHit hit;
        // Overlapping by 0.1 and moving away.
        CHECK(w.castGeometry(Box{}, at({ 0.9f, 0, 0 }), { 1, 0, 0 }, hit));
        CHECK(hit.shape == id);
        CHECK(hit.initialOverlap);
        CHECK(hit.fraction == 0.0f);
        CHECK(close(hit.normal, { 1, 0, 0 }, 1e-3f));
        CHECK(hit.point.x > 0.35f && hit.point.x < 0.55f);

        // No translation tests overlap.
        CHECK(w.castGeometry(Sphere{ 0.5f }, at({ 0.9f, 0, 0 }), {}, hit));
        CHECK(hit.initialOverlap);
        CHECK(!w.castGeometry(Sphere{ 0.5f }, at({ 1.1f, 0, 0 }), {}, hit));
    }

    // Either side of a mesh pushes out toward that side.
    World w;
    addTarget(w, kTargets - 1);
    CastHit hit;
    CHECK(w.castGeometry(Sphere{ 0.5f }, at({ 0.7f, 0, 0 }), { 0, 1, 0 }, hit));
    CHECK(hit.initialOverlap && close(hit.normal, { 1, 0, 0 }, 1e-3f));
    CHECK(w.castGeometry(Sphere{ 0.5f }, at({ 0.3f, 0, 0 }), { 0, 1, 0 }, hit));
    CHECK(hit.initialOverlap && close(hit.normal, { -1, 0, 0 }, 1e-3f));
}

TEST(cast_geometry_finds_the_closest_shape) {
    World w;
    addShape(w, Box{}, { 0, 0, 0 });
    addShape(w, Box{}, { 3, 0, 0 });
    const ShapeId nearest = addShape(w, Box{}, { 6, 0, 0 });
    CastHit hit;
    CHECK(w.castGeometry(Sphere{ 0.5f }, at({ 10, 0, 0 }), { -20, 0, 0 }, hit));
    CHECK(hit.shape == nearest);
    CHECK_NEAR(hit.fraction, 0.15f, 1e-4f);

    // A box turned 45 degrees about z touches with its edge.
    CHECK(w.castGeometry(Box{}, at({ 10, 0, 0 }, Quat::axisAngle({ 0, 0, 1 }, 0.78539816f)), { -20, 0, 0 }, hit));
    CHECK_NEAR(hit.fraction, (10.0f - 6.5f - 0.70710678f) / 20.0f, 1e-4f);
    CHECK(close(hit.normal, { 1, 0, 0 }, 1e-3f));

    // Coincident shapes go to the lower index.
    World c;
    const ShapeId first = addShape(c, Box{}, {});
    addShape(c, Box{}, {});
    CHECK(c.castGeometry(Sphere{ 0.5f }, at({ 0, 5, 0 }), { 0, -10, 0 }, hit));
    CHECK(hit.shape == first);
    CHECK(c.castGeometry(Sphere{ 0.5f }, at({ 0, 0.5f, 0 }), { 0, -10, 0 }, hit));
    CHECK(hit.shape == first && hit.initialOverlap);
}

TEST(cast_geometry_across_mesh_terrain) {
    World w;
    const ShapeId terrain = addShape(w, Mesh{ grid(10) }, {});
    CastHit hit;
    CHECK(w.castGeometry(Sphere{ 0.5f }, at({ -8, 3, -7.3f }), { 16, -6, 14 }, hit));
    CHECK(hit.shape == terrain);
    CHECK(!hit.initialOverlap);
    CHECK_NEAR(hit.fraction, 2.5f / 6.0f, 1e-4f);
    CHECK(close(hit.normal, { 0, 1, 0 }, 1e-3f));
    CHECK_NEAR(hit.point.y, 0.0f, 1e-3f);

    CHECK(w.castGeometry(Box{}, at({ 4, 3, 2 }, Quat::axisAngle({ 1, 1, 0 }, 0.3f)), { -9, -6, 1 }, hit));
    const float drop = 3.0f - 6.0f * hit.fraction;
    std::vector<ShapeId> found;
    CHECK(!w.overlapGeometry(Box{}, at({ 4 - 9 * hit.fraction, drop + 0.001f, 2 + hit.fraction }, Quat::axisAngle({ 1, 1, 0 }, 0.3f)), found));
    CHECK(w.overlapGeometry(Box{}, at({ 4 - 9 * hit.fraction, drop - 0.001f, 2 + hit.fraction }, Quat::axisAngle({ 1, 1, 0 }, 0.3f)), found));

    // Skimming over the top.
    CHECK(!w.castGeometry(Sphere{ 0.5f }, at({ -8, 0.52f, -7 }), { 16, 0, 14 }, hit));
    CHECK(w.castGeometry(Sphere{ 0.5f }, at({ -8, 0.48f, -7 }), { 16, 0, 14 }, hit));
    CHECK(hit.initialOverlap);
}

TEST(cast_ray_all_reports_every_hit_in_order) {
    World w;
    std::vector<ShapeId> ids;
    for (int k = 0; k < kProbes; ++k) { ids.push_back(addShape(w, probe(k), { 2.0f * k, 0, 0 })); }
    ids.push_back(addShape(w, Mesh{ wall(2.0f) }, { 10, 0, 0 }));
    ids.push_back(addShape(w, HalfSpace{}, { 12, 0, 0 }, Quat::fromTo({ 0, 0, 1 }, { -1, 0, 0 })));

    std::vector<RayHit> hits;
    CHECK(w.castRayAll({ -2, 0.1f, 0.1f }, { 20, 0, 0 }, hits));
    CHECK(hits.size() == ids.size());
    for (size_t i = 0; i < hits.size(); ++i) {
        CHECK(hits[i].shape == ids[i]);
        CHECK(hits[i].body == w.shapeBody(ids[i]));
        CHECK(hits[i].normal.x < -0.9f);
        CHECK(i == 0 || hits[i].fraction > hits[i - 1].fraction);
    }
    CHECK_NEAR(hits[1].fraction, 3.5f / 20.0f, 1e-5f);
    CHECK_NEAR(hits[6].fraction, 14.0f / 20.0f, 1e-5f);
    RayHit first;
    CHECK(w.castRay({ -2, 0.1f, 0.1f }, { 20, 0, 0 }, first));
    CHECK(first.shape == hits[0].shape && first.fraction == hits[0].fraction);

    // A ray starting inside a shape skips it, and one that falls short stops early.
    CHECK(w.castRayAll({ 0, 0.1f, 0.1f }, { 9, 0, 0 }, hits));
    CHECK(hits.size() == 4 && hits[0].shape == ids[1] && hits[3].shape == ids[4]);
    CHECK(!w.castRayAll({ 0, 0.1f, 0.1f }, { 1, 0, 0 }, hits));
    CHECK(hits.empty());

    // Coincident layers come in index order.
    World c;
    const Prism layer{ { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } }, -0.05f, 0.05f };
    std::vector<ShapeId> layers;
    for (int i = 0; i < 3; ++i) { layers.push_back(addShape(c, layer, { 0, 1, 0 })); }
    c.removeShape(layers[0]);
    layers.erase(layers.begin());
    layers.push_back(addShape(c, layer, { 0, 1, 0 }));
    CHECK(c.castRayAll({ 0.3f, 1.2f, 5 }, { 0, 0, -10 }, hits));
    CHECK(hits.size() == 3);
    for (size_t i = 1; i < hits.size(); ++i) {
        CHECK(hits[i].fraction == hits[0].fraction);
        CHECK(hits[i].shape.index > hits[i - 1].shape.index);
    }
}

TEST(queries_respect_masks) {
    World w;
    const ShapeId a = addShape(w, Sphere{ 0.5f }, { 0, 0, 0 }, {}, 1);
    const ShapeId b = addShape(w, Box{}, { 0.2f, 0, 0 }, {}, 2);
    const ShapeId c = addShape(w, Mesh{ wall(1.0f) }, { 0.1f, 0, 0 }, {}, 4);
    const ShapeId d = addShape(w, HalfSpace{}, { 0, 0.3f, 0 }, kUp, 8);
    std::vector<ShapeId> found;
    CHECK(w.overlapPoint({ 0.1f, 0, 0 }, found));
    CHECK(found.size() == 4);
    CHECK(w.overlapPoint({ 0.1f, 0, 0 }, found, 2) && found.size() == 1 && found[0] == b);
    CHECK(w.overlapPoint({ 0.1f, 0, 0 }, found, 9) && found.size() == 2 && found[0] == a && found[1] == d);
    CHECK(!w.overlapPoint({ 0.1f, 0, 0 }, found, 0));

    CHECK(w.overlapGeometry(Sphere{ 0.15f }, at({}), found) && found.size() == 4);
    CHECK(w.overlapGeometry(Sphere{ 0.15f }, at({}), found, 4) && found.size() == 1 && found[0] == c);
    CHECK(!w.overlapGeometry(Sphere{ 0.15f }, at({}), found, 16));

    CHECK(w.overlapAabb({ { -1, -1, -1 }, { 1, 1, 1 } }, found, 6) && found.size() == 2 && found[0] == b && found[1] == c);
    CHECK(!w.overlapAabb({ { -1, -1, -1 }, { 1, 1, 1 } }, found, 0));

    // Down onto the box at 0.44, the half-space at 0.46 and the sphere at 0.467. The wall stays 0.2 away.
    CastHit hit;
    const Transform from = at({ 0.4f, 5, 0.3f });
    CHECK(w.castGeometry(Sphere{ 0.1f }, from, { 0, -10, 0 }, hit));
    CHECK(hit.shape == b);
    CHECK_NEAR(hit.fraction, 0.44f, 1e-4f);
    CHECK(w.castGeometry(Sphere{ 0.1f }, from, { 0, -10, 0 }, hit, 9));
    CHECK(hit.shape == d);
    CHECK(w.castGeometry(Sphere{ 0.1f }, from, { 0, -10, 0 }, hit, 1));
    CHECK(hit.shape == a);
    CHECK(!w.castGeometry(Sphere{ 0.1f }, from, { 0, -10, 0 }, hit, 4));
    CHECK(!w.castGeometry(Sphere{ 0.1f }, from, { 0, -10, 0 }, hit, 0));

    std::vector<RayHit> hits;
    CHECK(w.castRayAll({ 0.3f, 5, 0.2f }, { 0, -10, 0 }, hits) && hits.size() == 3);
    CHECK(hits[0].shape == b && hits[1].shape == a && hits[2].shape == d);
    CHECK(w.castRayAll({ 0.3f, 5, 0.2f }, { 0, -10, 0 }, hits, 5) && hits.size() == 1 && hits[0].shape == a);
    CHECK(!w.castRayAll({ 0.3f, 5, 0.2f }, { 0, -10, 0 }, hits, 4));
}

TEST(queries_find_planar_bodies_and_outline_pieces) {
    const std::vector<Vec2> ell{ { 0, 0 }, { 3, 0 }, { 3, 1 }, { 1, 1 }, { 1, 3 }, { 0, 3 } };
    const std::vector<std::vector<Vec2>> pieces = convexDecompose({ ell });
    CHECK(pieces.size() >= 2);

    World w;
    addShape(w, HalfSpace{}, { 0, -1, 0 }, kUp);
    BodyDesc bd;
    bd.locks = LockPlanar;
    const BodyId body = w.createBody(bd);
    for (const std::vector<Vec2>& p : pieces) {
        ShapeDesc sd;
        sd.geometry = Prism{ p, -0.1f, 0.1f };
        w.addShape(body, sd);
    }

    std::vector<ShapeId> found;
    CHECK(!w.overlapPoint({ 2, 2, 0 }, found));
    CHECK(w.overlapPoint({ 2.5f, 0.5f, 0 }, found));
    for (const ShapeId& s : found) { CHECK(w.shapeBody(s) == body); }
    CHECK(w.overlapPoint({ 0.5f, 2.5f, 0.05f }, found));
    CHECK(!w.overlapPoint({ 0.5f, 2.5f, 0.15f }, found));

    // A layer at another depth passes.
    CHECK(w.overlapGeometry(Box{ { 0.2f, 0.2f, 0.2f } }, at({ 2, 0.5f, 0.25f }), found));
    CHECK(!w.overlapGeometry(Box{ { 0.2f, 0.2f, 0.2f } }, at({ 2, 0.5f, 0.35f }), found));
    CHECK(w.overlapGeometry(Box{ { 0.2f, 0.2f, 0.2f } }, at({ 2, -0.9f, 0.35f }), found));
    CHECK(found.size() == 1 && w.shapeBody(found[0]) != body);

    // Into the inner corner along the plane, and onto the face of the layer.
    CastHit hit;
    CHECK(w.castGeometry(Sphere{ 0.25f }, at({ 2, 2, 0 }), { -3, 0, 0 }, hit));
    CHECK(hit.body == body && !hit.initialOverlap);
    CHECK_NEAR(hit.fraction, 0.25f, 1e-4f);
    CHECK(close(hit.normal, { 1, 0, 0 }, 1e-3f));
    CHECK(w.castGeometry(Sphere{ 0.25f }, at({ 2, 0.5f, 5 }), { 0, 0, -10 }, hit));
    CHECK(hit.body == body);
    CHECK_NEAR(hit.fraction, 0.465f, 1e-4f);
    CHECK(close(hit.normal, { 0, 0, 1 }, 1e-3f));

    std::vector<RayHit> hits;
    CHECK(w.castRayAll({ -1, 0.5f, 0 }, { 5, 0, 0 }, hits));
    CHECK_NEAR(hits[0].fraction, 0.2f, 1e-5f);
    for (const RayHit& h : hits) { CHECK(h.body == body); }
}

TEST(queries_follow_moved_bodies) {
    World w;
    BodyDesc bd;
    bd.position = { 0, 5, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Sphere{ 0.5f };
    const ShapeId s = w.addShape(b, sd);
    std::vector<ShapeId> found;
    CHECK(w.overlapPoint({ 0, 5, 0 }, found));
    w.setTransform(b, { 10, 0, 0 }, {});
    CHECK(!w.overlapPoint({ 0, 5, 0 }, found));
    CHECK(w.overlapPoint({ 10, 0.4f, 0 }, found) && found[0] == s);
    for (int i = 0; i < 30; ++i) { w.step(1.0f / 60.0f); }
    const Vec3 p = w.position(b);
    CHECK(p.y < -1.0f);
    CHECK(w.overlapPoint(p, found) && found[0] == s);
    CastHit hit;
    CHECK(w.castGeometry(Box{}, at({ 10, 5, 0 }), { 0, -20, 0 }, hit) && hit.shape == s);

    // A copy answers the same.
    const World copy = w;
    std::vector<RayHit> a, c;
    CHECK(w.castRayAll({ 10, 5, 0 }, { 0, -20, 0 }, a));
    CHECK(copy.castRayAll({ 10, 5, 0 }, { 0, -20, 0 }, c));
    CHECK(a.size() == c.size() && a[0].fraction == c[0].fraction && a[0].shape == c[0].shape);
}

namespace {

struct Random {
    uint32_t state = 12345u;

    float next() {
        state = state * 1664525u + 1013904223u;
        return (float)(state >> 8) * (1.0f / 16777216.0f);
    }
    float range(float lo, float hi) { return lo + (hi - lo) * next(); }
    Vec3 vec(float r) { return { range(-r, r), range(-r, r), range(-r, r) }; }
    Quat turn() { return Quat::axisAngle(vec(1.0f), range(0.0f, 6.28f)); }
};

// Every shape kind at random poses, overlapping each other in places.
World mixedScene(Random& rng) {
    World w;
    for (int i = 0; i < 15; ++i) { addShape(w, probe(i % kProbes), rng.vec(4.0f), rng.turn()); }
    addShape(w, Mesh{ grid(3) }, rng.vec(2.0f), rng.turn());
    addShape(w, Mesh{ wall(2.0f) }, rng.vec(3.0f), rng.turn());
    addShape(w, HalfSpace{}, { 0, -5, 0 }, Quat::fromTo({ 0, 0, 1 }, normalize(Vec3{ 0.2f, 1, -0.1f })));
    return w;
}

bool contains(const std::vector<ShapeId>& ids, ShapeId id) {
    for (const ShapeId& s : ids) {
        if (s == id) { return true; }
    }
    return false;
}

}

TEST(cast_geometry_agrees_with_overlap_queries) {
    Random rng;
    int hits = 0, starts = 0, misses = 0;
    for (int scene = 0; scene < 4; ++scene) {
        const World w = mixedScene(rng);
        std::vector<ShapeId> found;
        for (int i = 0; i < 150; ++i) {
            const Geometry g = probe((int)(rng.next() * kProbes) % kProbes);
            const Transform from = at(rng.vec(6.0f), rng.turn());
            const Vec3 d = rng.vec(8.0f);
            const float travel = length(d);
            CastHit hit;
            if (!w.castGeometry(g, from, d, hit)) {
                ++misses;
                for (int k = 0; k <= 8; ++k) { CHECK(!w.overlapGeometry(g, at(from.p + d * (k / 8.0f), from.q), found)); }
                continue;
            }
            if (hit.initialOverlap) {
                ++starts;
                CHECK(hit.fraction == 0.0f);
                CHECK(w.overlapGeometry(g, from, found) && found[0] == hit.shape);
                continue;
            }
            ++hits;
            CHECK(!w.overlapGeometry(g, from, found));
            CHECK(hit.fraction > 0.0f && hit.fraction < 1.0f);
            CHECK(std::fabs(length(hit.normal) - 1.0f) < 1e-3f);
            CHECK(w.overlapGeometry(g, at(from.p + d * hit.fraction, from.q), found) && contains(found, hit.shape));
            // Clear of everything just before the hit, and into the hit shape just after, unless grazing.
            if (std::fabs(dot(d, hit.normal)) < 0.05f * travel) { continue; }
            const float before = std::fmax(0.0f, hit.fraction - 1e-3f / travel);
            CHECK(!w.overlapGeometry(g, at(from.p + d * before, from.q), found));
            const float after = hit.fraction + 2e-3f / travel;
            CHECK(w.overlapGeometry(g, at(from.p + d * after, from.q), found) && contains(found, hit.shape));
        }
    }
    CHECK(hits > 100 && starts > 50 && misses > 50);
}

TEST(overlap_queries_agree_with_rays_and_points) {
    Random rng;
    rng.state = 777u;
    for (int scene = 0; scene < 4; ++scene) {
        const World w = mixedScene(rng);
        std::vector<ShapeId> points, spheres;
        std::vector<RayHit> all;
        for (int i = 0; i < 300; ++i) {
            const Vec3 p = rng.vec(6.0f);
            w.overlapPoint(p, points);
            w.overlapGeometry(Sphere{ 0.0f }, at(p), spheres);
            CHECK(points.size() == spheres.size());
            for (size_t k = 0; k < points.size(); ++k) { CHECK(points[k] == spheres[k]); }
            for (size_t k = 1; k < points.size(); ++k) { CHECK(points[k].index > points[k - 1].index); }

            const Vec3 d = rng.vec(10.0f);
            RayHit closest;
            const bool any = w.castRay(p, d, closest);
            CHECK(w.castRayAll(p, d, all) == any);
            if (!any) { continue; }
            CHECK(all[0].shape == closest.shape && all[0].fraction == closest.fraction);
            for (size_t k = 1; k < all.size(); ++k) {
                CHECK(all[k].fraction >= all[k - 1].fraction);
                CHECK(all[k].fraction > all[k - 1].fraction || all[k].shape.index > all[k - 1].shape.index);
            }
            // A shape the ray enters is one the point at the hit touches.
            for (const RayHit& h : all) { CHECK(w.overlapGeometry(Sphere{ 1e-3f }, at(h.point), spheres) && contains(spheres, h.shape)); }
        }
    }
}