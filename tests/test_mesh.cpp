#include <cmath>
#include <memory>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kPi = 3.14159265358979f;

void run(World& w, int steps) {
    for (int i = 0; i < steps; ++i) { w.step(kDt); }
}

// A grid of cells in the xz plane, each split into two triangles facing up.
std::shared_ptr<const TriangleMesh> grid(int cells, float size, float (*height)(float, float) = nullptr) {
    std::vector<Vec3> v;
    std::vector<uint32_t> idx;
    const float half = 0.5f * size * (float)cells;
    for (int j = 0; j <= cells; ++j) {
        for (int i = 0; i <= cells; ++i) {
            const float x = -half + size * (float)i, z = -half + size * (float)j;
            v.push_back({ x, height ? height(x, z) : 0.0f, z });
        }
    }
    const uint32_t row = (uint32_t)cells + 1;
    for (uint32_t j = 0; j < (uint32_t)cells; ++j) {
        for (uint32_t i = 0; i < (uint32_t)cells; ++i) {
            const uint32_t a = j * row + i, b = a + 1, c = a + row, d = c + 1;
            idx.insert(idx.end(), { a, c, b, b, c, d });
        }
    }
    return std::make_shared<TriangleMesh>(v, idx);
}

BodyId addMesh(World& w, const std::shared_ptr<const TriangleMesh>& mesh, MotionType type = MotionType::Static,
               const Vec3& pos = {}, const Quat& rot = {}, float friction = 0.5f) {
    BodyDesc bd;
    bd.type = type;
    bd.position = pos;
    bd.rotation = rot;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Mesh{ mesh };
    sd.friction = friction;
    w.addShape(b, sd);
    return b;
}

BodyId addBody(World& w, const Vec3& pos, const Geometry& g, float friction = 0.5f, uint8_t locks = 0) {
    BodyDesc bd;
    bd.position = pos;
    bd.locks = locks;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = g;
    sd.friction = friction;
    w.addShape(b, sd);
    return b;
}

// The inside of a hemisphere of radius r, centered at the origin and open at the top.
std::shared_ptr<const TriangleMesh> bowl(float r, int rings, int segments) {
    std::vector<Vec3> v;
    std::vector<uint32_t> idx;
    v.push_back({ 0, -r, 0 });
    for (int i = 1; i <= rings; ++i) {
        const float a = 0.5f * kPi * (float)i / (float)rings;
        for (int k = 0; k < segments; ++k) {
            const float b = 2.0f * kPi * (float)k / (float)segments;
            v.push_back({ r * std::sin(a) * std::cos(b), -r * std::cos(a), r * std::sin(a) * std::sin(b) });
        }
    }
    auto at = [&](int ring, int k) { return (uint32_t)(1 + (ring - 1) * segments + (k % segments)); };
    for (int k = 0; k < segments; ++k) { idx.insert(idx.end(), { 0, at(1, k), at(1, k + 1) }); }
    for (int i = 1; i < rings; ++i) {
        for (int k = 0; k < segments; ++k) {
            idx.insert(idx.end(), { at(i, k), at(i + 1, k), at(i + 1, k + 1), at(i, k), at(i + 1, k + 1), at(i, k + 1) });
        }
    }
    return std::make_shared<TriangleMesh>(v, idx);
}

}

TEST(mesh_merges_vertices_and_drops_degenerate_triangles) {
    const std::vector<Vec3> v = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 }, { 2, 0, 0 } };
    const std::vector<uint32_t> idx = { 0, 2, 1, 3, 5, 4, 0, 1, 6, 0, 0, 2, 0, 1, 9 };
    const TriangleMesh mesh(v, idx);
    CHECK(mesh.triangleCount() == 2);
    CHECK(mesh.vertices().size() == 4);
    CHECK(mesh.indices().size() == 6);
}

TEST(box_rests_on_a_mesh_floor) {
    World w;
    addMesh(w, grid(8, 1.0f));
    const BodyId box = addBody(w, { 0.3f, 2, 0.2f }, Box{ { 0.5f, 0.5f, 0.5f } });
    run(w, 240);
    CHECK_NEAR(w.position(box).y, 0.5f, 0.01f);
    CHECK(!w.isAwake(box));
}

TEST(box_slides_across_triangle_seams_without_catching) {
    World w;
    w.setGravity({ 0, -9.81f, 0 });
    addMesh(w, grid(40, 0.25f), MotionType::Static, {}, {}, 0.0f);
    const BodyId box = addBody(w, { -3, 0.3f, 0.05f }, Box{ { 0.3f, 0.3f, 0.3f } }, 0.0f);
    run(w, 30);
    w.setLinearVelocity(box, { 4, 0, 0 });
    float lowest = 1.0f, highest = 0.0f, slowest = 4.0f, spin = 0.0f;
    for (int i = 0; i < 80; ++i) {
        w.step(kDt);
        lowest = std::fmin(lowest, w.position(box).y);
        highest = std::fmax(highest, w.position(box).y);
        slowest = std::fmin(slowest, w.linearVelocity(box).x);
        spin = std::fmax(spin, length(w.angularVelocity(box)));
    }
    CHECK(slowest > 3.99f);
    CHECK(lowest > 0.29f && highest < 0.305f);
    CHECK(spin < 0.01f);
}

TEST(ball_rolls_straight_across_triangle_seams) {
    World w;
    addMesh(w, grid(40, 0.25f));
    const BodyId ball = addBody(w, { -3, 0.3f, -1 }, Sphere{ 0.3f });
    w.setLinearVelocity(ball, { 3, 0, 1 });
    w.setAngularVelocity(ball, { 1 / 0.3f, 0, -3 / 0.3f });
    float highest = 0.0f;
    for (int i = 0; i < 90; ++i) {
        w.step(kDt);
        highest = std::fmax(highest, w.position(ball).y);
    }
    const Vec3 p = w.position(ball);
    CHECK(highest < 0.302f);
    CHECK_NEAR((p.z + 1.0f) / (p.x + 3.0f), 1.0f / 3.0f, 0.01f);
}

TEST(both_sides_of_a_triangle_collide) {
    World w;
    w.setGravity({});
    addMesh(w, std::make_shared<TriangleMesh>(std::vector<Vec3>{ { -2, 0, -2 }, { -2, 0, 2 }, { 2, 0, 0 } },
                                              std::vector<uint32_t>{ 0, 1, 2 }));
    const BodyId below = addBody(w, { 0, -1, 0 }, Sphere{ 0.25f });
    const BodyId above = addBody(w, { 0.2f, 1, 0 }, Box{ { 0.25f, 0.25f, 0.25f } });
    w.setLinearVelocity(below, { 0, 3, 0 });
    w.setLinearVelocity(above, { 0, -3, 0 });
    run(w, 60);
    CHECK(w.position(below).y < -0.24f);
    CHECK(w.position(above).y > 0.24f);
}

TEST(bodies_settle_in_a_mesh_bowl) {
    World w;
    addMesh(w, bowl(3.0f, 10, 32));
    std::vector<BodyId> ids;
    for (int i = 0; i < 16; ++i) {
        const Vec3 p{ 0.5f * (float)(i % 4) - 0.75f, 0.5f + 0.7f * (float)(i / 4), 0.4f * (float)(i % 3) - 0.4f };
        Geometry g = Sphere{ 0.3f };
        if (i % 3 == 1) { g = Box{ { 0.3f, 0.2f, 0.25f } }; }
        if (i % 3 == 2) { g = Cylinder{ 0.25f, 0.3f }; }
        ids.push_back(addBody(w, p, g));
    }
    run(w, 900);
    for (BodyId b : ids) {
        CHECK(length(w.position(b)) < 2.95f);
        CHECK(!w.isAwake(b));
    }
}

TEST(fast_ball_does_not_pass_through_a_mesh) {
    World w;
    addMesh(w, grid(2, 2.0f));
    const BodyId ball = addBody(w, { 0.3f, 3, 0.2f }, Sphere{ 0.1f });
    w.setLinearVelocity(ball, { 0, -300, 0 });
    run(w, 30);
    CHECK(w.position(ball).y > 0.0f);
}

TEST(ray_hits_either_side_of_a_mesh) {
    World w;
    const BodyId floor = addMesh(w, grid(4, 1.0f), MotionType::Static, { 0, 1, 0 });
    RayHit hit;
    CHECK(w.castRay({ 0.3f, 5, 0.1f }, { 0, -10, 0 }, hit));
    CHECK(hit.body == floor);
    CHECK_NEAR(hit.fraction, 0.4f, 1e-5f);
    CHECK_NEAR(hit.normal.y, 1.0f, 1e-5f);
    CHECK(w.castRay({ 0.3f, -3, 0.1f }, { 0, 10, 0 }, hit));
    CHECK_NEAR(hit.fraction, 0.4f, 1e-5f);
    CHECK_NEAR(hit.normal.y, -1.0f, 1e-5f);
    CHECK(!w.castRay({ 3, 5, 0 }, { 0, -10, 0 }, hit));
}

TEST(kinematic_mesh_carries_a_box) {
    World w;
    const BodyId platform = addMesh(w, grid(4, 1.0f), MotionType::Kinematic);
    const BodyId box = addBody(w, { 0, 0.5f, 0 }, Box{ { 0.5f, 0.5f, 0.5f } });
    run(w, 30);
    w.setLinearVelocity(platform, { 1, 0, 0 });
    run(w, 60);
    CHECK_NEAR(w.linearVelocity(box).x, 1.0f, 0.01f);
    CHECK_NEAR(w.position(platform).x - w.position(box).x, 0.1f, 0.01f);
    CHECK_NEAR(w.position(box).y, 0.5f, 0.01f);
}

TEST(box_rests_in_a_mesh_crease_on_two_contacts) {
    World w;
    const float s = std::sin(0.4f), c = std::cos(0.4f);
    const std::vector<Vec3> v = { { -3 * c, 3 * s, -2 }, { -3 * c, 3 * s, 2 }, { 0, 0, -2 },
                                  { 0, 0, 2 },           { 3 * c, 3 * s, -2 }, { 3 * c, 3 * s, 2 } };
    const BodyId valley = addMesh(w, std::make_shared<TriangleMesh>(v, std::vector<uint32_t>{ 0, 1, 2, 2, 1, 3, 2, 3, 4, 4, 3, 5 }));
    const BodyId box = addBody(w, { 0, 2, 0 }, Box{ { 0.3f, 0.3f, 0.3f } }, 0.6f);
    run(w, 300);
    int normals = 0;
    for (const ContactInfo& ci : w.contacts()) {
        if (ci.bodyA == valley || ci.bodyB == valley) { ++normals; }
    }
    CHECK(normals == 2);
    CHECK(length(w.linearVelocity(box)) < 0.01f);
    CHECK(w.position(box).y > 0.3f);
}

TEST(contacts_on_hilly_mesh_stay_shallow) {
    World w;
    addMesh(w, grid(16, 0.5f, [](float x, float z) { return 0.6f * std::sin(0.9f * x) * std::cos(0.7f * z); }));
    for (int i = 0; i < 16; ++i) {
        const Vec3 p{ -3.0f + 2.0f * (float)(i % 4), 1.5f + 0.5f * (float)(i / 4), -3.0f + 2.0f * (float)(i / 4) };
        Geometry g = Sphere{ 0.3f };
        if (i % 3 == 1) { g = Box{ { 0.3f, 0.2f, 0.25f } }; }
        if (i % 3 == 2) { g = Hull{ { { 0.3f, 0, 0 }, { -0.3f, 0, 0 }, { 0, 0.3f, 0 }, { 0, -0.3f, 0 }, { 0, 0, 0.3f }, { 0, 0, -0.3f } } }; }
        addBody(w, p, g);
    }
    float deepest = 0.0f;
    for (int step = 0; step < 300; ++step) {
        w.step(kDt);
        for (const ContactInfo& ci : w.contacts()) {
            for (int k = 0; k < ci.pointCount; ++k) { deepest = std::fmax(deepest, -ci.separations[k]); }
        }
    }
    CHECK(deepest < 0.05f);
}

TEST(planar_body_stays_in_its_plane_on_a_mesh) {
    World w;
    addMesh(w, grid(10, 0.5f, [](float x, float z) { return 0.2f * x * x + 0.1f * z; }));
    const BodyId box = addBody(w, { 1, 2, 0 }, Box{ { 0.4f, 0.2f, 0.3f } }, 0.5f, LockPlanar);
    run(w, 300);
    CHECK(w.position(box).z == 0.0f);
    CHECK(w.position(box).y > 0.15f && w.position(box).y < 0.3f);
    CHECK(!w.isAwake(box));
}

TEST(mesh_world_copy_steps_identically) {
    World w;
    addMesh(w, grid(12, 0.5f, [](float x, float z) { return 0.2f * std::sin(x) * std::cos(z); }));
    std::vector<BodyId> ids;
    for (int i = 0; i < 10; ++i) {
        const Vec3 p{ 0.4f * (float)(i % 4) - 0.6f, 1.0f + 0.8f * (float)i, 0.3f * (float)(i % 3) - 0.3f };
        ids.push_back(addBody(w, p, i % 2 ? Geometry{ Sphere{ 0.3f } } : Geometry{ Box{ { 0.3f, 0.2f, 0.4f } } }));
    }
    run(w, 60);
    World copy = w;
    run(w, 90);
    run(copy, 90);
    for (BodyId b : ids) {
        const Transform a = w.transform(b), c = copy.transform(b);
        CHECK(a.p.x == c.p.x && a.p.y == c.p.y && a.p.z == c.p.z);
        CHECK(a.q.x == c.q.x && a.q.y == c.q.y && a.q.z == c.q.z && a.q.w == c.q.w);
    }
}
