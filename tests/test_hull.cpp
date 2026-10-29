#include <cmath>
#include <cstdint>
#include <limits>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kPi = 3.14159265358979f;

std::vector<Vec3> spherePoints(int count, float radius) {
    std::vector<Vec3> points;
    const float golden = kPi * (3.0f - std::sqrt(5.0f));
    for (int i = 0; i < count; ++i) {
        const float y = 1.0f - 2.0f * ((float)i + 0.5f) / (float)count;
        const float r = std::sqrt(1.0f - y * y);
        points.push_back(Vec3{ r * std::cos(golden * (float)i), y, r * std::sin(golden * (float)i) } * radius);
    }
    return points;
}

std::vector<Vec3> cone(int sides, float radius, float height) {
    std::vector<Vec3> points{ { 0, height, 0 } };
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        points.push_back({ radius * std::cos(a), 0, radius * std::sin(a) });
    }
    return points;
}

// Every face is planar, faces outward, and leaves all vertices behind it. Euler's formula holds.
bool wellFormed(const HullMesh& h, float tol) {
    if (h.vertices.size() < 4 || h.faces.size() < 4) { return false; }
    size_t corners = 0;
    for (const std::vector<int>& f : h.faces) {
        if (f.size() < 3) { return false; }
        corners += f.size();
        Vec3 n;
        const Vec3& o = h.vertices[f[0]];
        for (size_t i = 1; i + 1 < f.size(); ++i) { n += cross(h.vertices[f[i]] - o, h.vertices[f[i + 1]] - o); }
        n = normalize(n);
        const float offset = dot(n, h.vertices[f[0]]);
        for (int v : f) {
            if (std::fabs(dot(n, h.vertices[v]) - offset) > tol) { return false; }
        }
        for (const Vec3& v : h.vertices) {
            if (dot(n, v) - offset > tol) { return false; }
        }
    }
    const int edges = (int)corners / 2;
    return (int)h.vertices.size() - edges + (int)h.faces.size() == 2;
}

struct Random {
    uint32_t s = 0x2545f491u;

    float next() {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return (float)(s & 0xffffff) / 8388608.0f - 1.0f;
    }
    Vec3 onSphere() {
        const float a = next() * kPi, b = next() * 0.5f * kPi;
        return { std::cos(a) * std::cos(b), std::sin(a) * std::cos(b), std::sin(b) };
    }
};

// Twice the tolerance computeHull works to.
float slack(const std::vector<Vec3>& points) {
    Vec3 lo = points[0], hi = points[0];
    for (const Vec3& p : points) {
        lo = minVec(lo, p);
        hi = maxVec(hi, p);
    }
    float reach = 0.0f;
    for (int k = 0; k < 3; ++k) { reach += std::fmax(std::fabs(lo[k]), std::fabs(hi[k])); }
    return 2.0f * std::fmax(1e-4f * length(hi - lo), 3.6e-7f * reach);
}

Vec3 faceNormal(const HullMesh& h, const std::vector<int>& f) {
    Vec3 n;
    const Vec3& o = h.vertices[f[0]];
    for (size_t i = 1; i + 1 < f.size(); ++i) { n += cross(h.vertices[f[i]] - o, h.vertices[f[i + 1]] - o); }
    return normalize(n);
}

// Every point lies behind every face.
bool holds(const HullMesh& h, const std::vector<Vec3>& points, float tol) {
    for (const std::vector<int>& f : h.faces) {
        const Vec3 n = faceNormal(h, f);
        const float offset = dot(n, h.vertices[f[0]]);
        for (const Vec3& p : points) {
            if (dot(n, p) - offset > tol) { return false; }
        }
    }
    return true;
}

// No face outline bends inward at a corner.
bool convexCorners(const HullMesh& h, float tol) {
    for (const std::vector<int>& f : h.faces) {
        const Vec3 n = faceNormal(h, f);
        const size_t c = f.size();
        for (size_t i = 0; i < c; ++i) {
            const Vec3& a = h.vertices[f[i]];
            const Vec3& b = h.vertices[f[(i + 1) % c]];
            const Vec3& d = h.vertices[f[(i + 2) % c]];
            if (dot(cross(b - a, d - b), n) < -tol * length(d - a)) { return false; }
        }
    }
    return true;
}

float deepestContact(const World& w) {
    float d = 0.0f;
    for (const ContactInfo& c : w.contacts()) {
        for (int i = 0; i < c.pointCount; ++i) { d = std::fmax(d, -c.separations[i]); }
    }
    return d;
}

BodyId addGround(World& w) {
    BodyDesc bd;
    bd.type = MotionType::Static;
    bd.rotation = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });
    const BodyId g = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = HalfSpace{};
    w.addShape(g, sd);
    return g;
}

}

TEST(hull_of_a_cube_cloud_is_the_cube) {
    std::vector<Vec3> points;
    for (int i = 0; i < 8; ++i) { points.push_back({ i & 1 ? 0.5f : -0.5f, i & 2 ? 0.5f : -0.5f, i & 4 ? 0.5f : -0.5f }); }
    for (int i = 0; i < 5; ++i) {
        const float t = -0.4f + 0.2f * (float)i;
        points.push_back({ t, 0.1f * t, -0.3f });
        points.push_back({ 0.5f, t, 0.25f * t });
        points.push_back({ t, -0.5f, t });
    }
    const HullMesh h = computeHull(points);
    CHECK(h.vertices.size() == 8);
    CHECK(h.faces.size() == 6);
    for (const std::vector<int>& f : h.faces) { CHECK(f.size() == 4); }
    CHECK(wellFormed(h, 1e-5f));
}

TEST(hull_mass_matches_box) {
    World w;
    BodyDesc bd;
    ShapeDesc sd;
    sd.geometry = Box{ { 0.5f, 1.0f, 1.5f } };
    const BodyId box = w.createBody(bd);
    w.addShape(box, sd);
    Hull hull;
    for (int i = 0; i < 8; ++i) { hull.points.push_back({ i & 1 ? 0.5f : -0.5f, i & 2 ? 1.0f : -1.0f, i & 4 ? 1.5f : -1.5f }); }
    sd.geometry = hull;
    const BodyId h = w.createBody(bd);
    w.addShape(h, sd);
    const MassProperties a = w.massProperties(box), b = w.massProperties(h);
    CHECK_NEAR(b.mass, a.mass, 1e-4f);
    CHECK_NEAR(b.inertia.r[0].x, a.inertia.r[0].x, 1e-3f);
    CHECK_NEAR(b.inertia.r[2].z, a.inertia.r[2].z, 1e-3f);
    CHECK(length(b.center) < 1e-5f);
}

TEST(hull_keeps_vertex_limit) {
    const std::vector<Vec3> points = spherePoints(2000, 1.0f);
    for (int limit : { 8, 32, 64 }) {
        const HullMesh h = computeHull(points, limit);
        CHECK((int)h.vertices.size() <= limit);
        CHECK((int)h.vertices.size() >= limit - 2);
        CHECK(wellFormed(h, 1e-4f));
    }
    const HullMesh fine = computeHull(points, 64);
    float nearest = 1.0f;
    for (const std::vector<int>& f : fine.faces) {
        Vec3 n;
        for (size_t i = 0; i < f.size(); ++i) { n += cross(fine.vertices[f[i]], fine.vertices[f[(i + 1) % f.size()]]); }
        nearest = std::fmin(nearest, dot(normalize(n), fine.vertices[f[0]]));
    }
    CHECK(nearest > 0.85f);
}

TEST(degenerate_clouds_become_thin_hulls) {
    const std::vector<std::vector<Vec3>> clouds = {
        { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0.5f, 0.5f, 0 } },
        { { 0, 0, 0 }, { 1, 1, 1 }, { 2, 2, 2 } },
        { { 3, 3, 3 } },
    };
    for (const std::vector<Vec3>& cloud : clouds) {
        const HullMesh h = computeHull(cloud);
        CHECK(wellFormed(h, 1e-4f));
    }
    CHECK(computeHull({}).vertices.empty());

    World w;
    BodyDesc bd;
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Hull{ clouds[0] };
    w.addShape(b, sd);
    CHECK(w.massProperties(b).mass > 0.0f);
}

TEST(cone_hull_rests_on_its_base) {
    World w;
    addGround(w);
    BodyDesc bd;
    bd.position = { 0, 1, 0 };
    const BodyId b = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Hull{ cone(24, 0.5f, 1.0f) };
    w.addShape(b, sd);
    for (int i = 0; i < 180; ++i) { w.step(kDt); }
    CHECK(!w.isAwake(b));
    CHECK(std::fabs(w.rotation(b).rotate({ 0, 1, 0 }).y) > 0.999f);
    CHECK(w.position(b).y < 0.01f && w.position(b).y > -0.01f);
}

TEST(tumbling_hulls_do_not_sink_into_ground) {
    float deepest = 0.0f;
    for (int i = 0; i < 16; ++i) {
        World w;
        addGround(w);
        const float t = (float)i;
        BodyDesc bd;
        bd.position = { 0, 2.0f + 0.25f * t, 0 };
        bd.rotation = Quat::axisAngle({ std::sin(1.3f * t), std::cos(2.1f * t), 0.5f }, 0.4f * t);
        bd.angularVelocity = Vec3{ std::sin(0.7f * t), std::cos(1.9f * t), std::sin(2.9f * t) } * 6.0f;
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Hull{ i % 2 ? cone(16, 0.4f, 0.8f) : spherePoints(300, 0.4f), 32 };
        w.addShape(b, sd);
        for (int f = 0; f < 720; ++f) {
            w.step(kDt);
            deepest = std::fmax(deepest, deepestContact(w));
        }
        CHECK(!w.isAwake(b));
    }
    CHECK(deepest < 0.005f);
}

TEST(deeply_overlapping_hulls_keep_their_contact) {
    int missing = 0;
    for (int i = 0; i < 400; ++i) {
        World w;
        w.setGravity({});
        const float t = (float)i;
        for (int k = 0; k < 2; ++k) {
            BodyDesc bd;
            bd.position = k ? normalize(Vec3{ std::sin(t), std::cos(1.7f * t), std::sin(2.3f * t) }) * 0.7f : Vec3{};
            const Vec3 axis{ std::cos(t + (float)k), 1, std::sin((3.1f + 0.6f * (float)k) * t) };
            bd.rotation = Quat::axisAngle(normalize(axis), (0.9f + 0.4f * (float)k) * t + (float)k);
            const BodyId b = w.createBody(bd);
            ShapeDesc sd;
            sd.geometry = Hull{ spherePoints(300, 0.4f), 32 };
            w.addShape(b, sd);
        }
        w.step(kDt);
        missing += w.contacts().empty();
    }
    CHECK(missing == 0);
}

TEST(hulls_stack_on_each_other) {
    World w;
    addGround(w);
    std::vector<Vec3> block;
    for (int i = 0; i < 8; ++i) { block.push_back({ i & 1 ? 0.5f : -0.5f, i & 2 ? 0.25f : -0.25f, i & 4 ? 0.5f : -0.5f }); }
    std::vector<BodyId> stack;
    for (int i = 0; i < 6; ++i) {
        BodyDesc bd;
        bd.position = { 0.02f * (float)(i % 2), 0.25f + 0.5f * (float)i, 0 };
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        sd.geometry = Hull{ block };
        w.addShape(b, sd);
        stack.push_back(b);
    }
    for (int i = 0; i < 300; ++i) { w.step(kDt); }
    CHECK_NEAR(w.position(stack.back()).y, 2.75f, 0.03f);
    CHECK(std::fabs(w.position(stack.back()).x) < 0.05f);
}

TEST(flattened_clouds_give_convex_hulls_holding_every_point) {
    Random r;
    for (int i = 0; i < 100; ++i) {
        std::vector<Vec3> points;
        for (int k = 0; k < 20 + 3 * i; ++k) {
            Vec3 p = r.onSphere();
            p.z *= 0.05f;
            points.push_back(p);
        }
        const HullMesh h = computeHull(points, 1000);
        CHECK(wellFormed(h, slack(points)));
        CHECK(holds(h, points, slack(points)));
    }
}

TEST(round_clouds_give_convex_hulls) {
    Random r;
    for (int i = 0; i < 300; ++i) {
        std::vector<Vec3> points;
        for (int k = 0; k < 4 + i % 61; ++k) { points.push_back(r.onSphere()); }
        const HullMesh h = computeHull(points, 1000);
        CHECK(wellFormed(h, slack(points)));
        CHECK(holds(h, points, slack(points)));
    }
}

TEST(merged_faces_of_box_clouds_stay_convex) {
    Random r;
    for (int i = 0; i < 200; ++i) {
        std::vector<Vec3> points;
        for (int k = 0; k < 8 + i % 100; ++k) {
            Vec3 p{ r.next(), r.next(), r.next() };
            const int side = k % 3;
            p[side] = p[side] < 0.0f ? -1.0f : 1.0f;
            points.push_back(p);
        }
        const HullMesh h = computeHull(points, 1000);
        CHECK(convexCorners(h, slack(points)));
        CHECK(wellFormed(h, slack(points)));
    }
}

TEST(flat_clouds_stay_inside_their_hulls) {
    Random r;
    for (int i = 0; i < 400; ++i) {
        std::vector<Vec3> points;
        for (int k = 0; k < 3 + i % 100; ++k) { points.push_back({ r.next(), r.next(), 0.0f }); }
        const HullMesh h = computeHull(points, 1000);
        CHECK(wellFormed(h, slack(points)));
        CHECK(holds(h, points, slack(points)));
    }
}

TEST(lone_points_far_out_give_small_hulls) {
    for (const Vec3& p : { Vec3{ 1e4f, 0, 0 }, Vec3{ -3e5f, 2e5f, 7e4f } }) {
        HullMesh h = computeHull({ p, p + Vec3{ 0, 1e-3f, 0 } });
        CHECK(!h.vertices.empty());
        for (Vec3& v : h.vertices) {
            v -= p;
            CHECK(length(v) < 1e-3f * length(p));
        }
        CHECK(wellFormed(h, slack(h.vertices)));
    }
}

TEST(non_finite_points_are_ignored) {
    std::vector<Vec3> points;
    for (int i = 0; i < 8; ++i) { points.push_back({ i & 1 ? 0.5f : -0.5f, i & 2 ? 0.5f : -0.5f, i & 4 ? 0.5f : -0.5f }); }
    const float nan = std::numeric_limits<float>::quiet_NaN(), inf = std::numeric_limits<float>::infinity();
    points.insert(points.begin() + 3, { nan, 0, 0 });
    points.push_back({ 0, -inf, 0 });
    const HullMesh h = computeHull(points);
    CHECK(h.vertices.size() == 8);
    CHECK(h.faces.size() == 6);
    CHECK(wellFormed(h, 1e-5f));
    CHECK(computeHull({ { nan, nan, nan } }).vertices.empty());
}

TEST(hull_far_from_its_origin_collides_like_one_at_the_origin) {
    std::vector<Vec3> slab;
    for (int i = 0; i < 8; ++i) { slab.push_back({ i & 1 ? 2.0f : -2.0f, i & 2 ? 0.0f : -0.2f, i & 4 ? 2.0f : -2.0f }); }
    const std::vector<Vec3> round = spherePoints(300, 0.5f);
    float height[2] = {};
    size_t corners[2] = {};
    for (int far = 0; far < 2; ++far) {
        const Vec3 offset = far ? Vec3{ 1e4f, -6e3f, 3e3f } : Vec3{};
        World w;
        BodyDesc gd;
        gd.type = MotionType::Static;
        gd.position = -offset;
        const BodyId g = w.createBody(gd);
        ShapeDesc gs;
        std::vector<Vec3> points;
        for (const Vec3& p : slab) { points.push_back(p + offset); }
        gs.geometry = Hull{ points };
        w.addShape(g, gs);
        BodyDesc bd;
        bd.position = Vec3{ 0, 0.6f, 0 } - offset;
        const BodyId b = w.createBody(bd);
        ShapeDesc sd;
        points.clear();
        for (const Vec3& p : round) { points.push_back(p + offset); }
        sd.geometry = Hull{ points };
        w.addShape(b, sd);
        corners[far] = computeHull(points).vertices.size();
        for (int i = 0; i < 120; ++i) {
            w.step(kDt);
            for (const ContactInfo& c : w.contacts()) { CHECK(std::fabs(c.normal.y) > 0.9999f); }
        }
        CHECK(deepestContact(w) < 0.005f);
        height[far] = w.centerOfMass(b).y;
    }
    CHECK(corners[0] == corners[1]);
    const HullMesh ball = computeHull(round);
    float inner = 0.5f;
    for (const std::vector<int>& f : ball.faces) { inner = std::fmin(inner, dot(faceNormal(ball, f), ball.vertices[f[0]])); }
    for (float y : height) { CHECK(y > inner - 0.005f && y < 0.505f); }
}
