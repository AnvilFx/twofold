#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>

#include "test.h"
#include "twofold/twofold.h"

using namespace twofold;

namespace {

constexpr float kDt = 1.0f / 60.0f;
constexpr float kPi = 3.14159265358979f;

struct Canvas {
    int width = 0;
    int height = 0;
    std::vector<uint8_t> alpha;

    AlphaImage image() const { return { alpha.data(), width, height, width, 1 }; }
};

// Antialiased fill of a signed distance field, negative inside.
template <typename F>
Canvas paint(int width, int height, F distance) {
    Canvas c{ width, height, std::vector<uint8_t>((size_t)width * height) };
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float d = distance(Vec2{ (float)x + 0.5f, (float)y + 0.5f });
            const float cover = std::fmin(std::fmax(0.5f - d, 0.0f), 1.0f);
            c.alpha[(size_t)y * width + x] = (uint8_t)std::lround(cover * 255.0f);
        }
    }
    return c;
}

float area(const std::vector<Vec2>& poly) {
    float a = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { a += cross(poly[i], poly[(i + 1) % poly.size()]); }
    return 0.5f * a;
}

float totalArea(const std::vector<std::vector<Vec2>>& pieces) {
    float a = 0.0f;
    for (const std::vector<Vec2>& p : pieces) { a += area(p); }
    return a;
}

bool convex(const std::vector<Vec2>& poly) {
    if (poly.size() < 3 || area(poly) <= 0.0f) { return false; }
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vec2 a = poly[i], b = poly[(i + 1) % poly.size()], c = poly[(i + 2) % poly.size()];
        if (cross(b - a, c - b) < -1e-4f * length(b - a) * length(c - b)) { return false; }
    }
    return true;
}

int coverCount(const std::vector<std::vector<Vec2>>& pieces, Vec2 p) {
    int n = 0;
    for (const std::vector<Vec2>& poly : pieces) {
        bool inside = true;
        for (size_t i = 0; i < poly.size() && inside; ++i) {
            inside = cross(poly[(i + 1) % poly.size()] - poly[i], p - poly[i]) > 0.0f;
        }
        n += inside ? 1 : 0;
    }
    return n;
}

float distanceToLoop(const std::vector<Vec2>& loop, Vec2 p) {
    float best = INFINITY;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Vec2 a = loop[i], d = loop[(i + 1) % loop.size()] - a;
        const float t = std::fmin(std::fmax(dot(p - a, d) / dot(d, d), 0.0f), 1.0f);
        best = std::fmin(best, length(p - a - d * t));
    }
    return best;
}

bool crossesItself(const std::vector<Vec2>& loop) {
    const size_t n = loop.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) { continue; }
            const Vec2 p = loop[i], q = loop[(i + 1) % n], r = loop[j], s = loop[(j + 1) % n];
            const float d1 = cross(q - p, r - p), d2 = cross(q - p, s - p);
            const float d3 = cross(s - r, p - r), d4 = cross(s - r, q - r);
            if (d1 * d2 < 0.0f && d3 * d4 < 0.0f) { return true; }
        }
    }
    return false;
}

float ring(Vec2 p, Vec2 center, float outer, float inner) {
    const float r = length(p - center);
    return std::fmax(r - outer, inner - r);
}

// Pixel pieces in a body frame centered on the image, y up.
BodyId addPieces(World& w, const BodyDesc& bd, const std::vector<std::vector<Vec2>>& pieces, Vec2 center, float scale) {
    const BodyId b = w.createBody(bd);
    for (const std::vector<Vec2>& piece : pieces) {
        Prism pr;
        for (const Vec2& p : piece) { pr.polygon.push_back({ (p.x - center.x) * scale, (center.y - p.y) * scale }); }
        pr.zMin = -0.05f;
        pr.zMax = 0.05f;
        ShapeDesc sd;
        sd.geometry = pr;
        w.addShape(b, sd);
    }
    return b;
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

float deepestContact(const World& w) {
    float d = 0.0f;
    for (const ContactInfo& c : w.contacts()) {
        for (int i = 0; i < c.pointCount; ++i) { d = std::fmax(d, -c.separations[i]); }
    }
    return d;
}

}

TEST(hard_square_traces_along_pixel_edges) {
    const Canvas c = paint(20, 20, [](Vec2 p) { return (p.x > 5 && p.x < 15 && p.y > 5 && p.y < 15) ? -1.0f : 1.0f; });
    const std::vector<std::vector<Vec2>> contours = traceAlpha(c.image());
    CHECK(contours.size() == 1);
    CHECK_NEAR(std::fabs(area(contours[0])), 99.5f, 1e-3f);
    for (const Vec2& p : contours[0]) { CHECK(p.x >= 4.999f && p.x <= 15.001f && p.y >= 4.999f && p.y <= 15.001f); }
    const std::vector<Vec2> square = simplifyContour(contours[0], 0.5f);
    CHECK(square.size() == 4);
    CHECK_NEAR(std::fabs(area(square)), 100.0f, 1e-3f);
}

TEST(antialiased_disc_traces_its_radius) {
    const Canvas c = paint(80, 80, [](Vec2 p) { return length(p - Vec2{ 40, 40 }) - 30.0f; });
    const std::vector<std::vector<Vec2>> contours = traceAlpha(c.image());
    CHECK(contours.size() == 1);
    for (const Vec2& p : contours[0]) { CHECK_NEAR(length(p - Vec2{ 40, 40 }), 30.0f, 0.1f); }
    CHECK_NEAR(std::fabs(area(contours[0])), kPi * 900.0f, kPi * 900.0f * 0.002f);
}

TEST(ring_traces_outline_and_hole_in_opposite_windings) {
    const Canvas c = paint(80, 80, [](Vec2 p) { return ring(p, { 40, 40 }, 30, 15); });
    const std::vector<std::vector<Vec2>> contours = traceAlpha(c.image());
    CHECK(contours.size() == 2);
    CHECK(area(contours[0]) * area(contours[1]) < 0.0f);
}

TEST(interleaved_channel_is_read_in_place) {
    const Canvas c = paint(40, 30, [](Vec2 p) { return length(p - Vec2{ 20, 15 }) - 10.0f; });
    std::vector<uint8_t> rgba((size_t)40 * 30 * 4 + 30 * 8, 7);
    for (int y = 0; y < 30; ++y) {
        for (int x = 0; x < 40; ++x) { rgba[(size_t)y * 168 + x * 4 + 3] = c.alpha[(size_t)y * 40 + x]; }
    }
    const std::vector<std::vector<Vec2>> plain = traceAlpha(c.image());
    const std::vector<std::vector<Vec2>> strided = traceAlpha({ rgba.data() + 3, 40, 30, 168, 4 });
    CHECK(plain.size() == 1 && strided.size() == 1 && plain[0].size() == strided[0].size());
    for (size_t i = 0; i < plain[0].size(); ++i) {
        CHECK(plain[0][i].x == strided[0][i].x && plain[0][i].y == strided[0][i].y);
    }
}

TEST(simplified_contour_stays_within_tolerance) {
    std::vector<Vec2> circle;
    for (int i = 0; i < 2000; ++i) {
        const float a = 2.0f * kPi * (float)i / 2000.0f;
        circle.push_back({ 100.0f * std::cos(a), 100.0f * std::sin(a) });
    }
    const std::vector<Vec2> simple = simplifyContour(circle, 0.25f);
    CHECK(simple.size() < 120);
    for (const Vec2& p : circle) { CHECK(distanceToLoop(simple, p) <= 0.2501f); }
}

TEST(concave_polygon_splits_into_convex_pieces) {
    const std::vector<Vec2> l{ { 0, 0 }, { 2, 0 }, { 2, 1 }, { 1, 1 }, { 1, 2 }, { 0, 2 } };
    const std::vector<std::vector<Vec2>> pieces = convexDecompose({ l });
    CHECK(pieces.size() == 2);
    for (const std::vector<Vec2>& p : pieces) { CHECK(convex(p)); }
    CHECK_NEAR(totalArea(pieces), 3.0f, 1e-5f);
}

TEST(holes_and_islands_fill_even_odd) {
    const std::vector<Vec2> outer{ { 0, 0 }, { 10, 0 }, { 10, 10 }, { 0, 10 } };
    const std::vector<Vec2> hole{ { 3, 3 }, { 7, 3 }, { 7, 7 }, { 3, 7 } };
    const std::vector<Vec2> island{ { 4, 4 }, { 6, 4 }, { 6, 6 }, { 4, 6 } };
    const std::vector<std::vector<Vec2>> pieces = convexDecompose({ island, outer, hole }, 6);
    CHECK_NEAR(totalArea(pieces), 88.0f, 1e-4f);
    for (const std::vector<Vec2>& p : pieces) { CHECK(convex(p) && p.size() <= 6); }
    CHECK(coverCount(pieces, { 1.3f, 0.7f }) == 1);
    CHECK(coverCount(pieces, { 3.5f, 5.2f }) == 0);
    CHECK(coverCount(pieces, { 5.1f, 5.2f }) == 1);
    CHECK(coverCount(pieces, { 8.6f, 5.3f }) == 1);
}

TEST(hole_bridges_to_a_vertex_on_its_own_row) {
    const std::vector<Vec2> outer{ { 0, 0 },  { 14, 0 }, { 10, 5 },  { 9, 6 },   { 16, 2 },  { 30, 0 },
                                   { 30, 4 }, { 20, 5 }, { 30, 6 }, { 30, 10 }, { 0, 10 } };
    const std::vector<Vec2> hole{ { 4, 4 }, { 8, 5 }, { 4, 6 } };
    const std::vector<std::vector<Vec2>> pieces = convexDecompose({ outer, hole });
    CHECK_NEAR(totalArea(pieces), area(outer) - std::fabs(area(hole)), 1e-4f);
    CHECK(coverCount(pieces, { 12.2f, 3.1f }) == 0);
    CHECK(coverCount(pieces, { 6.1f, 5.05f }) == 0);
    CHECK(coverCount(pieces, { 2.1f, 8.3f }) == 1);
}

TEST(random_star_polygons_decompose_exactly) {
    uint32_t s = 0x2545f491u;
    auto next = [&]() {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return (float)(s & 0xffffff) / 16777216.0f;
    };
    for (int trial = 0; trial < 200; ++trial) {
        const int n = 5 + trial % 60;
        std::vector<Vec2> poly;
        for (int i = 0; i < n; ++i) {
            const float a = 2.0f * kPi * ((float)i + 0.8f * next()) / (float)n;
            const float r = 1.0f + 9.0f * next();
            poly.push_back({ r * std::cos(a), r * std::sin(a) });
        }
        if (trial % 2) { std::reverse(poly.begin(), poly.end()); }
        const int limit = 3 + trial % 10;
        const std::vector<std::vector<Vec2>> pieces = convexDecompose({ poly }, limit);
        CHECK_NEAR(totalArea(pieces), std::fabs(area(poly)), std::fabs(area(poly)) * 1e-4f);
        for (const std::vector<Vec2>& p : pieces) { CHECK(convex(p) && (int)p.size() <= limit); }
    }
}

TEST(alpha_pieces_cover_the_shape_to_the_pixel) {
    auto shape = [](Vec2 p) {
        const float bar = std::fmax(std::fabs(p.x - 60.0f) - 6.0f, std::fabs(p.y - 50.0f) - 44.0f);
        return std::fmin(ring(p, { 60, 50 }, 36, 22), bar);
    };
    const Canvas c = paint(120, 100, shape);
    const std::vector<std::vector<Vec2>> pieces = alphaPieces(c.image(), 128, 0.25f, 8);
    for (const std::vector<Vec2>& p : pieces) { CHECK(convex(p) && p.size() <= 8); }
    int checked = 0;
    for (float y = 0.25f; y < 100.0f; y += 0.5f) {
        for (float x = 0.25f; x < 120.0f; x += 0.5f) {
            const float d = shape({ x, y });
            if (std::fabs(d) < 0.5f) { continue; }
            CHECK(coverCount(pieces, { x, y }) == (d < 0.0f ? 1 : 0));
            ++checked;
        }
    }
    CHECK(checked > 40000);
}

TEST(box_settles_inside_a_traced_cup) {
    const Canvas c = paint(100, 100, [](Vec2 p) {
        const float outer = std::fmax(std::fabs(p.x - 50.0f) - 40.0f, std::fabs(p.y - 50.0f) - 40.0f);
        const float hollow = std::fmax(std::fabs(p.x - 50.0f) - 20.0f, std::fabs(p.y - 30.0f) - 40.0f);
        return std::fmax(outer, -hollow);
    });
    World w;
    BodyDesc cup;
    cup.type = MotionType::Static;
    addPieces(w, cup, alphaPieces(c.image()), { 50, 100 }, 0.05f);
    BodyDesc bd;
    bd.position = { 0, 4, 0 };
    bd.locks = LockPlanar;
    const BodyId box = w.createBody(bd);
    ShapeDesc sd;
    sd.geometry = Box{ { 0.3f, 0.3f, 0.05f } };
    w.addShape(box, sd);
    for (int i = 0; i < 240; ++i) { w.step(kDt); }
    CHECK_NEAR(w.position(box).y, 1.8f, 0.002f);
    CHECK_NEAR(w.position(box).x, 0.0f, 0.01f);
    CHECK(!w.isAwake(box));
}

TEST(traced_ring_hangs_on_a_peg_through_its_hole) {
    const Canvas c = paint(100, 100, [](Vec2 p) { return ring(p, { 50, 50 }, 40, 24); });
    World w;
    BodyDesc peg;
    peg.type = MotionType::Static;
    peg.position = { 0, 2, 0 };
    const BodyId pb = w.createBody(peg);
    ShapeDesc ps;
    ps.geometry = Box{ { 0.1f, 0.1f, 0.05f } };
    w.addShape(pb, ps);
    BodyDesc bd;
    bd.position = { 0, 2.3f, 0 };
    bd.locks = LockPlanar;
    const BodyId ringBody = addPieces(w, bd, alphaPieces(c.image()), { 50, 50 }, 0.025f);
    float deepest = 0.0f;
    for (int i = 0; i < 300; ++i) {
        w.step(kDt);
        deepest = std::fmax(deepest, deepestContact(w));
    }
    CHECK_NEAR(w.position(ringBody).y, 1.5f, 0.01f);
    CHECK_NEAR(w.position(ringBody).x, 0.0f, 0.05f);
    CHECK(deepest < 0.01f);
}

TEST(traced_bodies_tumble_onto_ground_without_sinking) {
    const Canvas c = paint(100, 100, [](Vec2 p) {
        const float bar = std::fmax(std::fabs(p.x - 50.0f) - 5.0f, std::fabs(p.y - 50.0f) - 45.0f);
        return std::fmin(ring(p, { 50, 50 }, 30, 18), bar);
    });
    const std::vector<std::vector<Vec2>> pieces = alphaPieces(c.image());
    float deepest = 0.0f;
    for (int trial = 0; trial < 8; ++trial) {
        World w;
        addGround(w);
        BodyDesc bd;
        bd.position = { 0, 2.0f + 0.4f * (float)trial, 0 };
        bd.rotation = Quat::axisAngle({ 0, 0, 1 }, 0.7f * (float)trial);
        bd.angularVelocity = { 0, 0, (float)(trial % 3) - 1.0f };
        bd.locks = LockPlanar;
        const BodyId b = addPieces(w, bd, pieces, { 50, 50 }, 0.02f);
        for (int i = 0; i < 900; ++i) {
            w.step(kDt);
            deepest = std::fmax(deepest, deepestContact(w));
        }
        CHECK(!w.isAwake(b));
    }
    CHECK(deepest < 0.005f);
}

TEST(nearby_outlines_never_cross_when_simplified) {
    struct Disc {
        Vec2 center;
        float radius;
    };
    const std::vector<std::vector<Disc>> layouts{ { { { 28, 29 }, 5 }, { { 24, 21 }, 3 } }, { { { 25, 12 }, 5 }, { { 30, 9 }, 5 } } };
    for (size_t l = 0; l < layouts.size(); ++l) {
        const std::vector<Disc>& discs = layouts[l];
        const bool rings = l == 1;
        const Canvas c = paint(40, 40, [&](Vec2 p) {
            for (const Disc& d : discs) {
                const float r = length(p - d.center);
                if (r < d.radius && (!rings || r > d.radius - 2.0f)) { return -1.0f; }
            }
            return 1.0f;
        });
        for (const std::vector<Vec2>& contour : traceAlpha(c.image())) { CHECK(!crossesItself(simplifyContour(contour, 1.0f))); }
        const std::vector<std::vector<Vec2>> pieces = alphaPieces(c.image(), 128, 1.0f);
        for (const std::vector<Vec2>& p : pieces) { CHECK(convex(p)); }
        for (float y = 0.137f; y < 40.0f; y += 0.0391f) {
            for (float x = 0.211f; x < 40.0f; x += 0.0383f) { CHECK(coverCount(pieces, { x, y }) <= 1); }
        }
    }
}

// Four times the teeth should cost about four times as much, not sixteen.
TEST(comb_pieces_cost_grows_close_to_linearly) {
    auto comb = [](int teeth) {
        const int w = 4 * teeth + 4, h = 200;
        Canvas c{ w, h, std::vector<uint8_t>((size_t)w * h) };
        for (int y = 2; y < h - 2; ++y) {
            for (int x = 2; x < w - 2; ++x) {
                if (y > h - 20 || x % 4 < 2) { c.alpha[(size_t)y * w + x] = 255; }
            }
        }
        return c;
    };
    double ms[2] = { INFINITY, INFINITY };
    const int teeth[2] = { 500, 2000 };
    for (int k = 0; k < 2; ++k) {
        const Canvas c = comb(teeth[k]);
        for (int run = 0; run < 3; ++run) {
            const auto start = std::chrono::steady_clock::now();
            const std::vector<std::vector<Vec2>> pieces = alphaPieces(c.image());
            ms[k] = std::fmin(ms[k], std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
            CHECK((int)pieces.size() == 2 * teeth[k] + 4);
        }
    }
    CHECK(ms[1] < 8.0 * ms[0]);
}
