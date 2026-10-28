#include <algorithm>
#include <cfloat>

#include "internal.h"
#include "twofold/hull.h"

namespace twofold {

namespace {

constexpr float kPi = 3.14159265358979f;

std::vector<Vec2> counterClockwise(std::vector<Vec2> poly) {
    float area = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { area += cross(poly[i], poly[(i + 1) % poly.size()]); }
    if (area < 0.0f) { std::reverse(poly.begin(), poly.end()); }
    return poly;
}

// Below this a box side or prism edge is treated as having no length.
float thinLimit(float size, float reach) { return std::max({ 1e-6f * size, 3.0f * FLT_EPSILON * reach, 1e-8f }); }

// Half thickness given to a flat shape, as computeHull gives a flat cloud.
float slabFor(float size, float limit) { return std::max({ 0.005f * size, 1e-4f, 10.0f * limit }); }

bool turnsLeft(const Vec2& a, const Vec2& b, const Vec2& c, float limit) {
    return cross(b - a, c - b) > limit * length(c - a);
}

// Counter-clockwise, winding once, with every corner clear of the line through its neighbors.
bool cleanOutline(const std::vector<Vec2>& poly, float limit) {
    const size_t n = poly.size();
    if (n < 3) { return false; }
    float turning = 0.0f;
    for (size_t i = 0; i < n; ++i) {
        const Vec2& a = poly[i];
        const Vec2& b = poly[(i + 1) % n];
        const Vec2& c = poly[(i + 2) % n];
        if (!turnsLeft(a, b, c, limit)) { return false; }
        turning += std::atan2(cross(b - a, c - b), dot(b - a, c - b));
    }
    return turning < 3.0f * kPi;
}

// Convex outline of the points, dropping repeated and nearly collinear ones. A segment or a single point
// becomes a thin rectangle, or a square around the center.
std::vector<Vec2> solidOutline(std::vector<Vec2> pts, float limit, float slab, const Vec2& c) {
    std::sort(pts.begin(), pts.end(), [](const Vec2& a, const Vec2& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
    std::vector<Vec2> out;
    for (int pass = 0; pass < 2 && !pts.empty(); ++pass) {
        const size_t base = out.size();
        for (const Vec2& p : pts) {
            while (out.size() >= base + 2 && !turnsLeft(out[out.size() - 2], out.back(), p, limit)) { out.pop_back(); }
            out.push_back(p);
        }
        out.pop_back();
        std::reverse(pts.begin(), pts.end());
    }
    if (out.size() >= 3) { return out; }
    if (out.size() == 2 && length(out[1] - out[0]) > limit) {
        const Vec2 d = normalize(out[1] - out[0]);
        const Vec2 side = Vec2{ -d.y, d.x } * slab;
        return { out[0] - side, out[1] - side, out[1] + side, out[0] + side };
    }
    return { c + Vec2{ -slab, -slab }, c + Vec2{ slab, -slab }, c + Vec2{ slab, slab }, c + Vec2{ -slab, slab } };
}

float solidSize(float s) { return std::isfinite(s) ? std::fabs(s) : 0.0f; }

// Negative half extents flip and non-finite ones become zero. A side with no length is thickened.
Vec3 solidBox(Vec3 h) {
    for (int k = 0; k < 3; ++k) { h[k] = std::isfinite(h[k]) ? std::fabs(h[k]) : 0.0f; }
    const float size = 2.0f * length(h);
    const float limit = thinLimit(size, h.x + h.y + h.z);
    for (int k = 0; k < 3; ++k) {
        if (2.0f * h[k] <= limit) { h[k] = slabFor(size, limit); }
    }
    return h;
}

// Cleaned like a box, so a sphere with no radius gets the slab of a box with no size.
float solidRadius(float r) {
    r = solidSize(r);
    const float limit = thinLimit(2.0f * r, r);
    return 2.0f * r <= limit ? slabFor(2.0f * r, limit) : r;
}

void solidCylinder(Cylinder& c) {
    c.radius = solidSize(c.radius);
    c.halfHeight = solidSize(c.halfHeight);
    const float size = 2.0f * std::sqrt(c.radius * c.radius + c.halfHeight * c.halfHeight);
    const float limit = thinLimit(size, c.radius + c.halfHeight);
    if (2.0f * c.radius <= limit) { c.radius = slabFor(size, limit); }
    if (2.0f * c.halfHeight <= limit) { c.halfHeight = slabFor(size, limit); }
}

// Non-finite points are dropped. A polygon that is not a clean convex outline is replaced by its convex
// outline, and a prism with no height is thickened.
void solidPrism(Prism& pr) {
    std::vector<Vec2> pts;
    for (const Vec2& p : pr.polygon) {
        if (std::isfinite(p.x) && std::isfinite(p.y)) { pts.push_back(p); }
    }
    if (!std::isfinite(pr.zMin)) { pr.zMin = 0.0f; }
    if (!std::isfinite(pr.zMax)) { pr.zMax = 0.0f; }
    if (pr.zMax < pr.zMin) { std::swap(pr.zMin, pr.zMax); }
    pts = counterClockwise(std::move(pts));

    Vec2 lo{}, hi{};
    if (!pts.empty()) { lo = hi = pts[0]; }
    for (const Vec2& p : pts) {
        lo = { std::min(lo.x, p.x), std::min(lo.y, p.y) };
        hi = { std::max(hi.x, p.x), std::max(hi.y, p.y) };
    }
    const float width = length(hi - lo);
    const float wide = std::max(std::fabs(lo.x), std::fabs(hi.x)) + std::max(std::fabs(lo.y), std::fabs(hi.y));
    const float height = pr.zMax - pr.zMin;
    const float size = std::sqrt(width * width + height * height);
    const float limit = thinLimit(size, wide + std::max(std::fabs(pr.zMin), std::fabs(pr.zMax)));
    const float slab = slabFor(size, limit);

    const float flat = thinLimit(width, wide);
    if (!cleanOutline(pts, flat)) { pts = solidOutline(std::move(pts), flat, slab, (lo + hi) * 0.5f); }
    pr.polygon = std::move(pts);
    if (height <= limit) {
        const float mid = 0.5f * (pr.zMin + pr.zMax);
        pr.zMin = mid - slab;
        pr.zMax = mid + slab;
    }
}

}

Vec3 cylinderReach(const Vec3& axis, float radius, float halfHeight) {
    Vec3 r;
    for (int i = 0; i < 3; ++i) {
        r[i] = halfHeight * std::fabs(axis[i]) + radius * std::sqrt(maxf(1.0f - axis[i] * axis[i], 0.0f));
    }
    return r;
}

Vec3 cylinderSupport(float radius, float halfHeight, const Vec3& d) {
    const float z = d.z < 0.0f ? -halfHeight : halfHeight;
    const float side = std::sqrt(d.x * d.x + d.y * d.y);
    if (side == 0.0f) { return { 0.0f, 0.0f, z }; }
    return { d.x * (radius / side), d.y * (radius / side), z };
}

void buildShape(Shape& s) {
    s.poly = Polyhedron{};
    s.radius = 0.0f;
    s.halfHeight = 0.0f;
    s.centroid = Vec3{};
    if (Sphere* sp = std::get_if<Sphere>(&s.geometry)) {
        s.kind = ShapeKind::Sphere;
        sp->radius = solidRadius(sp->radius);
        s.radius = sp->radius;
        s.extent = sp->radius;
    } else if (Box* bx = std::get_if<Box>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        bx->halfExtents = solidBox(bx->halfExtents);
        s.poly = makeBox(bx->halfExtents);
    } else if (Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        s.kind = ShapeKind::Cylinder;
        solidCylinder(*cy);
        s.radius = cy->radius;
        s.halfHeight = cy->halfHeight;
        s.extent = std::sqrt(cy->radius * cy->radius + cy->halfHeight * cy->halfHeight);
    } else if (Prism* pr = std::get_if<Prism>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        solidPrism(*pr);
        s.poly = makePrism(pr->polygon, pr->zMin, pr->zMax);
    } else if (const Hull* hu = std::get_if<Hull>(&s.geometry)) {
        HullMesh mesh = computeHull(hu->points, hu->maxVertices);
        if (mesh.vertices.empty()) { mesh = computeHull({ Vec3{} }, hu->maxVertices); }
        s.kind = ShapeKind::Convex;
        s.poly = makePolyhedron(std::move(mesh.vertices), mesh.faces);
    } else if (std::holds_alternative<Mesh>(s.geometry)) {
        s.kind = ShapeKind::Mesh;
        s.extent = 0.0f;
        if (const TriangleMesh::Data* m = meshOf(s)) {
            for (const Vec3& v : m->vertices) { s.extent = std::max(s.extent, length(v)); }
            s.centroid = (m->nodes[0].lo + m->nodes[0].hi) * 0.5f;
        }
    } else {
        s.kind = ShapeKind::HalfSpace;
        s.extent = 0.0f;
    }
    s.minExtent = s.kind == ShapeKind::Cylinder ? std::min(s.radius, s.halfHeight) : s.radius;
    if (s.kind == ShapeKind::Convex) {
        float e = 0.0f;
        for (const Vec3& v : s.poly.vertices) { e = std::max(e, length(v)); }
        s.extent = e;
        s.centroid = s.poly.centroid;
        s.minExtent = e;
        for (size_t i = 0; i < s.poly.faces.size(); ++i) {
            s.minExtent = std::min(s.minExtent, s.poly.offsets[i] - dot(s.poly.normals[i], s.centroid));
        }
    }
}

MassProperties shapeMass(const Shape& s) {
    MassProperties mp;
    if (const Sphere* sp = std::get_if<Sphere>(&s.geometry)) {
        const float r = sp->radius;
        mp.mass = s.density * (4.0f / 3.0f) * kPi * r * r * r;
        const float i = 0.4f * mp.mass * r * r;
        mp.inertia = Mat3::diagonal({ i, i, i });
    } else if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        const float r = cy->radius, h = 2.0f * cy->halfHeight;
        mp.mass = s.density * kPi * r * r * h;
        const float side = mp.mass * (3.0f * r * r + h * h) / 12.0f;
        mp.inertia = Mat3::diagonal({ side, side, 0.5f * mp.mass * r * r });
    } else if (s.kind == ShapeKind::Convex) {
        mp = polyhedronMass(s.poly, s.density);
    } else {
        return mp;
    }
    const Mat3 r = s.local.q.toMat3();
    mp.inertia = r * mp.inertia * transpose(r);
    mp.center = s.local.apply(mp.center);
    return mp;
}

void shapeBounds(const Shape& s, const Transform& xf, Vec3& lo, Vec3& hi) {
    if (s.kind == ShapeKind::HalfSpace) {
        const float big = 3.0e38f;
        lo = { -big, -big, -big };
        hi = { big, big, big };
        return;
    }
    if (s.kind == ShapeKind::Sphere) {
        const Vec3 r{ s.radius, s.radius, s.radius };
        lo = xf.p - r;
        hi = xf.p + r;
        return;
    }
    if (s.kind == ShapeKind::Cylinder) {
        const Vec3 r = cylinderReach(xf.q.rotate({ 0, 0, 1 }), s.radius, s.halfHeight);
        lo = xf.p - r;
        hi = xf.p + r;
        return;
    }
    if (s.kind == ShapeKind::Mesh) {
        const TriangleMesh::Data* m = meshOf(s);
        if (!m) {
            lo = hi = xf.p;
            return;
        }
        const Vec3 c = xf.apply((m->nodes[0].lo + m->nodes[0].hi) * 0.5f);
        const Vec3 e = (m->nodes[0].hi - m->nodes[0].lo) * 0.5f;
        const Mat3 r = xf.q.toMat3();
        Vec3 reach;
        for (int i = 0; i < 3; ++i) {
            reach[i] = std::fabs(r.r[i].x) * e.x + std::fabs(r.r[i].y) * e.y + std::fabs(r.r[i].z) * e.z;
        }
        lo = c - reach;
        hi = c + reach;
        return;
    }
    lo = hi = xf.apply(s.poly.vertices[0]);
    for (const Vec3& v : s.poly.vertices) {
        const Vec3 w = xf.apply(v);
        lo = minVec(lo, w);
        hi = maxVec(hi, w);
    }
}

}
