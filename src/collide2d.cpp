#include <algorithm>
#include <cfloat>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kAligned = 0.9999f;

struct Profile {
    bool circle = false;
    Vec2 center;
    float radius = 0.0f;
    std::vector<Vec2> poly;
    float zLo = 0.0f, zHi = 0.0f;
};

struct Planar {
    Vec2 normal;
    Vec2 points[2];
    float separations[2] = {};
    int count = 0;
};

Vec2 project(const Vec3& p, const Vec3& u, const Vec3& v) { return { dot(p, u), dot(p, v) }; }

void setDepth(Profile& p, float center, float a, float b) {
    p.zLo = center + std::fmin(a, b);
    p.zHi = center + std::fmax(a, b);
}

void fixWinding(std::vector<Vec2>& poly) {
    float area = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { area += cross(poly[i], poly[(i + 1) % poly.size()]); }
    if (area < 0.0f) { std::reverse(poly.begin(), poly.end()); }
}

bool makeProfile(const Shape& s, const Transform& xf, const Vec3& n, const Vec3& u, const Vec3& v, Profile& p) {
    const float c = dot(xf.p, n);
    if (const Box* bx = std::get_if<Box>(&s.geometry)) {
        const Vec3 axes[3] = { xf.q.rotate({ 1, 0, 0 }), xf.q.rotate({ 0, 1, 0 }), xf.q.rotate({ 0, 0, 1 }) };
        for (int i = 0; i < 3; ++i) {
            const float a = dot(axes[i], n);
            if (std::fabs(a) < kAligned) { continue; }
            const Vec3 aj = axes[(i + 1) % 3] * bx->halfExtents[(i + 1) % 3];
            const Vec3 ak = axes[(i + 2) % 3] * bx->halfExtents[(i + 2) % 3];
            p.poly = { project(xf.p + aj + ak, u, v), project(xf.p - aj + ak, u, v),
                       project(xf.p - aj - ak, u, v), project(xf.p + aj - ak, u, v) };
            fixWinding(p.poly);
            setDepth(p, c, -bx->halfExtents[i], bx->halfExtents[i]);
            return true;
        }
        return false;
    }
    const Vec3 zw = xf.q.rotate({ 0, 0, 1 });
    const float a = dot(zw, n);
    if (std::fabs(a) < kAligned) { return false; }
    if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        p.circle = true;
        p.center = project(xf.p, u, v);
        p.radius = cy->radius;
        setDepth(p, c, -cy->halfHeight, cy->halfHeight);
        return true;
    }
    if (const Prism* pr = std::get_if<Prism>(&s.geometry)) {
        if (pr->polygon.size() < 3) { return false; }
        p.poly.clear();
        for (const Vec2& q : pr->polygon) { p.poly.push_back(project(xf.apply({ q.x, q.y, 0.0f }), u, v)); }
        fixWinding(p.poly);
        setDepth(p, c, a * pr->zMin, a * pr->zMax);
        return true;
    }
    return false;
}

Vec2 edgeNormal(const std::vector<Vec2>& poly, size_t i) {
    const Vec2 e = poly[(i + 1) % poly.size()] - poly[i];
    return normalize(Vec2{ e.y, -e.x });
}

struct SegmentClosest {
    Vec2 a, b;
    float s = 0.0f, t = 0.0f;
};

SegmentClosest closestOnSegments(const Vec2& p1, const Vec2& q1, const Vec2& p2, const Vec2& q2) {
    const Vec2 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    const float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    SegmentClosest c;
    if (a <= FLT_EPSILON && e > FLT_EPSILON) {
        c.t = std::clamp(f / e, 0.0f, 1.0f);
    } else if (a > FLT_EPSILON) {
        const float k = dot(d1, r);
        if (e <= FLT_EPSILON) {
            c.s = std::clamp(-k / a, 0.0f, 1.0f);
        } else {
            const float b = dot(d1, d2), denom = a * e - b * b;
            c.s = denom > 0.0f ? std::clamp((b * f - k * e) / denom, 0.0f, 1.0f) : 0.0f;
            c.t = (b * c.s + f) / e;
            if (c.t < 0.0f) {
                c.t = 0.0f;
                c.s = std::clamp(-k / a, 0.0f, 1.0f);
            } else if (c.t > 1.0f) {
                c.t = 1.0f;
                c.s = std::clamp((b - k) / a, 0.0f, 1.0f);
            }
        }
    }
    c.a = p1 + d1 * c.s;
    c.b = p2 + d2 * c.t;
    return c;
}

float maxSeparation(const std::vector<Vec2>& a, const std::vector<Vec2>& b, int& edge) {
    float best = -FLT_MAX;
    edge = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        const Vec2 n = edgeNormal(a, i);
        float s = FLT_MAX;
        for (const Vec2& q : b) { s = std::fmin(s, dot(n, q - a[i])); }
        if (s > best) { best = s; edge = (int)i; }
    }
    return best;
}

void collidePolygons(const std::vector<Vec2>& a, const std::vector<Vec2>& b, float margin, float slop, Planar& out) {
    int edgeA = 0, edgeB = 0;
    const float sa = maxSeparation(a, b, edgeA);
    if (sa > margin) { return; }
    const float sb = maxSeparation(b, a, edgeB);
    if (sb > margin) { return; }

    const bool flip = sb > sa + 0.1f * slop;
    const std::vector<Vec2>& ref = flip ? b : a;
    const std::vector<Vec2>& inc = flip ? a : b;
    const int i1 = flip ? edgeB : edgeA;
    const Vec2 v11 = ref[i1];
    const Vec2 v12 = ref[(i1 + 1) % ref.size()];
    const Vec2 normal = edgeNormal(ref, i1);

    size_t i2 = 0;
    float minDot = FLT_MAX;
    for (size_t i = 0; i < inc.size(); ++i) {
        const float d = dot(normal, edgeNormal(inc, i));
        if (d < minDot) { minDot = d; i2 = i; }
    }
    const Vec2 v21 = inc[i2];
    const Vec2 v22 = inc[(i2 + 1) % inc.size()];

    // Separated polygons nearest at two vertices touch along the line between them.
    if (std::fmax(sa, sb) > 0.1f * slop) {
        const SegmentClosest c = closestOnSegments(v11, v12, v21, v22);
        if ((c.s == 0.0f || c.s == 1.0f) && (c.t == 0.0f || c.t == 1.0f)) {
            const float dist = length(c.b - c.a);
            if (dist > margin || dist <= 0.0f) { return; }
            const Vec2 n = (c.b - c.a) * (1.0f / dist);
            out.normal = flip ? -n : n;
            out.points[0] = (c.a + c.b) * 0.5f;
            out.separations[0] = dist;
            out.count = 1;
            return;
        }
    }

    const Vec2 tangent = normalize(v12 - v11);
    const float lower1 = 0.0f;
    const float upper1 = dot(v12 - v11, tangent);
    const float upper2 = dot(v21 - v11, tangent);
    const float lower2 = dot(v22 - v11, tangent);
    const float span = upper2 - lower2;

    out.normal = flip ? -normal : normal;
    // An incident edge entirely past the reference edge cannot be clipped, so its deepest vertex stands in.
    if (upper2 < lower1 || lower2 > upper1) {
        size_t deepest = 0;
        float s = FLT_MAX;
        for (size_t i = 0; i < inc.size(); ++i) {
            const float d = dot(inc[i] - v11, normal);
            if (d < s) { s = d; deepest = i; }
        }
        if (s > margin) { return; }
        out.points[0] = inc[deepest] - normal * (0.5f * s);
        out.separations[0] = s;
        out.count = 1;
        return;
    }

    Vec2 vLower = v22, vUpper = v21;
    if (span > FLT_EPSILON) {
        if (lower2 < lower1) { vLower = v22 + (v21 - v22) * ((lower1 - lower2) / span); }
        if (upper2 > upper1) { vUpper = v22 + (v21 - v22) * ((upper1 - lower2) / span); }
    }

    for (const Vec2& p : { vLower, vUpper }) {
        const float s = dot(p - v11, normal);
        if (s > margin) { continue; }
        out.points[out.count] = p - normal * (0.5f * s);
        out.separations[out.count] = s;
        ++out.count;
    }
}

// Normal points from the polygon to the circle.
void collidePolygonCircle(const std::vector<Vec2>& poly, const Vec2& c, float r, float margin, Planar& out) {
    size_t edge = 0;
    float sep = -FLT_MAX;
    for (size_t i = 0; i < poly.size(); ++i) {
        const float s = dot(edgeNormal(poly, i), c - poly[i]);
        if (s > sep) { sep = s; edge = i; }
    }
    if (sep - r > margin) { return; }

    const Vec2 v1 = poly[edge];
    const Vec2 v2 = poly[(edge + 1) % poly.size()];
    Vec2 surface, n;
    if (sep > FLT_EPSILON && dot(c - v1, v2 - v1) < 0.0f) {
        surface = v1;
    } else if (sep > FLT_EPSILON && dot(c - v2, v1 - v2) < 0.0f) {
        surface = v2;
    } else {
        n = edgeNormal(poly, edge);
        out.normal = n;
        out.points[0] = ((c - n * sep) + (c - n * r)) * 0.5f;
        out.separations[0] = sep - r;
        out.count = 1;
        return;
    }
    const float dist = length(c - surface);
    const float s = dist - r;
    if (s > margin || dist <= 0.0f) { return; }
    n = (c - surface) * (1.0f / dist);
    out.normal = n;
    out.points[0] = (surface + (c - n * r)) * 0.5f;
    out.separations[0] = s;
    out.count = 1;
}

void collideCircles(const Vec2& ca, float ra, const Vec2& cb, float rb, float margin, Planar& out) {
    const Vec2 d = cb - ca;
    const float dist = length(d);
    const float s = dist - ra - rb;
    if (s > margin) { return; }
    const Vec2 n = dist > 0.0f ? d * (1.0f / dist) : Vec2{ 1, 0 };
    out.normal = n;
    out.points[0] = ca + n * (ra + 0.5f * s);
    out.separations[0] = s;
    out.count = 1;
}

}

bool collidePlanar(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
                   const Vec3& n, float margin, const Tolerances& tol, Manifold& m) {
    m.count = 0;
    const Vec3 u = perpendicular(n);
    const Vec3 v = cross(n, u);
    Profile pa, pb;
    if (!makeProfile(a, xa, n, u, v, pa) || !makeProfile(b, xb, n, u, v, pb)) { return false; }

    const float lo = std::fmax(pa.zLo, pb.zLo);
    const float hi = std::fmin(pa.zHi, pb.zHi);
    if (lo > hi) { return true; }
    const float zMid = 0.5f * (lo + hi);

    Planar out;
    if (pa.circle && pb.circle) {
        collideCircles(pa.center, pa.radius, pb.center, pb.radius, margin, out);
    } else if (pa.circle) {
        collidePolygonCircle(pb.poly, pa.center, pa.radius, margin, out);
        out.normal = -out.normal;
    } else if (pb.circle) {
        collidePolygonCircle(pa.poly, pb.center, pb.radius, margin, out);
    } else {
        collidePolygons(pa.poly, pb.poly, margin, tol.linearSlop, out);
    }

    m.normal = u * out.normal.x + v * out.normal.y;
    for (int i = 0; i < out.count; ++i) {
        m.points[i].position = u * out.points[i].x + v * out.points[i].y + n * zMid;
        m.points[i].separation = out.separations[i];
    }
    m.count = out.count;
    return true;
}

}
