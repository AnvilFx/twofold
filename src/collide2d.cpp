#include <algorithm>
#include <cfloat>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kAligned = 0.9999f;
constexpr size_t kInlinePoints = 32;

struct Polygon {
    const Vec2* points = nullptr;
    size_t count = 0;

    size_t size() const { return count; }
    const Vec2& operator[](size_t i) const { return points[i]; }
    const Vec2* begin() const { return points; }
    const Vec2* end() const { return points + count; }
};

enum class ProfileKind : uint8_t { Polygon, Circle, HalfPlane };

struct Profile {
    ProfileKind kind = ProfileKind::Polygon;
    Vec2 center;
    float radius = 0.0f;
    // The half-plane holds the points p at depth z with dot(normal, p) <= offset - tilt * z.
    Vec2 normal;
    float offset = 0.0f, tilt = 0.0f;
    float zLo = 0.0f, zHi = 0.0f;
    size_t count = 0;
    Vec2 inlinePoints[kInlinePoints];
    std::vector<Vec2> spill;

    Vec2* resize(size_t n) {
        count = n;
        if (n <= kInlinePoints) { return inlinePoints; }
        spill.resize(n);
        return spill.data();
    }
    Polygon polygon() const { return { count <= kInlinePoints ? inlinePoints : spill.data(), count }; }
};

struct Planar {
    Vec2 normal;
    Vec2 points[2];
    float separations[2] = {};
    int count = 0;
};

Vec2 project(const Vec3& p, const Vec3& u, const Vec3& v) { return { dot(p, u), dot(p, v) }; }

void setDepth(Profile& p, float center, float a, float b) {
    p.zLo = center + minf(a, b);
    p.zHi = center + maxf(a, b);
}

void fixWinding(Vec2* poly, size_t n) {
    float area = 0.0f;
    for (size_t i = 0; i < n; ++i) { area += cross(poly[i], poly[(i + 1) % n]); }
    if (area < 0.0f) { std::reverse(poly, poly + n); }
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
            Vec2* q = p.resize(4);
            q[0] = project(xf.p + aj + ak, u, v);
            q[1] = project(xf.p - aj + ak, u, v);
            q[2] = project(xf.p - aj - ak, u, v);
            q[3] = project(xf.p + aj - ak, u, v);
            fixWinding(q, 4);
            setDepth(p, c, -bx->halfExtents[i], bx->halfExtents[i]);
            return true;
        }
        return false;
    }
    const Vec3 zw = xf.q.rotate({ 0, 0, 1 });
    const float a = dot(zw, n);
    // A half-space standing on the plane cuts every depth along nearly the same line.
    if (s.kind == ShapeKind::HalfSpace) {
        if (a * a > 1.0f - kAligned * kAligned) { return false; }
        const Vec2 m = project(zw, u, v);
        const float inv = 1.0f / length(m);
        p.kind = ProfileKind::HalfPlane;
        p.normal = m * inv;
        p.offset = dot(zw, xf.p) * inv;
        p.tilt = a * inv;
        p.zLo = -FLT_MAX;
        p.zHi = FLT_MAX;
        return true;
    }
    if (std::fabs(a) < kAligned) { return false; }
    if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        p.kind = ProfileKind::Circle;
        p.center = project(xf.p, u, v);
        p.radius = cy->radius;
        setDepth(p, c, -cy->halfHeight, cy->halfHeight);
        return true;
    }
    if (const Prism* pr = std::get_if<Prism>(&s.geometry)) {
        if (pr->polygon.size() < 3) { return false; }
        Vec2* q = p.resize(pr->polygon.size());
        for (size_t i = 0; i < pr->polygon.size(); ++i) {
            q[i] = project(xf.apply({ pr->polygon[i].x, pr->polygon[i].y, 0.0f }), u, v);
        }
        fixWinding(q, pr->polygon.size());
        setDepth(p, c, a * pr->zMin, a * pr->zMax);
        return true;
    }
    return false;
}

Vec2 edgeNormal(const Polygon& poly, size_t i) {
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

float maxSeparation(const Polygon& a, const Polygon& b, int& edge) {
    float best = -FLT_MAX;
    edge = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        const Vec2 n = edgeNormal(a, i);
        float s = FLT_MAX;
        for (const Vec2& q : b) { s = minf(s, dot(n, q - a[i])); }
        if (s > best) { best = s; edge = (int)i; }
    }
    return best;
}

void collidePolygons(const Polygon& a, const Polygon& b, float margin, float slop, Planar& out) {
    int edgeA = 0, edgeB = 0;
    const float sa = maxSeparation(a, b, edgeA);
    if (sa > margin) { return; }
    const float sb = maxSeparation(b, a, edgeB);
    if (sb > margin) { return; }

    const bool flip = sb > sa + 0.1f * slop;
    const Polygon& ref = flip ? b : a;
    const Polygon& inc = flip ? a : b;
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
    if (maxf(sa, sb) > 0.1f * slop) {
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
void collidePolygonCircle(const Polygon& poly, const Vec2& c, float r, float margin, Planar& out) {
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

// Normal points out of the half-plane. The polygon touches with the edge that faces it most.
void collideHalfPlanePolygon(const Vec2& normal, float offset, const Polygon& poly, float margin, Planar& out) {
    size_t edge = 0;
    float minDot = FLT_MAX;
    for (size_t i = 0; i < poly.size(); ++i) {
        const float d = dot(normal, edgeNormal(poly, i));
        if (d < minDot) { minDot = d; edge = i; }
    }
    out.normal = normal;
    for (const Vec2& p : { poly[(edge + 1) % poly.size()], poly[edge] }) {
        const float s = dot(normal, p) - offset;
        if (s > margin) { continue; }
        out.points[out.count] = p - normal * (0.5f * s);
        out.separations[out.count] = s;
        ++out.count;
    }
}

void collideHalfPlaneCircle(const Vec2& normal, float offset, const Vec2& c, float r, float margin, Planar& out) {
    const float sep = dot(normal, c) - offset;
    if (sep - r > margin) { return; }
    out.normal = normal;
    out.points[0] = ((c - normal * sep) + (c - normal * r)) * 0.5f;
    out.separations[0] = sep - r;
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
    if (pa.kind == ProfileKind::HalfPlane && pb.kind == ProfileKind::HalfPlane) { return false; }

    const float lo = maxf(pa.zLo, pb.zLo);
    const float hi = minf(pa.zHi, pb.zHi);
    if (lo > hi) { return true; }
    const float zMid = 0.5f * (lo + hi);

    Planar out;
    if (pa.kind == ProfileKind::HalfPlane || pb.kind == ProfileKind::HalfPlane) {
        const bool flip = pb.kind == ProfileKind::HalfPlane;
        const Profile& h = flip ? pb : pa;
        const Profile& o = flip ? pa : pb;
        const float offset = h.offset - h.tilt * zMid;
        if (o.kind == ProfileKind::Circle) {
            collideHalfPlaneCircle(h.normal, offset, o.center, o.radius, margin, out);
        } else {
            collideHalfPlanePolygon(h.normal, offset, o.polygon(), margin, out);
        }
        if (flip) { out.normal = -out.normal; }
    } else if (pa.kind == ProfileKind::Circle && pb.kind == ProfileKind::Circle) {
        collideCircles(pa.center, pa.radius, pb.center, pb.radius, margin, out);
    } else if (pa.kind == ProfileKind::Circle) {
        collidePolygonCircle(pb.polygon(), pa.center, pa.radius, margin, out);
        out.normal = -out.normal;
    } else if (pb.kind == ProfileKind::Circle) {
        collidePolygonCircle(pa.polygon(), pb.center, pb.radius, margin, out);
    } else {
        collidePolygons(pa.polygon(), pb.polygon(), margin, tol.linearSlop, out);
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
