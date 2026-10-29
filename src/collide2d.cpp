#include <algorithm>
#include <cfloat>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kAligned = 0.9999f;
// How far the plane's normal may turn in a body's frame before its pose key and its shapes' profiles are
// made again.
constexpr float kProfileDrift = 1e-5f;

struct Polygon {
    const Vec2* points = nullptr;
    const Vec2* normals = nullptr;
    size_t count = 0;

    size_t size() const { return count; }
    const Vec2& operator[](size_t i) const { return points[i]; }
    const Vec2* begin() const { return points; }
    const Vec2* end() const { return points + count; }
};

struct Planar {
    Vec2 normal;
    Vec2 points[2];
    float separations[2] = {};
    int count = 0;
};

Vec2 project(const Vec3& p, const Vec3& u, const Vec3& v) { return { dot(p, u), dot(p, v) }; }

void setDepth(PlanarProfile& p, float center, float a, float b) {
    p.zLo = center + minf(a, b);
    p.zHi = center + maxf(a, b);
}

void fixWinding(std::vector<Vec2>& poly) {
    float area = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { area += cross(poly[i], poly[(i + 1) % poly.size()]); }
    if (area < 0.0f) { std::reverse(poly.begin(), poly.end()); }
}

// The profile of s seen along n, with xf, n, u and v all in one frame.
ProfileKind makeProfile(const Shape& s, const Transform& xf, const Vec3& n, const Vec3& u, const Vec3& v,
                        PlanarProfile& p) {
    const float c = dot(xf.p, n);
    std::vector<Vec2>& q = p.points;
    q.clear();
    if (const Box* bx = std::get_if<Box>(&s.geometry)) {
        const Vec3 axes[3] = { xf.q.rotate({ 1, 0, 0 }), xf.q.rotate({ 0, 1, 0 }), xf.q.rotate({ 0, 0, 1 }) };
        for (int i = 0; i < 3; ++i) {
            const float a = dot(axes[i], n);
            if (std::fabs(a) < kAligned) { continue; }
            const Vec3 aj = axes[(i + 1) % 3] * bx->halfExtents[(i + 1) % 3];
            const Vec3 ak = axes[(i + 2) % 3] * bx->halfExtents[(i + 2) % 3];
            q.push_back(project(xf.p + aj + ak, u, v));
            q.push_back(project(xf.p - aj + ak, u, v));
            q.push_back(project(xf.p - aj - ak, u, v));
            q.push_back(project(xf.p + aj - ak, u, v));
            setDepth(p, c, -bx->halfExtents[i], bx->halfExtents[i]);
            return ProfileKind::Polygon;
        }
        return ProfileKind::None;
    }
    const Vec3 zw = xf.q.rotate({ 0, 0, 1 });
    const float a = dot(zw, n);
    // A half-space standing on the plane cuts every depth along nearly the same line.
    if (s.kind == ShapeKind::HalfSpace) {
        if (a * a > 1.0f - kAligned * kAligned) { return ProfileKind::None; }
        const Vec2 m = project(zw, u, v);
        const float inv = 1.0f / length(m);
        p.normal = m * inv;
        p.offset = dot(zw, xf.p) * inv;
        p.tilt = a * inv;
        p.zLo = -FLT_MAX;
        p.zHi = FLT_MAX;
        return ProfileKind::HalfPlane;
    }
    if (std::fabs(a) < kAligned) { return ProfileKind::None; }
    if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        p.center = project(xf.p, u, v);
        p.radius = cy->radius;
        setDepth(p, c, -cy->halfHeight, cy->halfHeight);
        return ProfileKind::Circle;
    }
    if (const Prism* pr = std::get_if<Prism>(&s.geometry)) {
        if (pr->polygon.size() < 3) { return ProfileKind::None; }
        for (const Vec2& k : pr->polygon) { q.push_back(project(xf.apply({ k.x, k.y, 0.0f }), u, v)); }
        setDepth(p, c, a * pr->zMin, a * pr->zMax);
        return ProfileKind::Polygon;
    }
    return ProfileKind::None;
}

bool same(const Vec3& a, const Vec3& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }

// The profile of s in its body's frame for a plane whose normal there is n.
const PlanarProfile& profileOf(const Shape& s, const Vec3& n) {
    PlanarProfile& p = s.profile;
    if (same(n, p.key)) { return p; }
    p.key = n;
    p.u = perpendicular(n);
    p.v = cross(n, p.u);
    p.kind = makeProfile(s, s.local, n, p.u, p.v, p);
    p.normals.clear();
    if (p.kind == ProfileKind::Polygon) {
        fixWinding(p.points);
        Vec2 lo = p.points[0], hi = p.points[0];
        for (size_t i = 0; i < p.points.size(); ++i) {
            const Vec2 e = p.points[(i + 1) % p.points.size()] - p.points[i];
            p.normals.push_back(normalize(Vec2{ e.y, -e.x }));
            lo = { minf(lo.x, p.points[i].x), minf(lo.y, p.points[i].y) };
            hi = { maxf(hi.x, p.points[i].x), maxf(hi.y, p.points[i].y) };
        }
        p.center = (lo + hi) * 0.5f;
        p.extent = (hi - lo) * 0.5f;
    } else if (p.kind == ProfileKind::Circle) {
        p.extent = { p.radius, p.radius };
    }
    return p;
}

Polygon polygonOf(const PlanarProfile& p) { return { p.points.data(), p.normals.data(), p.points.size() }; }

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

// Stops once the polygons are known to be further apart than margin.
float maxSeparation(const Polygon& a, const Polygon& b, float margin, int& edge) {
    float best = -FLT_MAX;
    edge = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        const Vec2 n = a.normals[i];
        float s = FLT_MAX;
        for (const Vec2& q : b) {
            const float d = dot(n, q - a[i]);
            s = d < s ? d : s;
            if (s <= best) { break; }
        }
        if (s > best) {
            best = s;
            edge = (int)i;
            if (best > margin) { break; }
        }
    }
    return best;
}

void collidePolygons(const Polygon& a, const Polygon& b, float margin, float slop, Planar& out) {
    int edgeA = 0, edgeB = 0;
    const float sa = maxSeparation(a, b, margin, edgeA);
    if (sa > margin) { return; }
    const float sb = maxSeparation(b, a, margin, edgeB);
    if (sb > margin) { return; }

    const bool flip = sb > sa + 0.1f * slop;
    const Polygon& ref = flip ? b : a;
    const Polygon& inc = flip ? a : b;
    const int i1 = flip ? edgeB : edgeA;
    const Vec2 v11 = ref[i1];
    const Vec2 v12 = ref[(i1 + 1) % ref.size()];
    const Vec2 normal = ref.normals[i1];

    size_t i2 = 0;
    float minDot = FLT_MAX;
    for (size_t i = 0; i < inc.size(); ++i) {
        const float d = dot(normal, inc.normals[i]);
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
        const float s = dot(poly.normals[i], c - poly[i]);
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
        n = poly.normals[edge];
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
        const float d = dot(normal, poly.normals[i]);
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

bool collidePlanar(const Shape& a, const PlanarPose& bodyA, const Shape& b, const PlanarPose& bodyB, float margin,
                   const Tolerances& tol, Manifold& m) {
    m.count = 0;
    const PlanarProfile& pa = profileOf(a, bodyA.key);
    const PlanarProfile& pb = profileOf(b, bodyB.key);
    if (pa.kind == ProfileKind::None || pb.kind == ProfileKind::None) { return false; }
    if (pa.kind == ProfileKind::HalfPlane && pb.kind == ProfileKind::HalfPlane) { return false; }

    // The pair is worked out in a's 2D frame. b's profile comes over by one turn and shift, and its depths
    // by one lift.
    const Vec2 ta = bodyA.turn, tb = bodyB.turn;
    const Vec2 turn{ ta.x * tb.x + ta.y * tb.y, ta.x * tb.y - ta.y * tb.x };
    const Vec2 gap = bodyB.position - bodyA.position;
    const Vec2 shift{ ta.x * gap.x + ta.y * gap.y, ta.x * gap.y - ta.y * gap.x };
    const float lift = bodyB.depth - bodyA.depth;
    auto turned = [&](const Vec2& p) { return Vec2{ turn.x * p.x - turn.y * p.y, turn.y * p.x + turn.x * p.y }; };
    auto moved = [&](const Vec2& p) { return turned(p) + shift; };

    const bool halfB = pb.kind == ProfileKind::HalfPlane;
    const float lo = maxf(pa.zLo, halfB ? -FLT_MAX : pb.zLo + lift);
    const float hi = minf(pa.zHi, halfB ? FLT_MAX : pb.zHi + lift);
    if (lo > hi) { return true; }
    const float zMid = 0.5f * (lo + hi);
    // Profiles whose boxes are further apart than the margin can't touch.
    if (pa.kind != ProfileKind::HalfPlane && !halfB) {
        const Vec2 d = moved(pb.center) - pa.center;
        const float c = std::fabs(turn.x), s = std::fabs(turn.y);
        const Vec2 ea = pa.extent + Vec2{ margin, margin }, eb = pb.extent;
        if (std::fabs(d.x) > ea.x + c * eb.x + s * eb.y || std::fabs(d.y) > ea.y + s * eb.x + c * eb.y) { return true; }
        const Vec2 db{ turn.x * d.x + turn.y * d.y, turn.x * d.y - turn.y * d.x };
        if (std::fabs(db.x) > eb.x + c * ea.x + s * ea.y || std::fabs(db.y) > eb.y + s * ea.x + c * ea.y) { return true; }
    }

    thread_local std::vector<Vec2> points, normals;
    Polygon polyB;
    if (pb.kind == ProfileKind::Polygon) {
        points.resize(pb.points.size());
        normals.resize(pb.points.size());
        for (size_t i = 0; i < pb.points.size(); ++i) {
            points[i] = moved(pb.points[i]);
            normals[i] = turned(pb.normals[i]);
        }
        polyB = { points.data(), normals.data(), points.size() };
    }

    Planar out;
    if (pa.kind == ProfileKind::HalfPlane) {
        const float offset = pa.offset - pa.tilt * zMid;
        if (pb.kind == ProfileKind::Circle) {
            collideHalfPlaneCircle(pa.normal, offset, moved(pb.center), pb.radius, margin, out);
        } else {
            collideHalfPlanePolygon(pa.normal, offset, polyB, margin, out);
        }
    } else if (halfB) {
        const Vec2 normal = turned(pb.normal);
        const float offset = pb.offset + dot(normal, shift) - pb.tilt * (zMid - lift);
        if (pa.kind == ProfileKind::Circle) {
            collideHalfPlaneCircle(normal, offset, pa.center, pa.radius, margin, out);
        } else {
            collideHalfPlanePolygon(normal, offset, polygonOf(pa), margin, out);
        }
        out.normal = -out.normal;
    } else if (pa.kind == ProfileKind::Circle && pb.kind == ProfileKind::Circle) {
        collideCircles(pa.center, pa.radius, moved(pb.center), pb.radius, margin, out);
    } else if (pa.kind == ProfileKind::Circle) {
        collidePolygonCircle(polyB, pa.center, pa.radius, margin, out);
        out.normal = -out.normal;
    } else if (pb.kind == ProfileKind::Circle) {
        collidePolygonCircle(polygonOf(pa), moved(pb.center), pb.radius, margin, out);
    } else {
        collidePolygons(polygonOf(pa), polyB, margin, tol.linearSlop, out);
    }
    if (out.count == 0) { return true; }

    // The normal stays in the plane, so a contact never pushes layers apart in depth.
    auto toPlane = [&](const Vec2& p) { return Vec2{ ta.x * p.x - ta.y * p.y, ta.y * p.x + ta.x * p.y }; };
    const Vec2 normal = toPlane(out.normal);
    m.normal = bodyA.axisU * normal.x + bodyA.axisV * normal.y;
    const Vec3 lifted = bodyA.normal * (bodyA.depth + zMid);
    for (int i = 0; i < out.count; ++i) {
        const Vec2 p = toPlane(out.points[i]) + bodyA.position;
        m.points[i].position = bodyA.axisU * p.x + bodyA.axisV * p.y + lifted;
        m.points[i].separation = out.separations[i];
    }
    m.count = out.count;
    return true;
}

void placeInPlane(PlanarPose& pose, const Transform& xf, const Vec3& n) {
    const Vec3 key = xf.q.unrotate(n);
    if (!same(pose.normal, n) || lengthSquared(key - pose.key) > kProfileDrift * kProfileDrift) {
        pose.normal = n;
        pose.key = key;
        pose.u = perpendicular(key);
        pose.v = cross(key, pose.u);
        pose.axisU = perpendicular(n);
        pose.axisV = cross(n, pose.axisU);
    }
    const Vec3 u = xf.q.rotate(pose.u);
    pose.turn = normalize(Vec2{ dot(u, pose.axisU), dot(u, pose.axisV) });
    pose.position = { dot(xf.p, pose.axisU), dot(xf.p, pose.axisV) };
    pose.depth = dot(xf.p, n);
}

}
