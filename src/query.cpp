#include <cfloat>
#include <cmath>
#include <utility>
#include <vector>

#include "internal.h"

namespace twofold {

namespace {

bool raySphere(float radius, const Vec3& o, const Vec3& d, float maxFraction, float& fraction, Vec3& normal) {
    const float a = dot(d, d);
    const float b = dot(o, d);
    const float c = dot(o, o) - radius * radius;
    if (c < 0.0f || a == 0.0f) { return false; }
    const float disc = b * b - a * c;
    if (disc < 0.0f) { return false; }
    const float t = (-b - std::sqrt(disc)) / a;
    if (t < 0.0f || t > maxFraction) { return false; }
    fraction = t;
    normal = normalize(o + d * t);
    return true;
}

bool rayHalfSpace(const Vec3& o, const Vec3& d, float maxFraction, float& fraction, Vec3& normal) {
    if (o.z <= 0.0f || d.z >= 0.0f) { return false; }
    const float t = -o.z / d.z;
    if (t > maxFraction) { return false; }
    fraction = t;
    normal = { 0, 0, 1 };
    return true;
}

bool rayCylinder(float radius, float halfHeight, const Vec3& o, const Vec3& d, float maxFraction, float& fraction,
                 Vec3& normal) {
    const float rr = radius * radius;
    const float c = o.x * o.x + o.y * o.y - rr;
    if (c <= 0.0f && std::fabs(o.z) <= halfHeight) { return false; }
    float best = maxFraction;
    bool hit = false;
    if (std::fabs(o.z) > halfHeight && o.z * d.z < 0.0f) {
        const float cap = o.z > 0.0f ? halfHeight : -halfHeight;
        const float t = (cap - o.z) / d.z;
        const float x = o.x + d.x * t, y = o.y + d.y * t;
        if (t <= best && x * x + y * y <= rr) {
            best = t;
            normal = { 0.0f, 0.0f, o.z > 0.0f ? 1.0f : -1.0f };
            hit = true;
        }
    }
    const float a = d.x * d.x + d.y * d.y;
    const float b = o.x * d.x + o.y * d.y;
    if (c > 0.0f && a > 0.0f && b < 0.0f) {
        const float disc = b * b - a * c;
        if (disc >= 0.0f) {
            const float t = (-b - std::sqrt(disc)) / a;
            if (t <= best && std::fabs(o.z + d.z * t) <= halfHeight) {
                best = t;
                normal = Vec3{ o.x + d.x * t, o.y + d.y * t, 0.0f } * (1.0f / radius);
                hit = true;
            }
        }
    }
    if (hit) { fraction = best; }
    return hit;
}

bool rayPolyhedron(const Polyhedron& poly, const Vec3& o, const Vec3& d, float maxFraction, float& fraction, Vec3& normal) {
    float lower = 0.0f, upper = maxFraction;
    int index = -1;
    for (int i = 0; i < (int)poly.faces.size(); ++i) {
        const float num = poly.offsets[i] - dot(poly.normals[i], o);
        const float den = dot(poly.normals[i], d);
        if (den == 0.0f) {
            if (num < 0.0f) { return false; }
        } else if (den < 0.0f && num < lower * den) {
            lower = num / den;
            index = i;
        } else if (den > 0.0f && num < upper * den) {
            upper = num / den;
        }
        if (upper < lower) { return false; }
    }
    if (index < 0) { return false; }
    fraction = lower;
    normal = poly.normals[index];
    return true;
}

constexpr int kMaxSweepSteps = 96;

// GJK can run out of iterations on a cylinder and stop above the distance, by an amount that depends on which
// shape comes first. Both orders give an upper bound, so near touch the smaller one is kept.
DistanceResult queryDistance(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, float touch) {
    DistanceResult r = distance(a, xa, b, xb);
    if (r.distance <= touch || r.distance > 100.0f * touch || (!a.cylinder && !b.cylinder)) { return r; }
    const DistanceResult s = distance(b, xb, a, xa);
    if (s.distance < r.distance) { r = { s.pointB, s.pointA, -s.normal, s.distance }; }
    return r;
}

// The first fraction up to maxFraction at which a, moving from xa along d, comes within touch of b by GJK
// distance, the measure the overlap queries use. Each advance is the larger of the gap along the GJK normal
// over the closing speed, which never passes the first touch, and the GJK distance over the speed, which GJK
// overestimates less than it tilts its normal on curved shapes. A step that lands inside is bisected back.
Impact sweepConvex(const Convex& a, const Transform& xa, const Vec3& d, float maxFraction, const Convex& b,
                   const Transform& xb, float touch, float& fraction) {
    const float speed = length(d);
    const float target = 0.5f * touch;
    float t = 0.0f, clear = 0.0f, inside = -1.0f;
    for (int step = 0; step < kMaxSweepSteps; ++step) {
        const Transform x{ xa.p + d * t, xa.q };
        const DistanceResult r = queryDistance(a, x, b, xb, touch);
        if (step == 0 && r.distance <= touch) { return Impact::Touching; }
        if (r.distance <= touch && (r.distance > 0.0f || (t - clear) * speed <= target)) {
            fraction = t;
            return Impact::Hit;
        }
        if (r.distance <= touch) {
            inside = t;
        } else {
            clear = t;
            if (t >= maxFraction) { return Impact::Miss; }
        }
        if (inside >= 0.0f) {
            t = 0.5f * (clear + inside);
            continue;
        }
        const float closing = dot(d, r.normal);
        if (closing <= 0.0f) { return Impact::Miss; }
        float gap = r.distance;
        if (!a.halfSpace && !b.halfSpace) {
            gap = minf(gap, dot(r.normal, convexSupport(b, xb, -r.normal) - convexSupport(a, x, r.normal)) - a.radius -
                                b.radius);
        }
        t = minf(t + maxf((gap - target) / closing, (r.distance - target) / speed), maxFraction);
    }
    if (inside < 0.0f) { return Impact::Miss; }
    fraction = inside;
    return Impact::Hit;
}

// The point of b nearest a and the normal there facing a. False when they overlap.
bool nearest(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, Vec3& point, Vec3& normal) {
    const DistanceResult r = distance(a, xa, b, xb);
    if (r.distance <= 0.0f) { return false; }
    point = r.pointB;
    normal = -r.normal;
    return true;
}

// The deepest contact point between the shape and q, midway across the gap, and its normal pushing q out.
// Returns its separation, or FLT_MAX with a zero normal when the narrowphase finds none.
float deepestContact(const Shape& s, const Transform& xf, const Shape& q, const Transform& xq, const Tolerances& tol,
                     Vec3& point, Vec3& normal) {
    point = xq.p;
    normal = Vec3{};
    float deepest = FLT_MAX;
    auto take = [&](const Manifold& m) {
        for (int k = 0; k < m.count; ++k) {
            if (m.points[k].separation >= deepest) { continue; }
            deepest = m.points[k].separation;
            point = m.points[k].position;
            normal = m.normal;
        }
    };
    const Motion rest{};
    if (s.kind == ShapeKind::Mesh) {
        std::vector<Manifold> found;
        collideMesh(s, xf, q, xq, tol.speculative, tol, rest, found);
        for (const Manifold& m : found) { take(m); }
    } else {
        Manifold m;
        collide(s, xf, q, xq, tol.speculative, tol, rest, m);
        take(m);
    }
    return deepest;
}

}

bool rayShape(const Shape& s, const Transform& xf, const Vec3& origin, const Vec3& d, float maxFraction,
              float& fraction, Vec3& normal) {
    const Vec3 o = xf.applyInverse(origin);
    const Vec3 dl = xf.q.unrotate(d);
    bool hit = false;
    switch (s.kind) {
        case ShapeKind::Sphere: hit = raySphere(s.radius, o, dl, maxFraction, fraction, normal); break;
        case ShapeKind::HalfSpace: hit = rayHalfSpace(o, dl, maxFraction, fraction, normal); break;
        case ShapeKind::Convex: hit = rayPolyhedron(s.poly, o, dl, maxFraction, fraction, normal); break;
        case ShapeKind::Cylinder: hit = rayCylinder(s.radius, s.halfHeight, o, dl, maxFraction, fraction, normal); break;
        case ShapeKind::Mesh:
            if (const TriangleMesh::Data* m = meshOf(s)) { hit = rayMesh(*m, o, dl, maxFraction, fraction, normal); }
            break;
    }
    if (hit) { normal = xf.q.rotate(normal); }
    return hit;
}

bool segmentHitsBox(const Vec3& lo, const Vec3& hi, const Vec3& origin, const Vec3& d, float maxFraction) {
    float lower = 0.0f, upper = maxFraction;
    for (int i = 0; i < 3; ++i) {
        if (d[i] == 0.0f) {
            if (origin[i] < lo[i] || origin[i] > hi[i]) { return false; }
            continue;
        }
        const float inv = 1.0f / d[i];
        float t0 = (lo[i] - origin[i]) * inv, t1 = (hi[i] - origin[i]) * inv;
        if (t0 > t1) { std::swap(t0, t1); }
        lower = maxf(lower, t0);
        upper = minf(upper, t1);
        if (lower > upper) { return false; }
    }
    return true;
}

bool shapeContains(const Shape& s, const Transform& xf, const Vec3& point, float tolerance) {
    const Vec3 p = xf.applyInverse(point);
    switch (s.kind) {
        case ShapeKind::Sphere: return length(p) <= s.radius + tolerance;
        case ShapeKind::HalfSpace: return p.z <= tolerance;
        case ShapeKind::Cylinder: {
            const float r = s.radius + tolerance;
            return std::fabs(p.z) <= s.halfHeight + tolerance && p.x * p.x + p.y * p.y <= r * r;
        }
        case ShapeKind::Convex:
            for (size_t i = 0; i < s.poly.faces.size(); ++i) {
                if (dot(s.poly.normals[i], p) - s.poly.offsets[i] > tolerance) { return false; }
            }
            return true;
        case ShapeKind::Mesh: {
            const TriangleMesh::Data* m = meshOf(s);
            if (!m) { return false; }
            const Convex single{ &p, 1, 0.0f, false };
            const Vec3 reach{ tolerance, tolerance, tolerance };
            bool inside = false;
            m->query(p - reach, p + reach, [&](int t) {
                if (inside) { return; }
                Vec3 c[3];
                m->corners(t, c);
                inside = !apart({ c, 3, 0.0f, false }, Transform{}, single, Transform{}, tolerance);
            });
            return inside;
        }
    }
    return false;
}

bool shapesOverlap(const Shape& s, const Transform& xf, const Shape& q, const Transform& xq, float tolerance) {
    const Convex cq = convexOf(q);
    if (s.kind != ShapeKind::Mesh) { return queryDistance(cq, xq, convexOf(s), xf, tolerance).distance <= tolerance; }
    const TriangleMesh::Data* m = meshOf(s);
    if (!m) { return false; }
    Vec3 lo, hi, meshLo, meshHi;
    shapeBounds(q, xq, lo, hi);
    const Vec3 reach{ tolerance, tolerance, tolerance };
    boundsInMesh(xf, lo - reach, hi + reach, meshLo, meshHi);
    bool found = false;
    m->query(meshLo, meshHi, [&](int t) {
        if (found) { return; }
        Vec3 c[3];
        m->corners(t, c);
        const Vec3 triLo = minVec(minVec(c[0], c[1]), c[2]), triHi = maxVec(maxVec(c[0], c[1]), c[2]);
        if (triLo.x > meshHi.x || triHi.x < meshLo.x || triLo.y > meshHi.y || triHi.y < meshLo.y ||
            triLo.z > meshHi.z || triHi.z < meshLo.z) {
            return;
        }
        found = queryDistance(cq, xq, { c, 3, 0.0f, false }, xf, tolerance).distance <= tolerance;
    });
    return found;
}

Impact castShape(const Shape& s, const Transform& xf, const Shape& q, const Transform& xq, const Vec3& d, float touch,
                 const Tolerances& tol, float& fraction, Vec3& point, Vec3& normal) {
    const Convex cq = convexOf(q);
    Impact result = Impact::Miss;
    float best = fraction;
    Vec3 tri[3];
    if (s.kind != ShapeKind::Mesh) {
        result = sweepConvex(cq, xq, d, best, convexOf(s), xf, touch, best);
    } else if (const TriangleMesh::Data* m = meshOf(s)) {
        // The bounds of q in the mesh's frame, slid through the triangle tree.
        Vec3 lo, hi, meshLo, meshHi;
        shapeBounds(q, xq, lo, hi);
        boundsInMesh(xf, lo, hi, meshLo, meshHi);
        const Vec3 center = (meshLo + meshHi) * 0.5f;
        const Vec3 reach = (meshHi - meshLo) * 0.5f + Vec3{ touch, touch, touch };
        const Vec3 dm = xf.q.unrotate(d);
        int stack[64];
        int top = 0;
        stack[top++] = 0;
        while (top > 0 && result != Impact::Touching) {
            const TriangleMesh::Data::Node& n = m->nodes[stack[--top]];
            if (!segmentHitsBox(n.lo - reach, n.hi + reach, center, dm, best)) { continue; }
            if (n.count == 0) {
                stack[top++] = n.first + 1;
                stack[top++] = n.first;
                continue;
            }
            for (int i = n.first; i < n.first + n.count && result != Impact::Touching; ++i) {
                Vec3 c[3];
                m->corners(i, c);
                const Impact r = sweepConvex(cq, xq, d, best, { c, 3, 0.0f, false }, xf, touch, best);
                if (r == Impact::Miss) { continue; }
                result = r;
                for (int k = 0; k < 3; ++k) { tri[k] = c[k]; }
            }
        }
    }
    if (result == Impact::Miss) { return result; }
    if (result == Impact::Touching) {
        fraction = 0.0f;
        deepestContact(s, xf, q, xq, tol, point, normal);
        return result;
    }
    fraction = best;
    // A normal across a gap this small is mostly rounding, so the contact is found a little before the hit.
    const Transform before{ xq.p + d * maxf(0.0f, best - 0.1f * tol.linearSlop / length(d)), xq.q };
    const float separation = deepestContact(s, xf, q, before, tol, point, normal);
    if (separation < FLT_MAX) {
        point -= normal * (0.5f * separation);
    } else if (!nearest(cq, before, s.kind == ShapeKind::Mesh ? Convex{ tri, 3, 0.0f, false } : convexOf(s), xf, point,
                        normal)) {
        point = xq.p + d * best;
        normal = Vec3{};
    }
    return result;
}

}
