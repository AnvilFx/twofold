#include <algorithm>
#include <cfloat>
#include <utility>

#include "internal.h"

namespace twofold {

namespace {

constexpr int kMaxIterations = 32;
constexpr int kMaxAdvances = 64;

struct Vertex {
    Vec3 a, b, w;
};

// Simplex weights are found in double. Near flat features the simplex can be a sliver many times longer
// than its distance to the origin, and float weights tilt the normal enough to lose millimeters along a face.
struct Point {
    double x = 0.0, y = 0.0, z = 0.0;

    Point() = default;
    Point(const Vec3& v) : x(v.x), y(v.y), z(v.z) {}
    Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Point operator+(const Point& b) const { return { x + b.x, y + b.y, z + b.z }; }
    Point operator-(const Point& b) const { return { x - b.x, y - b.y, z - b.z }; }
    Point operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3 vec() const { return { (float)x, (float)y, (float)z }; }
};

double dot(const Point& a, const Point& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

Point cross(const Point& a, const Point& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

struct Simplex {
    Vertex v[4];
    double u[4] = {};
    int count = 0;
};

Vec3 support(const Convex& c, const Transform& xf, const Vec3& d) {
    const Vec3 local = xf.q.unrotate(d);
    if (c.cylinder) { return xf.apply(cylinderSupport(c.rim, c.halfHeight, local)); }
    int best = 0;
    float most = dot(c.points[0], local);
    for (int i = 1; i < c.count; ++i) {
        const float s = dot(c.points[i], local);
        if (s > most) {
            most = s;
            best = i;
        }
    }
    return xf.apply(c.points[best]);
}

void segmentWeights(const Point& a, const Point& b, double& ua, double& ub) {
    const Point e = b - a;
    const double ee = dot(e, e);
    const double t = ee > 0.0 ? std::clamp(-dot(a, e) / ee, 0.0, 1.0) : 0.0;
    ua = 1.0 - t;
    ub = t;
}

// Ericson's closest point on a triangle, as barycentric weights for the point nearest the origin.
void triangleWeights(const Point& a, const Point& b, const Point& c, double u[3]) {
    u[0] = u[1] = u[2] = 0.0;
    const Point ab = b - a, ac = c - a;
    const double d1 = -dot(ab, a), d2 = -dot(ac, a);
    if (d1 <= 0.0 && d2 <= 0.0) {
        u[0] = 1.0;
        return;
    }
    const double d3 = -dot(ab, b), d4 = -dot(ac, b);
    if (d3 >= 0.0 && d4 <= d3) {
        u[1] = 1.0;
        return;
    }
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        segmentWeights(a, b, u[0], u[1]);
        return;
    }
    const double d5 = -dot(ab, c), d6 = -dot(ac, c);
    if (d6 >= 0.0 && d5 <= d6) {
        u[2] = 1.0;
        return;
    }
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        segmentWeights(a, c, u[0], u[2]);
        return;
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
        segmentWeights(b, c, u[1], u[2]);
        return;
    }
    const double sum = va + vb + vc;
    if (sum > DBL_MIN) {
        u[0] = va / sum;
        u[1] = vb / sum;
        u[2] = vc / sum;
        return;
    }
    const Point* p[3] = { &a, &b, &c };
    double best = DBL_MAX;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        double ui, uj;
        segmentWeights(*p[i], *p[j], ui, uj);
        const Point q = *p[i] * ui + *p[j] * uj;
        if (dot(q, q) < best) {
            best = dot(q, q);
            u[0] = u[1] = u[2] = 0.0;
            u[i] = ui;
            u[j] = uj;
        }
    }
}

// False when the origin lies inside the tetrahedron. A face too flat to tell its inside from rounding is
// treated as facing the origin.
bool tetrahedronWeights(const Point w[4], double u[4]) {
    static const int faces[4][4] = { { 0, 1, 2, 3 }, { 0, 3, 1, 2 }, { 0, 2, 3, 1 }, { 1, 3, 2, 0 } };
    bool outside = false;
    double best = DBL_MAX;
    for (const int* f : faces) {
        const Point n = cross(w[f[1]] - w[f[0]], w[f[2]] - w[f[0]]);
        const Point apex = w[f[3]] - w[f[0]];
        const double height = dot(apex, n);
        if (-dot(w[f[0]], n) * height > 0.0 && height * height > 1e-8 * dot(n, n) * dot(apex, apex)) { continue; }
        outside = true;
        double t[3];
        triangleWeights(w[f[0]], w[f[1]], w[f[2]], t);
        const Point q = w[f[0]] * t[0] + w[f[1]] * t[1] + w[f[2]] * t[2];
        if (dot(q, q) < best) {
            best = dot(q, q);
            u[0] = u[1] = u[2] = u[3] = 0.0;
            for (int i = 0; i < 3; ++i) { u[f[i]] = t[i]; }
        }
    }
    return outside;
}

// Reduces the simplex to the vertices that carry the point nearest the origin. False on overlap.
bool solve(Simplex& s) {
    double u[4] = { 1.0, 0.0, 0.0, 0.0 };
    if (s.count == 2) {
        segmentWeights(s.v[0].w, s.v[1].w, u[0], u[1]);
    } else if (s.count == 3) {
        triangleWeights(s.v[0].w, s.v[1].w, s.v[2].w, u);
    } else if (s.count == 4) {
        const Point w[4] = { s.v[0].w, s.v[1].w, s.v[2].w, s.v[3].w };
        if (!tetrahedronWeights(w, u)) { return false; }
    }
    int kept = 0;
    for (int i = 0; i < s.count; ++i) {
        if (u[i] <= 0.0) { continue; }
        s.v[kept] = s.v[i];
        s.u[kept] = u[i];
        ++kept;
    }
    s.count = kept;
    return true;
}

Vec3 closest(const Simplex& s) {
    Point v;
    for (int i = 0; i < s.count; ++i) { v = v + Point(s.v[i].w) * s.u[i]; }
    return v.vec();
}

DistanceResult planeDistance(const Transform& plane, const Convex& c, const Transform& xf) {
    DistanceResult r;
    const Vec3 n = plane.q.rotate({ 0, 0, 1 });
    const Vec3 deepest = support(c, xf, -n);
    const float height = dot(n, deepest - plane.p);
    r.normal = n;
    r.pointA = deepest - n * height;
    r.pointB = deepest - n * c.radius;
    r.distance = height - c.radius;
    return r;
}

}

Convex convexOf(const Shape& s) {
    Convex c;
    if (s.kind == ShapeKind::HalfSpace) {
        c.halfSpace = true;
    } else if (s.kind == ShapeKind::Sphere) {
        c.points = &s.centroid;
        c.count = 1;
        c.radius = s.radius;
    } else if (s.kind == ShapeKind::Cylinder) {
        c.cylinder = true;
        c.rim = s.radius;
        c.halfHeight = s.halfHeight;
    } else {
        c.points = s.poly.vertices.data();
        c.count = (int)s.poly.vertices.size();
    }
    return c;
}

namespace {

enum class Verdict : uint8_t { Open, Within, Beyond };

// GJK from start, leaving the simplex and its point v nearest the origin. With a limit of zero or more, stops
// as soon as the gap between the cores is known to be within the limit or beyond it.
Verdict search(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, const Vec3& start,
               float limit, Simplex& s, Vec3& v, bool& overlap) {
    const bool bounded = limit >= 0.0f;
    if (lengthSquared(start) > 0.0f) {
        s.v[0].a = support(a, xa, start);
        s.v[0].b = support(b, xb, -start);
    } else {
        s.v[0].a = a.cylinder ? xa.apply({ a.rim, 0.0f, a.halfHeight }) : xa.apply(a.points[0]);
        s.v[0].b = b.cylinder ? xb.apply({ b.rim, 0.0f, b.halfHeight }) : xb.apply(b.points[0]);
    }
    s.v[0].w = s.v[0].b - s.v[0].a;
    s.u[0] = 1.0;
    s.count = 1;
    v = s.v[0].w;
    overlap = false;
    for (int iter = 0; iter < kMaxIterations; ++iter) {
        const float vv = dot(v, v);
        if (vv < 1.0e-14f) {
            overlap = true;
            return bounded ? Verdict::Within : Verdict::Open;
        }
        if (bounded && vv <= limit * limit) { return Verdict::Within; }
        Vertex n;
        n.a = support(a, xa, v);
        n.b = support(b, xb, -v);
        n.w = n.b - n.a;
        const float vw = dot(v, n.w);
        if (bounded && vw > 0.0f && vw * vw > limit * limit * vv) { return Verdict::Beyond; }
        if (vv - vw <= 1.0e-6f * vv) { break; }
        bool repeated = false;
        for (int i = 0; i < s.count; ++i) { repeated = repeated || lengthSquared(s.v[i].w - n.w) < 1.0e-14f; }
        if (repeated) { break; }

        const Simplex previous = s;
        s.v[s.count++] = n;
        if (!solve(s)) {
            overlap = true;
            return bounded ? Verdict::Within : Verdict::Open;
        }
        const Vec3 next = closest(s);
        if (dot(next, next) >= vv) {
            s = previous;
            break;
        }
        v = next;
    }
    return Verdict::Open;
}

}

bool apart(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, float limit) {
    if (a.halfSpace || b.halfSpace || limit < 0.0f) { return distance(a, xa, b, xb).distance > limit; }
    Simplex s;
    Vec3 v;
    bool overlap;
    const Verdict verdict = search(a, xa, b, xb, {}, limit + a.radius + b.radius, s, v, overlap);
    if (verdict != Verdict::Open) { return verdict == Verdict::Beyond; }
    return length(v) - a.radius - b.radius > limit;
}

DistanceResult distance(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, const Vec3& start) {
    if (a.halfSpace || b.halfSpace) {
        if (a.halfSpace && b.halfSpace) { return { xa.p, xb.p, { 0, 0, 1 }, FLT_MAX }; }
        if (a.halfSpace) { return planeDistance(xa, b, xb); }
        DistanceResult r = planeDistance(xb, a, xa);
        std::swap(r.pointA, r.pointB);
        r.normal = -r.normal;
        return r;
    }

    Simplex s;
    Vec3 v;
    bool overlap;
    search(a, xa, b, xb, start, -1.0f, s, v, overlap);

    Point pointA, pointB;
    for (int i = 0; i < s.count; ++i) {
        pointA = pointA + Point(s.v[i].a) * s.u[i];
        pointB = pointB + Point(s.v[i].b) * s.u[i];
    }
    DistanceResult r;
    r.pointA = pointA.vec();
    r.pointB = pointB.vec();
    const float gap = overlap ? 0.0f : length(v);
    r.normal = gap > 0.0f ? v * (1.0f / gap) : Vec3{ 0, 0, 1 };
    r.distance = gap - a.radius - b.radius;
    r.pointA += r.normal * a.radius;
    r.pointB -= r.normal * b.radius;
    return r;
}

Transform Sweep::at(float t) const {
    const Quat q = angle > 0.0f ? Quat::axisAngle(axis, angle * t) * q0 : q0;
    const Vec3 c = c0 + (c1 - c0) * t;
    return { c - q.rotate(localCenter), q };
}

Sweep makeSweep(const Vec3& localCenter, const Vec3& c0, const Quat& q0, const Vec3& c1, const Quat& q1) {
    Sweep s;
    s.localCenter = localCenter;
    s.c0 = c0;
    s.c1 = c1;
    s.q0 = q0;
    Quat d = normalize(q1 * q0.conjugate());
    if (d.w < 0.0f) { d = -d; }
    const float sine = length(d.vec());
    if (sine > 1.0e-7f) {
        s.axis = d.vec() * (1.0f / sine);
        s.angle = 2.0f * std::atan2(sine, d.w);
    }
    return s;
}

Impact timeOfImpact(const Convex& a, const Transform& localA, const Sweep& sweepA, float reachA, const Convex& b,
                    const Transform& localB, const Sweep& sweepB, float reachB, float target, float& t) {
    const float tolerance = 0.25f * target;
    const Vec3 approach = (sweepA.c1 - sweepA.c0) - (sweepB.c1 - sweepB.c0);
    const float spin = sweepA.angle * reachA + sweepB.angle * reachB;
    t = 0.0f;
    for (int iter = 0; iter < kMaxAdvances; ++iter) {
        const Transform xa = sweepA.at(t) * localA, xb = sweepB.at(t) * localB;
        const DistanceResult r = distance(a, xa, b, xb);
        if (r.distance <= target + tolerance) { return t > 0.0f ? Impact::Hit : Impact::Touching; }
        // GJK can stop above the true distance. The gap along its normal never is.
        float gap = r.distance;
        if (!a.halfSpace && !b.halfSpace) {
            gap = minf(gap, dot(r.normal, support(b, xb, -r.normal) - support(a, xa, r.normal)) - a.radius - b.radius);
        }
        // No point of either shape closes the gap along the normal faster than this.
        const float closing = dot(approach, r.normal) + spin;
        if (closing <= 0.0f) { return Impact::Miss; }
        if (gap <= target) { return t > 0.0f ? Impact::Hit : Impact::Touching; }
        t += (gap - target) / closing;
        if (t >= 1.0f) { return Impact::Miss; }
    }
    return Impact::Hit;
}

Vec3 convexSupport(const Convex& c, const Transform& xf, const Vec3& d) { return support(c, xf, d); }

}
