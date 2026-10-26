#include <cfloat>

#include "internal.h"

namespace twofold {

namespace {

struct WorldPoly {
    const Polyhedron* src = nullptr;
    std::vector<Vec3> v;
    std::vector<Vec3> n;
    std::vector<float> d;
    Vec3 centroid;
};

WorldPoly toWorld(const Polyhedron& p, const Transform& xf) {
    WorldPoly w;
    w.src = &p;
    w.v.reserve(p.vertices.size());
    for (const Vec3& v : p.vertices) { w.v.push_back(xf.apply(v)); }
    w.n.reserve(p.faces.size());
    w.d.reserve(p.faces.size());
    for (const Polyhedron::Face& f : p.faces) {
        const Vec3 n = xf.q.rotate(f.normal);
        w.n.push_back(n);
        w.d.push_back(dot(n, w.v[p.indices[f.first]]));
    }
    w.centroid = xf.apply(p.centroid);
    return w;
}

struct FaceQuery {
    int index = -1;
    float separation = -FLT_MAX;
};

struct EdgeQuery {
    int a = -1, b = -1;
    float separation = -FLT_MAX;
    Vec3 axis;
};

FaceQuery queryFaces(const WorldPoly& a, const WorldPoly& b) {
    FaceQuery q;
    for (size_t i = 0; i < a.n.size(); ++i) {
        float s = FLT_MAX;
        for (const Vec3& v : b.v) { s = std::fmin(s, dot(a.n[i], v) - a.d[i]); }
        if (s > q.separation) { q.separation = s; q.index = (int)i; }
    }
    return q;
}

// Whether the Gauss map arcs a-b and c-d intersect.
bool isMinkowskiFace(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    const Vec3 bxa = cross(b, a);
    const Vec3 dxc = cross(d, c);
    const float cba = dot(c, bxa), dba = dot(d, bxa);
    const float adc = dot(a, dxc), bdc = dot(b, dxc);
    return cba * dba < 0.0f && adc * bdc < 0.0f && cba * bdc > 0.0f;
}

EdgeQuery queryEdges(const WorldPoly& a, const WorldPoly& b) {
    EdgeQuery q;
    const std::vector<Polyhedron::Edge>& ea = a.src->edges;
    const std::vector<Polyhedron::Edge>& eb = b.src->edges;
    for (size_t i = 0; i < ea.size(); ++i) {
        const Vec3 pa = a.v[ea[i].v0];
        const Vec3 da = a.v[ea[i].v1] - pa;
        const Vec3 ua = a.n[ea[i].f0], va = a.n[ea[i].f1];
        for (size_t j = 0; j < eb.size(); ++j) {
            const Vec3 ub = b.n[eb[j].f0], vb = b.n[eb[j].f1];
            if (!isMinkowskiFace(ua, va, -ub, -vb)) { continue; }
            const Vec3 pb = b.v[eb[j].v0];
            const Vec3 db = b.v[eb[j].v1] - pb;
            Vec3 axis = cross(da, db);
            const float len = length(axis);
            if (len < 1e-5f * length(da) * length(db)) { continue; }
            axis = axis / len;
            if (dot(axis, pa - a.centroid) < 0.0f) { axis = -axis; }
            const float s = dot(axis, pb - pa);
            if (s > q.separation) {
                q.separation = s;
                q.a = (int)i;
                q.b = (int)j;
                q.axis = axis;
            }
        }
    }
    return q;
}

void clipPolygon(std::vector<Vec3>& poly, const Vec3& n, const Vec3& origin) {
    if (poly.empty()) { return; }
    std::vector<Vec3> out;
    out.reserve(poly.size() + 1);
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vec3& p = poly[i];
        const Vec3& q = poly[(i + 1) % poly.size()];
        const float dp = dot(n, p - origin);
        const float dq = dot(n, q - origin);
        if (dp <= 0.0f) { out.push_back(p); }
        if ((dp <= 0.0f) != (dq <= 0.0f)) { out.push_back(p + (q - p) * (dp / (dp - dq))); }
    }
    poly.swap(out);
}

void faceContact(const WorldPoly& ref, int refFace, const WorldPoly& inc, bool flip, float margin, Manifold& m) {
    const Vec3 refN = ref.n[refFace];
    const float refD = ref.d[refFace];

    int incFace = 0;
    float minDot = FLT_MAX;
    for (size_t i = 0; i < inc.n.size(); ++i) {
        const float d = dot(inc.n[i], refN);
        if (d < minDot) { minDot = d; incFace = (int)i; }
    }

    const Polyhedron::Face& fi = inc.src->faces[incFace];
    std::vector<Vec3> poly;
    for (int k = 0; k < fi.count; ++k) { poly.push_back(inc.v[inc.src->indices[fi.first + k]]); }

    const Polyhedron::Face& fr = ref.src->faces[refFace];
    for (int k = 0; k < fr.count; ++k) {
        const Vec3 a = ref.v[ref.src->indices[fr.first + k]];
        const Vec3 b = ref.v[ref.src->indices[fr.first + (k + 1) % fr.count]];
        clipPolygon(poly, cross(b - a, refN), a);
    }

    std::vector<ManifoldPoint> pts;
    for (const Vec3& p : poly) {
        const float s = dot(refN, p) - refD;
        if (s <= margin) { pts.push_back({ p - refN * (0.5f * s), s }); }
    }
    reduceManifold(pts, flip ? -refN : refN, m);
}

void closestOnSegments(const Vec3& p1, const Vec3& q1, const Vec3& p2, const Vec3& q2, Vec3& c1, Vec3& c2) {
    const Vec3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    const float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    float s = 0.0f, t = 0.0f;
    if (a <= FLT_EPSILON && e <= FLT_EPSILON) {
        c1 = p1;
        c2 = p2;
        return;
    }
    if (a <= FLT_EPSILON) {
        t = std::fmin(std::fmax(f / e, 0.0f), 1.0f);
    } else {
        const float c = dot(d1, r);
        if (e <= FLT_EPSILON) {
            s = std::fmin(std::fmax(-c / a, 0.0f), 1.0f);
        } else {
            const float b = dot(d1, d2);
            const float denom = a * e - b * b;
            s = denom != 0.0f ? std::fmin(std::fmax((b * f - c * e) / denom, 0.0f), 1.0f) : 0.0f;
            t = (b * s + f) / e;
            if (t < 0.0f) {
                t = 0.0f;
                s = std::fmin(std::fmax(-c / a, 0.0f), 1.0f);
            } else if (t > 1.0f) {
                t = 1.0f;
                s = std::fmin(std::fmax((b - c) / a, 0.0f), 1.0f);
            }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
}

void collideConvex(const Shape& sa, const Transform& xa, const Shape& sb, const Transform& xb,
                   float margin, const Tolerances& tol, Manifold& m) {
    const WorldPoly a = toWorld(sa.poly, xa);
    const WorldPoly b = toWorld(sb.poly, xb);

    const FaceQuery fa = queryFaces(a, b);
    if (fa.separation > margin) { return; }
    const FaceQuery fb = queryFaces(b, a);
    if (fb.separation > margin) { return; }
    const EdgeQuery eq = queryEdges(a, b);
    if (eq.separation > margin) { return; }

    const float absTol = 0.5f * tol.linearSlop;
    if (eq.a >= 0 && eq.separation > 0.9f * std::fmax(fa.separation, fb.separation) + absTol) {
        const Polyhedron::Edge& ea = sa.poly.edges[eq.a];
        const Polyhedron::Edge& eb = sb.poly.edges[eq.b];
        Vec3 c1, c2;
        closestOnSegments(a.v[ea.v0], a.v[ea.v1], b.v[eb.v0], b.v[eb.v1], c1, c2);
        m.normal = eq.axis;
        m.points[0] = { (c1 + c2) * 0.5f, eq.separation };
        m.count = 1;
    } else if (fb.separation > 0.98f * fa.separation + absTol) {
        faceContact(b, fb.index, a, true, margin, m);
    } else {
        faceContact(a, fa.index, b, false, margin, m);
    }
}

// Normal points from the polyhedron to the sphere.
void collideSpherePoly(const Shape& poly, const Transform& xp, const Vec3& c, float r, float margin, Manifold& m) {
    const WorldPoly p = toWorld(poly.poly, xp);
    int best = -1;
    float bestSep = -FLT_MAX;
    bool outside = false;
    for (size_t i = 0; i < p.n.size(); ++i) {
        const float s = dot(p.n[i], c) - p.d[i];
        if (s > r + margin) { return; }
        if (s > 0.0f) { outside = true; }
        if (s > bestSep) { bestSep = s; best = (int)i; }
    }
    if (best < 0) { return; }

    if (!outside) {
        const Vec3 n = p.n[best];
        const float sep = bestSep - r;
        m.normal = n;
        m.points[0] = { c - n * (0.5f * (bestSep + r)), sep };
        m.count = 1;
        return;
    }

    Vec3 closest;
    float closestSq = FLT_MAX;
    for (size_t i = 0; i < p.n.size(); ++i) {
        const float s = dot(p.n[i], c) - p.d[i];
        if (s <= 0.0f) { continue; }
        const Polyhedron::Face& f = poly.poly.faces[i];
        const Vec3 proj = c - p.n[i] * s;
        bool inside = true;
        for (int k = 0; k < f.count; ++k) {
            const Vec3 a = p.v[poly.poly.indices[f.first + k]];
            const Vec3 b = p.v[poly.poly.indices[f.first + (k + 1) % f.count]];
            if (dot(cross(b - a, p.n[i]), proj - a) > 0.0f) {
                inside = false;
                const Vec3 ab = b - a;
                const float t = std::fmin(std::fmax(dot(c - a, ab) / std::fmax(dot(ab, ab), FLT_MIN), 0.0f), 1.0f);
                const Vec3 q = a + ab * t;
                const float dsq = lengthSquared(c - q);
                if (dsq < closestSq) { closestSq = dsq; closest = q; }
            }
        }
        if (inside && s * s < closestSq) {
            closestSq = s * s;
            closest = proj;
        }
    }
    const float dist = std::sqrt(closestSq);
    const float sep = dist - r;
    if (sep > margin || dist <= 0.0f) { return; }
    const Vec3 n = (c - closest) / dist;
    m.normal = n;
    m.points[0] = { (closest + c - n * r) * 0.5f, sep };
    m.count = 1;
}

void collideSpheres(const Vec3& ca, float ra, const Vec3& cb, float rb, float margin, Manifold& m) {
    const Vec3 d = cb - ca;
    const float dist = length(d);
    const float sep = dist - ra - rb;
    if (sep > margin) { return; }
    const Vec3 n = dist > 0.0f ? d / dist : Vec3{ 0, 1, 0 };
    m.normal = n;
    m.points[0] = { ca + n * (ra + 0.5f * sep), sep };
    m.count = 1;
}

// Normal points out of the half-space.
void collideHalfSpace(const Transform& xh, const Shape& other, const Transform& xo, float margin, Manifold& m) {
    const Vec3 n = xh.q.rotate({ 0, 0, 1 });
    const float offset = dot(n, xh.p);
    if (other.kind == ShapeKind::Sphere) {
        const float d = dot(n, xo.p) - offset;
        const float sep = d - other.radius;
        if (sep > margin) { return; }
        m.normal = n;
        m.points[0] = { xo.p - n * (0.5f * (d + other.radius)), sep };
        m.count = 1;
        return;
    }
    std::vector<ManifoldPoint> pts;
    for (const Vec3& v : other.poly.vertices) {
        const Vec3 w = xo.apply(v);
        const float s = dot(n, w) - offset;
        if (s <= margin) { pts.push_back({ w - n * (0.5f * s), s }); }
    }
    reduceManifold(pts, n, m);
}

}

void reduceManifold(std::vector<ManifoldPoint>& pts, const Vec3& normal, Manifold& m) {
    m.normal = normal;
    m.count = 0;
    if (pts.size() <= 4) {
        for (const ManifoldPoint& p : pts) { m.points[m.count++] = p; }
        return;
    }

    size_t i0 = 0;
    for (size_t i = 1; i < pts.size(); ++i) {
        if (pts[i].separation < pts[i0].separation) { i0 = i; }
    }
    const Vec3 p0 = pts[i0].position;

    size_t i1 = i0;
    float best = -1.0f;
    for (size_t i = 0; i < pts.size(); ++i) {
        const float d = lengthSquared(pts[i].position - p0);
        if (d > best) { best = d; i1 = i; }
    }
    const Vec3 p1 = pts[i1].position;

    size_t i2 = i0;
    float area = 0.0f;
    for (size_t i = 0; i < pts.size(); ++i) {
        const float a = dot(cross(p1 - p0, pts[i].position - p0), normal);
        if (std::fabs(a) > std::fabs(area)) { area = a; i2 = i; }
    }
    const Vec3 p2 = pts[i2].position;
    const float sign = area < 0.0f ? -1.0f : 1.0f;

    m.points[m.count++] = pts[i0];
    if (i1 != i0) { m.points[m.count++] = pts[i1]; }
    if (i2 != i0 && i2 != i1) { m.points[m.count++] = pts[i2]; }
    if (m.count < 3) { return; }

    size_t i3 = pts.size();
    float outside = 0.0f;
    const Vec3 tri[3] = { p0, p1, p2 };
    for (size_t i = 0; i < pts.size(); ++i) {
        if (i == i0 || i == i1 || i == i2) { continue; }
        for (int e = 0; e < 3; ++e) {
            const Vec3& a = tri[e];
            const Vec3& b = tri[(e + 1) % 3];
            const float s = -sign * dot(cross(b - a, pts[i].position - a), normal);
            if (s > outside) { outside = s; i3 = i; }
        }
    }
    if (i3 < pts.size()) { m.points[m.count++] = pts[i3]; }
}

void collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
             float margin, const Tolerances& tol, Manifold& m) {
    m.count = 0;
    const ShapeKind ka = a.kind, kb = b.kind;
    if (ka == ShapeKind::HalfSpace && kb == ShapeKind::HalfSpace) { return; }

    if (ka == ShapeKind::HalfSpace) {
        collideHalfSpace(xa, b, xb, margin, m);
    } else if (kb == ShapeKind::HalfSpace) {
        collideHalfSpace(xb, a, xa, margin, m);
        m.normal = -m.normal;
    } else if (ka == ShapeKind::Sphere && kb == ShapeKind::Sphere) {
        collideSpheres(xa.p, a.radius, xb.p, b.radius, margin, m);
    } else if (ka == ShapeKind::Sphere) {
        collideSpherePoly(b, xb, xa.p, a.radius, margin, m);
        m.normal = -m.normal;
    } else if (kb == ShapeKind::Sphere) {
        collideSpherePoly(a, xa, xb.p, b.radius, margin, m);
    } else {
        collideConvex(a, xa, b, xb, margin, tol, m);
    }
}

}
