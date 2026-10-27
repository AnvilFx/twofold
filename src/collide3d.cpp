#include <cfloat>
#include <utility>

#include "internal.h"

namespace twofold {

namespace {

// Edge pairs closer to parallel than this sine give an unstable axis and are left to the face queries.
constexpr float kParallelEdges = 0.05f;

// A polyhedron seen from the frame of the pair's first shape.
struct FramePoly {
    const Polyhedron* src = nullptr;
    const Vec3* v = nullptr;
    const Vec3* n = nullptr;
    const float* d = nullptr;
    const Vec3* arc = nullptr;
    int vertexCount = 0;
    int faceCount = 0;
    Vec3 centroid, lo, hi;
};

struct FrameStore {
    std::vector<Vec3> v, n, arc;
    std::vector<float> d;
};

FramePoly ownFrame(const Polyhedron& p) {
    FramePoly f;
    f.src = &p;
    f.v = p.vertices.data();
    f.n = p.normals.data();
    f.d = p.offsets.data();
    f.arc = p.arcs.data();
    f.vertexCount = (int)p.vertices.size();
    f.faceCount = (int)p.faces.size();
    f.centroid = p.centroid;
    f.lo = p.lo;
    f.hi = p.hi;
    return f;
}

FramePoly otherFrame(const Polyhedron& p, const Transform& xf, FrameStore& store) {
    const Mat3 r = xf.q.toMat3();
    store.v.resize(p.vertices.size());
    store.n.resize(p.faces.size());
    store.d.resize(p.faces.size());
    store.arc.resize(p.edges.size());
    FramePoly f;
    f.src = &p;
    f.vertexCount = (int)p.vertices.size();
    f.faceCount = (int)p.faces.size();
    f.lo = f.hi = r * p.vertices[0] + xf.p;
    for (size_t i = 0; i < p.vertices.size(); ++i) {
        const Vec3 w = r * p.vertices[i] + xf.p;
        store.v[i] = w;
        f.lo = minVec(f.lo, w);
        f.hi = maxVec(f.hi, w);
    }
    for (size_t i = 0; i < p.faces.size(); ++i) {
        store.n[i] = r * p.normals[i];
        store.d[i] = dot(store.n[i], store.v[p.indices[p.faces[i].first]]);
    }
    for (size_t i = 0; i < p.edges.size(); ++i) { store.arc[i] = r * p.arcs[i]; }
    f.v = store.v.data();
    f.n = store.n.data();
    f.d = store.d.data();
    f.arc = store.arc.data();
    f.centroid = r * p.centroid + xf.p;
    return f;
}

int supportFace(const FramePoly& p, const Vec3& n) {
    int best = 0;
    for (int i = 1; i < p.faceCount; ++i) {
        if (dot(p.n[i], n) > dot(p.n[best], n)) { best = i; }
    }
    return best;
}

struct FaceQuery {
    int index = -1;
    float separation = -FLT_MAX;
};

struct EdgeQuery {
    int a = -1, b = -1;
    float separation = -FLT_MAX;
    Vec3 axis;
    bool pruned = false;
};

// Starts with the face of a toward b, and tests each face against the vertex of b deepest toward a
// first, so most faces are ruled out by one vertex. Stops once a face separates by more than limit.
FaceQuery queryFaces(const FramePoly& a, const FramePoly& b, float limit) {
    const Vec3 toward = b.centroid - a.centroid;
    int deep = 0;
    for (int k = 1; k < b.vertexCount; ++k) {
        if (dot(b.v[k], toward) < dot(b.v[deep], toward)) { deep = k; }
    }
    const int first = supportFace(a, toward);

    FaceQuery q;
    for (int j = -1; j < a.faceCount; ++j) {
        if (j == first) { continue; }
        const int i = j < 0 ? first : j;
        const Vec3 n = a.n[i];
        const float bound = q.separation + a.d[i];
        float s = dot(n, b.v[deep]);
        if (s > bound) {
            for (int k = 0; k < b.vertexCount; ++k) {
                s = std::fmin(s, dot(n, b.v[k]));
                if (s <= bound) { break; }
            }
        }
        s -= a.d[i];
        if (s > q.separation) {
            q.separation = s;
            q.index = i;
            if (s > limit) { break; }
        }
    }
    return q;
}

void edgesInBox(const FramePoly& p, const Vec3& lo, const Vec3& hi, std::vector<int>& out) {
    out.clear();
    const std::vector<Polyhedron::Edge>& edges = p.src->edges;
    for (size_t i = 0; i < edges.size(); ++i) {
        const Vec3 s = p.v[edges[i].v0], t = p.v[edges[i].v1];
        if (std::fmin(s.x, t.x) > hi.x || std::fmax(s.x, t.x) < lo.x) { continue; }
        if (std::fmin(s.y, t.y) > hi.y || std::fmax(s.y, t.y) < lo.y) { continue; }
        if (std::fmin(s.z, t.z) > hi.z || std::fmax(s.z, t.z) < lo.z) { continue; }
        out.push_back((int)i);
    }
}

// An edge pair is a candidate when its arcs cross on the Gauss map of the Minkowski difference.
// The deciding pair meets within reach of both shapes, so only edges reaching into the overlap of
// their bounds, grown by reach, are paired.
EdgeQuery queryEdges(const FramePoly& a, const FramePoly& b, float reach, std::vector<int>& nearA,
                     std::vector<int>& nearB) {
    EdgeQuery q;
    const Vec3 grow{ reach, reach, reach };
    const Vec3 lo = maxVec(a.lo, b.lo) - grow, hi = minVec(a.hi, b.hi) + grow;
    edgesInBox(a, lo, hi, nearA);
    edgesInBox(b, lo, hi, nearB);

    const std::vector<Polyhedron::Edge>& ea = a.src->edges;
    const std::vector<Polyhedron::Edge>& eb = b.src->edges;
    q.pruned = nearA.size() < ea.size() || nearB.size() < eb.size();
    for (int i : nearA) {
        const Vec3 ua = a.n[ea[i].f0], va = a.n[ea[i].f1], ga = a.arc[i];
        for (int j : nearB) {
            const float cba = dot(b.n[eb[j].f0], ga), dba = dot(b.n[eb[j].f1], ga);
            if (cba * dba >= 0.0f) { continue; }
            const Vec3 gb = b.arc[j];
            const float adc = dot(ua, gb), bdc = dot(va, gb);
            if (adc * bdc >= 0.0f || cba * bdc >= 0.0f) { continue; }

            const Vec3 pa = a.v[ea[i].v0];
            const Vec3 pb = b.v[eb[j].v0];
            Vec3 axis = cross(a.v[ea[i].v1] - pa, b.v[eb[j].v1] - pb);
            const float len = length(axis);
            if (len < kParallelEdges * ea[i].length * eb[j].length) { continue; }
            axis = axis / len;
            if (dot(axis, pa - a.centroid) < 0.0f) { axis = -axis; }
            const float s = dot(axis, pb - pa);
            if (s > q.separation) {
                q.separation = s;
                q.a = i;
                q.b = j;
                q.axis = axis;
            }
        }
    }
    return q;
}

void clipPolygon(std::vector<Vec3>& poly, std::vector<Vec3>& scratch, const Vec3& n, const Vec3& origin) {
    if (poly.empty()) { return; }
    scratch.clear();
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vec3& p = poly[i];
        const Vec3& q = poly[(i + 1) % poly.size()];
        const float dp = dot(n, p - origin);
        const float dq = dot(n, q - origin);
        if (dp <= 0.0f) { scratch.push_back(p); }
        if ((dp <= 0.0f) != (dq <= 0.0f)) { scratch.push_back(p + (q - p) * (dp / (dp - dq))); }
    }
    poly.swap(scratch);
}

struct ClipBuffers {
    std::vector<Vec3> poly;
    std::vector<Vec3> scratch;
    std::vector<ManifoldPoint> points;
};

// Points are found in the frame of the pair, then placed in the world by xf.
void faceContact(const FramePoly& ref, int refFace, const FramePoly& inc, bool flip, const Transform& xf,
                 float margin, float band, const Motion& motion, ClipBuffers& buf, Manifold& m) {
    const Vec3 refN = ref.n[refFace];
    const float refD = ref.d[refFace];
    const Vec3 worldN = xf.q.rotate(refN);

    int incFace = 0;
    float minDot = FLT_MAX;
    for (int i = 0; i < inc.faceCount; ++i) {
        const float d = dot(inc.n[i], refN);
        if (d < minDot) { minDot = d; incFace = i; }
    }

    const Polyhedron::Face& fi = inc.src->faces[incFace];
    buf.poly.clear();
    for (int k = 0; k < fi.count; ++k) { buf.poly.push_back(inc.v[inc.src->indices[fi.first + k]]); }

    const Polyhedron::Face& fr = ref.src->faces[refFace];
    for (int k = 0; k < fr.count; ++k) {
        const Vec3 a = ref.v[ref.src->indices[fr.first + k]];
        const Vec3 b = ref.v[ref.src->indices[fr.first + (k + 1) % fr.count]];
        clipPolygon(buf.poly, buf.scratch, cross(b - a, refN), a);
    }

    buf.points.clear();
    for (const Vec3& p : buf.poly) {
        const float s = dot(refN, p) - refD;
        if (s > margin) { continue; }
        const Vec3 w = xf.apply(p);
        buf.points.push_back({ w - worldN * (0.5f * s), s, motion.predict(w, worldN, s) });
    }
    // Small faces of a deep pair can miss each other along the normal, leaving the deepest vertex.
    if (buf.points.empty()) {
        Vec3 deepest = inc.v[0];
        for (int k = 1; k < inc.vertexCount; ++k) {
            if (dot(refN, inc.v[k]) < dot(refN, deepest)) { deepest = inc.v[k]; }
        }
        const float s = dot(refN, deepest) - refD;
        if (s <= margin) {
            const Vec3 w = xf.apply(deepest);
            buf.points.push_back({ w - worldN * (0.5f * s), s, motion.predict(w, worldN, s) });
        }
    }
    reduceManifold(buf.points, flip ? -worldN : worldN, band, m);
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

// Clips the face of b toward a against the face of a toward b, seen along the edge axis, so an edge
// lying across a face or another edge rests on more than one point.
void edgeContact(const FramePoly& a, const FramePoly& b, const EdgeQuery& eq, const Transform& xf, float margin,
                 float band, const Motion& motion, ClipBuffers& buf, Manifold& m) {
    const Vec3 n = eq.axis;
    const Vec3 worldN = xf.q.rotate(n);
    const int fa = supportFace(a, n);
    const float cosA = dot(a.n[fa], n);
    buf.points.clear();
    if (cosA > 0.0f) {
        const Polyhedron::Face& fi = b.src->faces[supportFace(b, -n)];
        buf.poly.clear();
        for (int k = 0; k < fi.count; ++k) { buf.poly.push_back(b.v[b.src->indices[fi.first + k]]); }
        const Polyhedron::Face& fr = a.src->faces[fa];
        for (int k = 0; k < fr.count; ++k) {
            const Vec3 p = a.v[a.src->indices[fr.first + k]];
            const Vec3 q = a.v[a.src->indices[fr.first + (k + 1) % fr.count]];
            clipPolygon(buf.poly, buf.scratch, cross(q - p, n), p);
        }
        // No point needs a deeper push than the separating axis allows.
        for (const Vec3& p : buf.poly) {
            const float s = std::fmax((dot(a.n[fa], p) - a.d[fa]) / cosA, eq.separation);
            if (s > margin) { continue; }
            const Vec3 w = xf.apply(p);
            buf.points.push_back({ w - worldN * (0.5f * s), s, motion.predict(w, worldN, s) });
        }
    }
    if (buf.points.size() >= 2) {
        reduceManifold(buf.points, worldN, band, m);
        return;
    }

    const Polyhedron::Edge& ea = a.src->edges[eq.a];
    const Polyhedron::Edge& eb = b.src->edges[eq.b];
    Vec3 c1, c2;
    closestOnSegments(a.v[ea.v0], a.v[ea.v1], b.v[eb.v0], b.v[eb.v1], c1, c2);
    m.normal = worldN;
    m.points[0] = { xf.apply((c1 + c2) * 0.5f), eq.separation };
    m.count = 1;
}

struct ConvexBuffers {
    FrameStore frame;
    std::vector<int> nearA, nearB;
    ClipBuffers clip;
};

void collideConvex(const Shape& sa, const Transform& xa, const Shape& sb, const Transform& xb,
                   float margin, const Tolerances& tol, const Motion& motion, Manifold& m) {
    thread_local ConvexBuffers buf;
    const FramePoly a = ownFrame(sa.poly);
    const FramePoly b = otherFrame(sb.poly, { xa.q.unrotate(xb.p - xa.p), xa.q.conjugate() * xb.q }, buf.frame);

    const FaceQuery fa = queryFaces(a, b, margin);
    if (fa.separation > margin) { return; }
    const FaceQuery fb = queryFaces(b, a, margin);
    if (fb.separation > margin) { return; }
    // An edge axis only wins with less overlap than the best face axis.
    const float reach = std::fmax(margin, -std::fmax(fa.separation, fb.separation));
    const EdgeQuery eq = queryEdges(a, b, reach, buf.nearA, buf.nearB);
    if (eq.separation > margin) { return; }
    // Pruned edges can hold the only axis separating shapes that are further apart than margin.
    if (eq.pruned && distance(convexOf(sa), xa, convexOf(sb), xb).distance > margin) { return; }

    const float absTol = 0.5f * tol.linearSlop;
    if (eq.a >= 0 && eq.separation > 0.9f * std::fmax(fa.separation, fb.separation) + absTol) {
        edgeContact(a, b, eq, xa, margin, tol.linearSlop, motion, buf.clip, m);
    } else if (fb.separation > 0.98f * fa.separation + absTol) {
        faceContact(b, fb.index, a, true, xa, margin, tol.linearSlop, motion.swapped(), buf.clip, m);
    } else {
        faceContact(a, fa.index, b, false, xa, margin, tol.linearSlop, motion, buf.clip, m);
    }
}

// Normal points from the polyhedron to the sphere.
void collideSpherePoly(const Shape& poly, const Transform& xp, const Vec3& center, float r, float margin, Manifold& m) {
    const Polyhedron& p = poly.poly;
    const Vec3 c = xp.applyInverse(center);
    int best = -1;
    float bestSep = -FLT_MAX;
    bool outside = false;
    for (size_t i = 0; i < p.faces.size(); ++i) {
        const float s = dot(p.normals[i], c) - p.offsets[i];
        if (s > r + margin) { return; }
        if (s > 0.0f) { outside = true; }
        if (s > bestSep) { bestSep = s; best = (int)i; }
    }
    if (best < 0) { return; }

    if (!outside) {
        const Vec3 n = p.normals[best];
        m.normal = xp.q.rotate(n);
        m.points[0] = { xp.apply(c - n * (0.5f * (bestSep + r))), bestSep - r };
        m.count = 1;
        return;
    }

    Vec3 closest;
    float closestSq = FLT_MAX;
    for (size_t i = 0; i < p.faces.size(); ++i) {
        const Polyhedron::Face& f = p.faces[i];
        const Vec3 fn = p.normals[i];
        const float s = dot(fn, c) - p.offsets[i];
        if (s <= 0.0f) { continue; }
        const Vec3 proj = c - fn * s;
        bool inside = true;
        for (int k = 0; k < f.count; ++k) {
            const Vec3 a = p.vertices[p.indices[f.first + k]];
            const Vec3 b = p.vertices[p.indices[f.first + (k + 1) % f.count]];
            if (dot(cross(b - a, fn), proj - a) > 0.0f) {
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
    m.normal = xp.q.rotate(n);
    m.points[0] = { xp.apply((closest + c - n * r) * 0.5f), sep };
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
void collideHalfSpace(const Transform& xh, const Shape& other, const Transform& xo, float margin, float band,
                      const Motion& motion, Manifold& m) {
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
    thread_local std::vector<ManifoldPoint> pts;
    pts.clear();
    const Vec3 localN = xo.q.unrotate(n);
    const float localOffset = offset - dot(n, xo.p);
    for (const Vec3& v : other.poly.vertices) {
        const float s = dot(localN, v) - localOffset;
        if (s > margin) { continue; }
        const Vec3 w = xo.apply(v);
        pts.push_back({ w - n * (0.5f * s), s, motion.predict(w, n, s) });
    }
    reduceManifold(pts, n, band, m);
}

}

// Ranks points by predicted separation. Spread is chosen among points within band of the deepest,
// then free slots take the deepest of the rest.
void reduceManifold(std::vector<ManifoldPoint>& pts, const Vec3& normal, float band, Manifold& m) {
    m.normal = normal;
    m.count = 0;
    if (pts.size() <= 4) {
        for (const ManifoldPoint& p : pts) { m.points[m.count++] = p; }
        return;
    }

    size_t i0 = 0;
    for (size_t i = 1; i < pts.size(); ++i) {
        if (pts[i].predicted < pts[i0].predicted) { i0 = i; }
    }
    const Vec3 p0 = pts[i0].position;
    const float limit = pts[i0].predicted + band;
    size_t chosen[4] = { i0, i0, i0, i0 };
    auto used = [&](size_t i) { return i == chosen[0] || i == chosen[1] || i == chosen[2] || i == chosen[3]; };
    m.points[m.count++] = pts[i0];

    size_t i1 = pts.size();
    float best = 0.0f;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (pts[i].predicted > limit) { continue; }
        const float d = lengthSquared(pts[i].position - p0);
        if (d > best) { best = d; i1 = i; }
    }

    size_t i2 = pts.size();
    float area = 0.0f;
    if (i1 < pts.size()) {
        chosen[1] = i1;
        m.points[m.count++] = pts[i1];
        const Vec3 p1 = pts[i1].position;
        for (size_t i = 0; i < pts.size(); ++i) {
            if (used(i) || pts[i].predicted > limit) { continue; }
            const float a = dot(cross(p1 - p0, pts[i].position - p0), normal);
            if (std::fabs(a) > std::fabs(area)) { area = a; i2 = i; }
        }
    }

    if (i2 < pts.size()) {
        chosen[2] = i2;
        m.points[m.count++] = pts[i2];
        const float sign = area < 0.0f ? -1.0f : 1.0f;
        const Vec3 tri[3] = { p0, m.points[1].position, pts[i2].position };
        size_t i3 = pts.size();
        float outside = 0.0f;
        for (size_t i = 0; i < pts.size(); ++i) {
            if (used(i) || pts[i].predicted > limit) { continue; }
            for (int e = 0; e < 3; ++e) {
                const Vec3& a = tri[e];
                const Vec3& b = tri[(e + 1) % 3];
                const float s = -sign * dot(cross(b - a, pts[i].position - a), normal);
                if (s > outside) { outside = s; i3 = i; }
            }
        }
        if (i3 < pts.size()) {
            chosen[3] = i3;
            m.points[m.count++] = pts[i3];
        }
    }

    while (m.count < 4) {
        size_t next = pts.size();
        for (size_t i = 0; i < pts.size(); ++i) {
            if (used(i)) { continue; }
            if (next == pts.size() || pts[i].predicted < pts[next].predicted) { next = i; }
        }
        if (next == pts.size()) { break; }
        chosen[m.count] = next;
        m.points[m.count++] = pts[next];
    }
}

void collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
             float margin, const Tolerances& tol, const Motion& motion, Manifold& m) {
    m.count = 0;
    const ShapeKind ka = a.kind, kb = b.kind;
    if (ka == ShapeKind::HalfSpace && kb == ShapeKind::HalfSpace) { return; }

    if (ka == ShapeKind::HalfSpace) {
        collideHalfSpace(xa, b, xb, margin, tol.linearSlop, motion, m);
    } else if (kb == ShapeKind::HalfSpace) {
        collideHalfSpace(xb, a, xa, margin, tol.linearSlop, motion.swapped(), m);
        m.normal = -m.normal;
    } else if (ka == ShapeKind::Sphere && kb == ShapeKind::Sphere) {
        collideSpheres(xa.p, a.radius, xb.p, b.radius, margin, m);
    } else if (ka == ShapeKind::Sphere) {
        collideSpherePoly(b, xb, xa.p, a.radius, margin, m);
        m.normal = -m.normal;
    } else if (kb == ShapeKind::Sphere) {
        collideSpherePoly(a, xa, xb.p, b.radius, margin, m);
    } else {
        collideConvex(a, xa, b, xb, margin, tol, motion, m);
    }
}

}
