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

// Rim points standing in for a cap.
constexpr int kCapPoints = 24;
// A polyhedron face takes part as a face when its normal is this close to the contact normal.
constexpr float kFlatFace = 0.9998f;

// A cylinder or a polyhedron in the frame of the pair.
struct Solid {
    bool cylinder = false;
    Vec3 center, axis, side;
    float radius = 0.0f, halfHeight = 0.0f;
    FramePoly poly;
};

Solid cylinderSolid(const Shape& s, const Transform& xf) {
    Solid c;
    c.cylinder = true;
    c.center = xf.p;
    c.axis = xf.q.rotate({ 0, 0, 1 });
    c.side = xf.q.rotate({ 1, 0, 0 });
    c.radius = s.radius;
    c.halfHeight = s.halfHeight;
    return c;
}

Solid solidOf(const Shape& s, const Transform& xf, FrameStore& store) {
    if (s.kind == ShapeKind::Cylinder) { return cylinderSolid(s, xf); }
    Solid p;
    p.poly = otherFrame(s.poly, xf, store);
    return p;
}

float extentAlong(const Solid& s, const Vec3& d) {
    if (s.cylinder) {
        const float a = dot(s.axis, d);
        return dot(s.center, d) + s.halfHeight * std::fabs(a) + s.radius * std::sqrt(std::fmax(1.0f - a * a, 0.0f));
    }
    float most = dot(s.poly.v[0], d);
    for (int k = 1; k < s.poly.vertexCount; ++k) { most = std::fmax(most, dot(s.poly.v[k], d)); }
    return most;
}

Vec3 supportPoint(const Solid& s, const Vec3& d) {
    if (s.cylinder) {
        const float a = dot(s.axis, d);
        const Vec3 radial = d - s.axis * a;
        const float len = length(radial);
        const Vec3 rim = len > 1e-6f ? radial * (s.radius / len) : Vec3{};
        return s.center + s.axis * (a < 0.0f ? -s.halfHeight : s.halfHeight) + rim;
    }
    int best = 0;
    for (int k = 1; k < s.poly.vertexCount; ++k) {
        if (dot(s.poly.v[k], d) > dot(s.poly.v[best], d)) { best = k; }
    }
    return s.poly.v[best];
}

// The part of a solid farthest along d: a face with its outward normal, an edge or a point. A cylinder
// gives a cap within 45 degrees of its axis, and otherwise the line of its side.
void featureAlong(const Solid& s, const Vec3& d, float band, std::vector<Vec3>& out, Vec3& normal) {
    out.clear();
    if (s.cylinder) {
        const float a = dot(s.axis, d);
        const Vec3 radial = d - s.axis * a;
        const float len = length(radial);
        if (std::fabs(a) > len) {
            normal = a < 0.0f ? -s.axis : s.axis;
            // A tilted cap starts at its deepest point, so that point is exact. A level one keeps the same
            // points as it turns, so its contacts carry over between steps.
            const Vec3 u = 2.0f * s.radius * len > band ? radial / len : s.side;
            const Vec3 v = cross(normal, u);
            const Vec3 c = s.center + normal * s.halfHeight;
            for (int k = 0; k < kCapPoints; ++k) {
                const float t = 2.0f * 3.14159265358979f * (float)k / (float)kCapPoints;
                out.push_back(c + (u * std::cos(t) + v * std::sin(t)) * s.radius);
            }
        } else {
            const Vec3 p = s.center + radial * (s.radius / len);
            out.push_back(p - s.axis * s.halfHeight);
            out.push_back(p + s.axis * s.halfHeight);
        }
        return;
    }
    const FramePoly& p = s.poly;
    const int f = supportFace(p, d);
    if (dot(p.n[f], d) >= kFlatFace) {
        normal = p.n[f];
        const Polyhedron::Face& face = p.src->faces[f];
        for (int k = 0; k < face.count; ++k) { out.push_back(p.v[p.src->indices[face.first + k]]); }
        return;
    }
    int first = 0, second = -1;
    for (int k = 1; k < p.vertexCount; ++k) {
        if (dot(p.v[k], d) > dot(p.v[first], d)) {
            second = first;
            first = k;
        } else if (second < 0 || dot(p.v[k], d) > dot(p.v[second], d)) {
            second = k;
        }
    }
    out.push_back(p.v[first]);
    if (second >= 0 && dot(p.v[second], d) >= dot(p.v[first], d) - band) { out.push_back(p.v[second]); }
}

struct Separation {
    Vec3 normal;
    float gap = -FLT_MAX;
};

void tryAxis(const Solid& a, const Solid& b, const Vec3& d, Separation& best) {
    const float len = length(d);
    if (len < 1e-6f) { return; }
    const Vec3 n = d / len;
    const float gap = -extentAlong(a, n) - extentAlong(b, -n);
    if (gap > best.gap) {
        best.normal = n;
        best.gap = gap;
    }
}

// No direction separates the solids by more than they are apart, so the best direction along which a flat
// or straight pair of features can touch is exact whenever such a pair touches.
Separation featureAxes(const Solid& a, const Solid& b) {
    Separation best;
    for (int i = 0; !a.cylinder && i < a.poly.faceCount; ++i) { tryAxis(a, b, a.poly.n[i], best); }
    for (int i = 0; !b.cylinder && i < b.poly.faceCount; ++i) { tryAxis(a, b, -b.poly.n[i], best); }
    for (const Solid* c : { &a, &b }) {
        if (!c->cylinder) { continue; }
        tryAxis(a, b, c->axis, best);
        tryAxis(a, b, -c->axis, best);
    }
    if (a.cylinder && b.cylinder) {
        const Vec3 c = cross(a.axis, b.axis);
        if (length(c) > kParallelEdges) {
            tryAxis(a, b, c, best);
            tryAxis(a, b, -c, best);
        } else {
            const Vec3 d = b.center - a.center;
            tryAxis(a, b, d - a.axis * dot(d, a.axis), best);
        }
        return best;
    }
    const Solid& cyl = a.cylinder ? a : b;
    const FramePoly& p = a.cylinder ? b.poly : a.poly;
    for (const Polyhedron::Edge& e : p.src->edges) {
        const Vec3 edge = p.v[e.v1] - p.v[e.v0];
        Vec3 c = cross(edge, cyl.axis);
        if (length(c) <= kParallelEdges * e.length) {
            const Vec3 d = p.v[e.v0] - cyl.center;
            c = d - cyl.axis * dot(d, cyl.axis);
        }
        tryAxis(a, b, c, best);
        tryAxis(a, b, -c, best);
    }
    return best;
}

// GJK directions between the shapes and between their cores, cylinders shrunk into a rounded core. Near
// touching, the first is poorly placed and the second takes over.
int distanceAxes(Convex ca, const Transform& xa, Convex cb, const Transform& xb, const Tolerances& tol, Vec3 n[2]) {
    int count = 0;
    const DistanceResult r = distance(ca, xa, cb, xb);
    if (r.distance > 0.0f) { n[count++] = r.normal; }
    if (r.distance > tol.linearSlop) { return count; }
    float shrunk = 0.0f;
    for (Convex* c : { &ca, &cb }) {
        if (!c->cylinder) { continue; }
        const float core = 0.25f * std::fmin(c->rim, c->halfHeight);
        c->rim -= core;
        c->halfHeight -= core;
        c->radius += core;
        shrunk += core;
    }
    const DistanceResult core = distance(ca, xa, cb, xb);
    if (core.distance + shrunk > 0.0f) { n[count++] = core.normal; }
    return count;
}

Vec3 centerOf(const Solid& s) {
    if (s.cylinder) { return s.center; }
    Vec3 c;
    for (int k = 0; k < s.poly.vertexCount; ++k) { c = c + s.poly.v[k]; }
    return c / (float)s.poly.vertexCount;
}

// Tilts from start while that separates the solids more. No direction separates them by more than they are
// apart, so each step only comes closer to the true one.
Separation climb(const Solid& a, const Solid& b, Separation start) {
    constexpr float kDiagonal = 0.70710678f;
    const float turns[8][2] = { { 1, 0 }, { kDiagonal, kDiagonal }, { 0, 1 }, { -kDiagonal, kDiagonal },
                                { -1, 0 }, { -kDiagonal, -kDiagonal }, { 0, -1 }, { kDiagonal, -kDiagonal } };
    Separation best = start;
    float step = 0.1f;
    for (int i = 0; i < 64 && step > 1e-3f; ++i) {
        const Vec3 u = perpendicular(best.normal), v = cross(best.normal, u);
        Separation next = best;
        for (const float* t : turns) { tryAxis(a, b, best.normal + (u * t[0] + v * t[1]) * step, next); }
        if (next.gap > best.gap) {
            best = next;
        } else {
            step *= 0.5f;
        }
    }
    return best;
}

// False when the solids are further apart than margin. A GJK direction only wins where no flat or straight
// feature pair comes close, such as a rim on an edge. In a deep overlap the GJK directions are approximate or
// missing, so the best axis climbs toward the direction a curved pair touches along.
bool separate(const Solid& a, const Solid& b, const Convex& ca, const Transform& xa, const Convex& cb,
              const Transform& xb, const Quat& toFrame, float margin, const Tolerances& tol, Separation& best) {
    best = featureAxes(a, b);
    if (best.gap > margin) { return false; }
    Vec3 n[2];
    const int count = distanceAxes(ca, xa, cb, xb, tol, n);
    Separation curved;
    for (int i = 0; i < count; ++i) { tryAxis(a, b, toFrame.unrotate(n[i]), curved); }
    if (curved.gap > best.gap + 0.1f * tol.linearSlop) { best = curved; }
    if (best.gap < -2.0f * tol.speculative) {
        // The gap has several peaks. A rim pressed into the other solid peaks near its cap normal or its side.
        Separation climbed = climb(a, b, best);
        const Vec3 d = centerOf(b) - centerOf(a);
        for (const Solid* s : { &a, &b }) {
            if (!s->cylinder) { continue; }
            Separation cap, side;
            tryAxis(a, b, s->axis, cap);
            tryAxis(a, b, -s->axis, cap);
            tryAxis(a, b, d - s->axis * dot(d, s->axis), side);
            for (const Separation& start : { cap, side }) {
                if (start.gap == -FLT_MAX) { continue; }
                const Separation peak = climb(a, b, start);
                if (peak.gap > climbed.gap) { climbed = peak; }
            }
        }
        if (climbed.gap > best.gap + 0.1f * tol.linearSlop) { best = climbed; }
    }
    return best.gap <= margin;
}

void clipSegment(std::vector<Vec3>& seg, const Vec3& n, const Vec3& origin) {
    const float dp = dot(n, seg[0] - origin), dq = dot(n, seg[1] - origin);
    if (dp > 0.0f && dq > 0.0f) {
        seg.clear();
    } else if (dp > 0.0f) {
        seg[0] = seg[0] + (seg[1] - seg[0]) * (dp / (dp - dq));
    } else if (dq > 0.0f) {
        seg[1] = seg[1] + (seg[0] - seg[1]) * (dq / (dq - dp));
    }
}

struct RoundBuffers {
    FrameStore frame;
    std::vector<Vec3> fa, fb, turn;
    ClipBuffers clip;
};

// The side lines one cap point either side of the one along d. A spinning cylinder brings them into
// contact within a step.
void turnLines(const Solid& s, const Vec3& d, std::vector<Vec3>& out) {
    out.clear();
    const Vec3 radial = normalize(d - s.axis * dot(s.axis, d));
    const Vec3 across = cross(s.axis, radial);
    const float t = 2.0f * 3.14159265358979f / (float)kCapPoints;
    for (float sign : { -1.0f, 1.0f }) {
        const Vec3 p = s.center + (radial * std::cos(t) + across * (sign * std::sin(t))) * s.radius;
        out.push_back(p - s.axis * s.halfHeight);
        out.push_back(p + s.axis * s.halfHeight);
    }
}

// Contacts of solids that separate by gap along n, from a to b. Whichever feature is a face clips the
// other, and a face pair is clipped by the face that best faces along n.
void solidContact(const Solid& a, const Solid& b, const Vec3& n, float gap, const Transform& xf, float margin,
                  float band, const Motion& motion, RoundBuffers& buf, Manifold& m) {
    Vec3 na, nb;
    featureAlong(a, n, band, buf.fa, na);
    featureAlong(b, -n, band, buf.fb, nb);
    const Vec3 worldN = xf.q.rotate(n);
    std::vector<ManifoldPoint>& pts = buf.clip.points;
    pts.clear();
    // No point needs a deeper push than the separating direction allows.
    auto add = [&](const Vec3& onA, const Vec3& onB, float s) {
        s = std::fmax(s, gap);
        if (s > margin) { return; }
        const Vec3 w = xf.apply((onA + onB) * 0.5f);
        pts.push_back({ w, s, motion.predict(w, worldN, s) });
    };

    const bool faceA = buf.fa.size() >= 3, faceB = buf.fb.size() >= 3;
    const bool refA = faceA && (!faceB || dot(na, n) >= -dot(nb, n));
    auto clipToFace = [&](const Vec3* feature, size_t count) {
        const std::vector<Vec3>& ref = refA ? buf.fa : buf.fb;
        const Vec3 refN = refA ? na : nb;
        const Vec3 out = refA ? n : -n;
        const float facing = dot(refN, out);
        std::vector<Vec3>& inc = buf.clip.poly;
        inc.assign(feature, feature + count);
        for (size_t k = 0; k < ref.size() && !inc.empty(); ++k) {
            const Vec3 side = cross(ref[(k + 1) % ref.size()] - ref[k], refN);
            if (inc.size() >= 3) {
                clipPolygon(inc, buf.clip.scratch, side, ref[k]);
            } else if (inc.size() == 2) {
                clipSegment(inc, side, ref[k]);
            } else if (dot(side, inc[0] - ref[k]) > 0.0f) {
                inc.clear();
            }
        }
        for (const Vec3& q : inc) {
            const float s = dot(refN, q - ref[0]) / facing;
            const Vec3 onRef = q - out * s;
            if (refA) {
                add(onRef, q, s);
            } else {
                add(q, onRef, s);
            }
        }
    };
    // Lines closer to parallel than kParallelEdges touch along their overlap, others at their closest points.
    auto linePair = [&](const Vec3& a0, const Vec3& a1, const Vec3& b0, const Vec3& b1) {
        const Vec3 da = a1 - a0, db = b1 - b0;
        const float la = dot(da, da), lb = dot(db, db);
        if (la > 0.0f && lb > 0.0f && length(cross(da, db)) <= kParallelEdges * std::sqrt(la * lb)) {
            const float t0 = dot(b0 - a0, da) / la, t1 = dot(b1 - a0, da) / la;
            const float lo = std::fmax(std::fmin(t0, t1), 0.0f), hi = std::fmin(std::fmax(t0, t1), 1.0f);
            if (lo <= hi) {
                for (float t : { lo, hi }) {
                    const Vec3 onA = a0 + da * t;
                    const Vec3 onB = b0 + db * std::fmin(std::fmax(dot(onA - b0, db) / lb, 0.0f), 1.0f);
                    add(onA, onB, dot(onB - onA, n));
                }
                return;
            }
        }
        // Lines tilted against n reach deeper at their ends than where they pass closest.
        Vec3 onA, onB;
        closestOnSegments(a0, a1, b0, b1, onA, onB);
        add(onA, onB, std::fmin(dot(b0, n), dot(b1, n)) - std::fmax(dot(a0, n), dot(a1, n)));
    };

    const Solid& incident = refA ? b : a;
    const std::vector<Vec3>& incFeature = refA ? buf.fb : buf.fa;
    if (faceA || faceB) {
        clipToFace(incFeature.data(), incFeature.size());
    } else if (buf.fa.size() == 2 && buf.fb.size() == 2) {
        linePair(buf.fa[0], buf.fa[1], buf.fb[0], buf.fb[1]);
    }
    if (pts.empty()) {
        const Vec3 onB = buf.fa.size() == 1 && buf.fb.size() != 1 ? buf.fa[0] + n * gap : supportPoint(b, -n);
        add(onB - n * gap, onB, gap);
    }
    if ((faceA || faceB) && incident.cylinder && incFeature.size() == 2) {
        turnLines(incident, refA ? -n : n, buf.turn);
        for (int k = 0; k < 2; ++k) { clipToFace(&buf.turn[2 * k], 2); }
    } else if (!faceA && !faceB && buf.fa.size() == 2 && buf.fb.size() == 2) {
        if (a.cylinder) {
            turnLines(a, n, buf.turn);
            for (int k = 0; k < 2; ++k) { linePair(buf.turn[2 * k], buf.turn[2 * k + 1], buf.fb[0], buf.fb[1]); }
        }
        if (b.cylinder) {
            turnLines(b, -n, buf.turn);
            for (int k = 0; k < 2; ++k) { linePair(buf.fa[0], buf.fa[1], buf.turn[2 * k], buf.turn[2 * k + 1]); }
        }
    }
    reduceManifold(pts, worldN, band, m, 0.5f * band);
}

// Collides two shapes in the world frame, at least one of them a cylinder.
void collideRound(const Shape& sa, const Transform& xa, const Shape& sb, const Transform& xb, float margin,
                  const Tolerances& tol, const Motion& motion, Manifold& m) {
    thread_local RoundBuffers buf;
    const Solid a = solidOf(sa, xa, buf.frame);
    const Solid b = solidOf(sb, xb, buf.frame);
    Separation best;
    if (!separate(a, b, convexOf(sa), xa, convexOf(sb), xb, Quat{}, margin, tol, best)) { return; }
    solidContact(a, b, best.normal, best.gap, Transform{}, margin, tol.linearSlop, motion, buf, m);
}

// Normal points from the cylinder to the sphere.
void collideSphereCylinder(const Shape& cyl, const Transform& xc, const Vec3& center, float r, float margin,
                           Manifold& m) {
    const Vec3 c = xc.applyInverse(center);
    const float rho = std::sqrt(c.x * c.x + c.y * c.y);
    Vec3 n, surface;
    float sep;
    if (rho <= cyl.radius && std::fabs(c.z) <= cyl.halfHeight) {
        const float side = cyl.radius - rho, cap = cyl.halfHeight - std::fabs(c.z);
        if (side < cap) {
            n = rho > 0.0f ? Vec3{ c.x / rho, c.y / rho, 0.0f } : Vec3{ 1, 0, 0 };
            surface = c + n * side;
            sep = -side - r;
        } else {
            n = { 0.0f, 0.0f, c.z < 0.0f ? -1.0f : 1.0f };
            surface = c + n * cap;
            sep = -cap - r;
        }
    } else {
        const float scale = rho > cyl.radius ? cyl.radius / rho : 1.0f;
        surface = { c.x * scale, c.y * scale, std::fmin(std::fmax(c.z, -cyl.halfHeight), cyl.halfHeight) };
        const float dist = length(c - surface);
        if (dist <= 0.0f) { return; }
        n = (c - surface) / dist;
        sep = dist - r;
    }
    if (sep > margin) { return; }
    m.normal = xc.q.rotate(n);
    m.points[0] = { xc.apply((surface + c - n * r) * 0.5f), sep };
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
    if (other.kind == ShapeKind::Cylinder) {
        thread_local std::vector<Vec3> feature;
        Vec3 unused;
        featureAlong(cylinderSolid(other, xo), -n, band, feature, unused);
        for (const Vec3& w : feature) {
            const float s = dot(n, w) - offset;
            if (s > margin) { continue; }
            pts.push_back({ w - n * (0.5f * s), s, motion.predict(w, n, s) });
        }
        reduceManifold(pts, n, band, m, 0.5f * band);
        return;
    }
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

// Triangle manifolds whose normals are this close merge into one.
constexpr float kSameNormal = 0.996f;

// A triangle as a polyhedron with a front and a back face, so the face queries and clipping above apply.
struct TriangleFrame {
    Polyhedron poly;

    TriangleFrame() {
        poly.vertices.resize(3);
        poly.indices = { 0, 1, 2, 0, 2, 1 };
        poly.faces = { { 0, 3 }, { 3, 3 } };
        poly.normals.resize(2);
        poly.offsets.resize(2);
        poly.edges.resize(3);
        for (int k = 0; k < 3; ++k) {
            poly.edges[k].v0 = k;
            poly.edges[k].v1 = (k + 1) % 3;
            poly.edges[k].f0 = 0;
            poly.edges[k].f1 = 1;
        }
    }

    void set(const Vec3 c[3], const Vec3& n) {
        for (int k = 0; k < 3; ++k) {
            poly.vertices[k] = c[k];
            poly.edges[k].length = length(c[(k + 1) % 3] - c[k]);
        }
        poly.normals[0] = n;
        poly.normals[1] = -n;
        poly.offsets[0] = dot(n, c[0]);
        poly.offsets[1] = -poly.offsets[0];
        poly.centroid = (c[0] + c[1] + c[2]) * (1.0f / 3.0f);
        poly.lo = minVec(minVec(c[0], c[1]), c[2]);
        poly.hi = maxVec(maxVec(c[0], c[1]), c[2]);
    }
};

// Edge queries against a triangle. A triangle edge's arc on the Gauss map is the half circle from the
// front normal through the edge's outward direction to the back normal.
EdgeQuery queryTriangleEdges(const FramePoly& a, const FramePoly& b, float reach, std::vector<int>& nearB) {
    EdgeQuery q;
    const Vec3 grow{ reach, reach, reach };
    edgesInBox(b, maxVec(a.lo, b.lo) - grow, minVec(a.hi, b.hi) + grow, nearB);
    const std::vector<Polyhedron::Edge>& eb = b.src->edges;
    q.pruned = nearB.size() < eb.size();
    const Vec3 n = a.n[0];
    for (int k = 0; k < 3; ++k) {
        const Vec3 pa = a.v[k];
        const Vec3 e = a.v[(k + 1) % 3] - pa;
        const Vec3 out = cross(e, n);
        const float lengthA = a.src->edges[k].length;
        for (int j : nearB) {
            const Vec3 u = b.n[eb[j].f0], v = b.n[eb[j].f1];
            const float du = dot(u, e), dv = dot(v, e);
            if (du * dv >= 0.0f) { continue; }
            // Where b's arc crosses the plane of a's arc. The pair forms a face when a's normal there is
            // the opposite direction.
            const Vec3 p = (u * dv - v * du) * (dv > 0.0f ? 1.0f : -1.0f);
            if (dot(p, out) >= 0.0f) { continue; }

            const Vec3 pb = b.v[eb[j].v0];
            Vec3 axis = cross(e, b.v[eb[j].v1] - pb);
            const float len = length(axis);
            if (len < kParallelEdges * lengthA * eb[j].length) { continue; }
            axis = axis / len;
            if (dot(axis, p) > 0.0f) { axis = -axis; }
            const float s = dot(axis, pb - pa);
            if (s > q.separation) {
                q.separation = s;
                q.a = k;
                q.b = j;
                q.axis = axis;
            }
        }
    }
    return q;
}

// The point of the triangle closest to p, and its feature: corner k is k, edge k is 3 + k, the face is 6.
Vec3 closestOnTriangle(const Vec3 c[3], const Vec3& p, int& feature) {
    const Vec3 ab = c[1] - c[0], ac = c[2] - c[0], ap = p - c[0];
    const float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f) {
        feature = 0;
        return c[0];
    }
    const Vec3 bp = p - c[1];
    const float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) {
        feature = 1;
        return c[1];
    }
    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        feature = 3;
        return c[0] + ab * (d1 / (d1 - d3));
    }
    const Vec3 cp = p - c[2];
    const float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) {
        feature = 2;
        return c[2];
    }
    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        feature = 5;
        return c[0] + ac * (d2 / (d2 - d6));
    }
    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && d4 - d3 >= 0.0f && d5 - d6 >= 0.0f) {
        feature = 4;
        return c[1] + (c[2] - c[1]) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    feature = 6;
    const float denom = 1.0f / (va + vb + vc);
    return c[0] + ab * (vb * denom) + ac * (vc * denom);
}

// The edge or corner of the triangle farthest along d, numbered as above. Only the face has the face's
// normal, so any other d meets an edge or corner: the edge when its far corner is within tolerance.
int triangleFeature(const Vec3 c[3], const Vec3& d, float tolerance) {
    const float h[3] = { dot(c[0], d), dot(c[1], d), dot(c[2], d) };
    const int top = h[1] > h[0] ? (h[2] > h[1] ? 2 : 1) : (h[2] > h[0] ? 2 : 0);
    const int next = (top + 1) % 3, prev = (top + 2) % 3;
    if (h[next] >= h[prev]) { return h[next] >= h[top] - tolerance ? 3 + top : top; }
    return h[prev] >= h[top] - tolerance ? 3 + prev : top;
}

// Whether p lies over the triangle along its normal n, within tolerance of its edges.
bool overTriangle(const Vec3 c[3], const Vec3& n, const Vec3& p, float tolerance) {
    for (int k = 0; k < 3; ++k) {
        const Vec3 e = c[(k + 1) % 3] - c[k];
        if (dot(cross(e, n), p - c[k]) > tolerance * length(e)) { return false; }
    }
    return true;
}

// Whether a contact on the feature may keep a normal other than the face's, seen from side 0 (front) or 1.
bool ownNormal(uint8_t edges, int feature, int side) {
    const int shift = 3 * side;
    if (feature >= 3 && feature < 6) { return (edges >> (shift + feature - 3)) & 1; }
    if (feature < 3) { return ((edges >> (shift + feature)) & 1) || ((edges >> (shift + (feature + 2) % 3)) & 1); }
    return true;
}

// A body beside the triangle can be on one side while the normal leans to the other. A feature with a
// normal of its own on either side keeps the normal.
bool keepsNormal(uint8_t edges, int feature, int side, const Vec3& normal, const Vec3& n) {
    return ownNormal(edges, feature, side) || ownNormal(edges, feature, dot(normal, n) >= 0.0f ? 0 : 1);
}

struct MeshBuffers {
    TriangleFrame tri;
    FrameStore frame;
    std::vector<int> near;
    ClipBuffers clip;
    RoundBuffers round;
    std::vector<Manifold> found;
    std::vector<int> owner;
    std::vector<std::vector<ManifoldPoint>> groups;
};

}

// Ranks points by predicted separation. Spread is chosen among points within band of the deepest,
// then free slots take the deepest of the rest.
void reduceManifold(std::vector<ManifoldPoint>& pts, const Vec3& normal, float band, Manifold& m, float tie) {
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
    for (size_t i = 0; i < i0; ++i) {
        if (pts[i].predicted <= pts[i0].predicted + tie) {
            i0 = i;
            break;
        }
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

void collideMesh(const Shape& mesh, const Transform& xm, const Shape& other, const Transform& xo, float margin,
                 const Tolerances& tol, const Motion& motion, std::vector<Manifold>& out) {
    out.clear();
    const TriangleMesh::Data* d = meshOf(mesh);
    if (!d || other.kind == ShapeKind::HalfSpace || other.kind == ShapeKind::Mesh) { return; }
    thread_local MeshBuffers buf;
    buf.found.clear();
    const Transform rel{ xm.q.unrotate(xo.p - xm.p), xm.q.conjugate() * xo.q };
    const float absTol = 0.5f * tol.linearSlop;

    if (other.kind == ShapeKind::Sphere) {
        const Vec3 center = rel.p;
        const float r = other.radius;
        const Vec3 reach{ r + margin, r + margin, r + margin };
        d->query(center - reach, center + reach, [&](int t) {
            Vec3 c[3];
            d->corners(t, c);
            const Vec3 n = d->normals[t];
            int feature;
            const Vec3 q = closestOnTriangle(c, center, feature);
            const float h = dot(n, center - c[0]);
            const int side = h >= 0.0f ? 0 : 1;
            const float dist = length(center - q);
            Vec3 normal, surface;
            float sep;
            const bool own = feature != 6 && dist > 0.0f && ownNormal(d->edges[t], feature, side);
            if (!own && dist - std::fabs(h) > tol.linearSlop) { return; }
            if (!own) {
                normal = side == 0 ? n : -n;
                sep = std::fabs(h) - r;
                surface = center - normal * std::fabs(h);
            } else {
                normal = (center - q) / dist;
                sep = dist - r;
                surface = q;
            }
            if (sep > margin) { return; }
            Manifold m;
            m.normal = xm.q.rotate(normal);
            const Vec3 w = xm.apply((surface + center - normal * r) * 0.5f);
            m.points[0] = { w, sep, motion.predict(w, m.normal, sep) };
            m.count = 1;
            buf.found.push_back(m);
        });
    } else if (other.kind == ShapeKind::Cylinder) {
        const Solid b = cylinderSolid(other, rel);
        const Vec3 reach = cylinderReach(b.axis, b.radius, b.halfHeight) + Vec3{ margin, margin, margin };
        const Convex convexB = convexOf(other);
        d->query(b.center - reach, b.center + reach, [&](int t) {
            Vec3 c[3];
            d->corners(t, c);
            const Vec3 n = d->normals[t];
            buf.tri.set(c, n);
            Solid a;
            a.poly = ownFrame(buf.tri.poly);
            const Convex triangle{ buf.tri.poly.vertices.data(), 3, 0.0f, false };
            Separation best;
            if (!separate(a, b, triangle, xm, convexB, xo, xm.q, margin, tol, best)) { return; }
            Manifold m;
            solidContact(a, b, best.normal, best.gap, xm, margin, tol.linearSlop, motion, buf.round, m);
            if (m.count == 0) { return; }

            // An edge or corner without a normal of its own leaves the push to the face.
            if (std::fabs(dot(best.normal, n)) < 0.9999f) {
                Separation front, back;
                tryAxis(a, b, n, front);
                tryAxis(a, b, -n, back);
                const int side = front.gap >= back.gap ? 0 : 1;
                if (!keepsNormal(d->edges[t], triangleFeature(c, best.normal, tol.linearSlop), side, best.normal, n)) {
                    const Separation& face = side == 0 ? front : back;
                    m.count = 0;
                    solidContact(a, b, face.normal, face.gap, xm, margin, tol.linearSlop, motion, buf.round, m);
                    for (int k = 0; k < m.count; ++k) {
                        const Vec3 p = xm.applyInverse(m.points[k].position);
                        if (m.points[k].separation < best.gap - tol.linearSlop && !overTriangle(c, n, p, tol.linearSlop)) {
                            return;
                        }
                    }
                    if (m.count == 0) { return; }
                }
            }
            buf.found.push_back(m);
        });
    } else {
        const FramePoly b = otherFrame(other.poly, rel, buf.frame);
        const Convex convexB = convexOf(other);
        const Vec3 grow{ margin, margin, margin };
        d->query(b.lo - grow, b.hi + grow, [&](int t) {
            Vec3 c[3];
            d->corners(t, c);
            const Vec3 n = d->normals[t];
            buf.tri.set(c, n);
            const FramePoly a = ownFrame(buf.tri.poly);

            const FaceQuery fa = queryFaces(a, b, margin);
            if (fa.separation > margin) { return; }
            const FaceQuery fb = queryFaces(b, a, margin);
            if (fb.separation > margin) { return; }
            const float reach = std::fmax(margin, -std::fmax(fa.separation, fb.separation));
            const EdgeQuery eq = queryTriangleEdges(a, b, reach, buf.near);
            if (eq.separation > margin) { return; }
            if (eq.pruned) {
                const Convex triangle{ buf.tri.poly.vertices.data(), 3, 0.0f, false };
                if (distance(triangle, xm, convexB, xo).distance > margin) { return; }
            }

            Manifold m;
            if (eq.a >= 0 && eq.separation > 0.9f * std::fmax(fa.separation, fb.separation) + absTol) {
                edgeContact(a, b, eq, xm, margin, tol.linearSlop, motion, buf.clip, m);
            } else if (fb.separation > 0.98f * fa.separation + absTol) {
                faceContact(b, fb.index, a, true, xm, margin, tol.linearSlop, motion.swapped(), buf.clip, m);
            } else {
                faceContact(a, fa.index, b, false, xm, margin, tol.linearSlop, motion, buf.clip, m);
            }
            if (m.count == 0) { return; }

            // Off the face normal on an edge or corner that has no normal of its own, the neighbor takes
            // that direction, so this triangle only pushes along its face. Where the face is deeper than
            // the shapes overlap, the contact belongs to the neighbor.
            const Vec3 normal = xm.q.unrotate(m.normal);
            if (std::fabs(dot(normal, n)) < 0.9999f) {
                if (!keepsNormal(d->edges[t], triangleFeature(c, normal, tol.linearSlop), fa.index, normal, n)) {
                    const float floor = std::fmax(std::fmax(fa.separation, fb.separation), eq.separation) - tol.linearSlop;
                    m.count = 0;
                    faceContact(a, fa.index, b, false, xm, margin, tol.linearSlop, motion, buf.clip, m);
                    for (int k = 0; k < m.count; ++k) {
                        const Vec3 p = xm.applyInverse(m.points[k].position);
                        if (m.points[k].separation < floor && !overTriangle(c, n, p, tol.linearSlop)) { return; }
                    }
                    if (m.count == 0) { return; }
                }
            }
            for (int k = 0; k < m.count; ++k) {
                ManifoldPoint& p = m.points[k];
                p.predicted = motion.predict(p.position, m.normal, p.separation);
            }
            buf.found.push_back(m);
        });
    }

    // Triangles facing the same way share one manifold, under the normal of the deepest of them.
    if (buf.groups.size() < buf.found.size()) { buf.groups.resize(buf.found.size()); }
    std::vector<int>& owner = buf.owner;
    owner.assign(buf.found.size(), -1);
    int groupCount = 0;
    for (size_t i = 0; i < buf.found.size(); ++i) {
        if (owner[i] >= 0) { continue; }
        std::vector<ManifoldPoint>& pts = buf.groups[groupCount];
        pts.clear();
        size_t lead = i;
        float deepest = FLT_MAX;
        for (size_t j = i; j < buf.found.size(); ++j) {
            if (owner[j] >= 0 || dot(buf.found[j].normal, buf.found[i].normal) < kSameNormal) { continue; }
            owner[j] = groupCount;
            for (int k = 0; k < buf.found[j].count; ++k) {
                pts.push_back(buf.found[j].points[k]);
                if (buf.found[j].points[k].separation < deepest) {
                    deepest = buf.found[j].points[k].separation;
                    lead = j;
                }
            }
        }
        // A curved patch is deepest at one end, so spread is chosen among all penetrating points.
        Manifold m;
        float deep = FLT_MAX;
        for (const ManifoldPoint& p : pts) { deep = std::fmin(deep, p.predicted); }
        reduceManifold(pts, buf.found[lead].normal, std::fmax(tol.linearSlop, -deep), m);
        out.push_back(m);
        ++groupCount;
    }
}

void collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
             float margin, const Tolerances& tol, const Motion& motion, Manifold& m) {
    m.count = 0;
    const ShapeKind ka = a.kind, kb = b.kind;
    if (ka == ShapeKind::HalfSpace && kb == ShapeKind::HalfSpace) { return; }
    if (ka == ShapeKind::Mesh || kb == ShapeKind::Mesh) { return; }

    if (ka == ShapeKind::HalfSpace) {
        collideHalfSpace(xa, b, xb, margin, tol.linearSlop, motion, m);
    } else if (kb == ShapeKind::HalfSpace) {
        collideHalfSpace(xb, a, xa, margin, tol.linearSlop, motion.swapped(), m);
        m.normal = -m.normal;
    } else if (ka == ShapeKind::Sphere && kb == ShapeKind::Sphere) {
        collideSpheres(xa.p, a.radius, xb.p, b.radius, margin, m);
    } else if (ka == ShapeKind::Sphere) {
        if (kb == ShapeKind::Cylinder) {
            collideSphereCylinder(b, xb, xa.p, a.radius, margin, m);
        } else {
            collideSpherePoly(b, xb, xa.p, a.radius, margin, m);
        }
        m.normal = -m.normal;
    } else if (kb == ShapeKind::Sphere) {
        if (ka == ShapeKind::Cylinder) {
            collideSphereCylinder(a, xa, xb.p, b.radius, margin, m);
        } else {
            collideSpherePoly(a, xa, xb.p, b.radius, margin, m);
        }
    } else if (ka == ShapeKind::Cylinder || kb == ShapeKind::Cylinder) {
        collideRound(a, xa, b, xb, margin, tol, motion, m);
    } else {
        collideConvex(a, xa, b, xb, margin, tol, motion, m);
    }
}

}
