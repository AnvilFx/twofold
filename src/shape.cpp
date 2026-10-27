#include <algorithm>

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

}

Vec3 cylinderReach(const Vec3& axis, float radius, float halfHeight) {
    Vec3 r;
    for (int i = 0; i < 3; ++i) {
        r[i] = halfHeight * std::fabs(axis[i]) + radius * std::sqrt(std::fmax(1.0f - axis[i] * axis[i], 0.0f));
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
    if (const Sphere* sp = std::get_if<Sphere>(&s.geometry)) {
        s.kind = ShapeKind::Sphere;
        s.radius = sp->radius;
        s.extent = sp->radius;
    } else if (const Box* bx = std::get_if<Box>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        s.poly = makeBox(bx->halfExtents);
    } else if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        s.kind = ShapeKind::Cylinder;
        s.radius = cy->radius;
        s.halfHeight = cy->halfHeight;
        s.extent = std::sqrt(cy->radius * cy->radius + cy->halfHeight * cy->halfHeight);
    } else if (Prism* pr = std::get_if<Prism>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        pr->polygon = counterClockwise(pr->polygon);
        if (pr->zMax < pr->zMin) { std::swap(pr->zMin, pr->zMax); }
        s.poly = makePrism(pr->polygon, pr->zMin, pr->zMax);
    } else if (const Hull* hu = std::get_if<Hull>(&s.geometry)) {
        HullMesh mesh = computeHull(hu->points, hu->maxVertices);
        if (mesh.vertices.empty()) {
            s.kind = ShapeKind::Sphere;
            s.extent = 0.0f;
        } else {
            s.kind = ShapeKind::Convex;
            s.poly = makePolyhedron(std::move(mesh.vertices), mesh.faces);
        }
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
