#include <algorithm>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kPi = 3.14159265358979f;
constexpr int kCylinderSides = 24;

std::vector<Vec2> counterClockwise(std::vector<Vec2> poly) {
    float area = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { area += cross(poly[i], poly[(i + 1) % poly.size()]); }
    if (area < 0.0f) { std::reverse(poly.begin(), poly.end()); }
    return poly;
}

std::vector<Vec2> regularPolygon(float radius, int sides) {
    std::vector<Vec2> poly;
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        poly.push_back({ radius * std::cos(a), radius * std::sin(a) });
    }
    return poly;
}

}

void buildShape(Shape& s) {
    s.poly = Polyhedron{};
    s.radius = 0.0f;
    if (const Sphere* sp = std::get_if<Sphere>(&s.geometry)) {
        s.kind = ShapeKind::Sphere;
        s.radius = sp->radius;
        s.extent = sp->radius;
    } else if (const Box* bx = std::get_if<Box>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        s.poly = makeBox(bx->halfExtents);
    } else if (const Cylinder* cy = std::get_if<Cylinder>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        s.poly = makePrism(regularPolygon(cy->radius, kCylinderSides), -cy->halfHeight, cy->halfHeight);
    } else if (Prism* pr = std::get_if<Prism>(&s.geometry)) {
        s.kind = ShapeKind::Convex;
        pr->polygon = counterClockwise(pr->polygon);
        if (pr->zMax < pr->zMin) { std::swap(pr->zMin, pr->zMax); }
        s.poly = makePrism(pr->polygon, pr->zMin, pr->zMax);
    } else {
        s.kind = ShapeKind::HalfSpace;
        s.extent = 0.0f;
    }
    if (s.kind == ShapeKind::Convex) {
        float e = 0.0f;
        for (const Vec3& v : s.poly.vertices) { e = std::max(e, length(v)); }
        s.extent = e;
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
    lo = hi = xf.apply(s.poly.vertices[0]);
    for (const Vec3& v : s.poly.vertices) {
        const Vec3 w = xf.apply(v);
        lo = minVec(lo, w);
        hi = maxVec(hi, w);
    }
}

}
