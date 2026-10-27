#include <cmath>

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

}
