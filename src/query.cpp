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

bool rayPolyhedron(const Polyhedron& poly, const Vec3& o, const Vec3& d, float maxFraction, float& fraction, Vec3& normal) {
    float lower = 0.0f, upper = maxFraction;
    int index = -1;
    for (int i = 0; i < (int)poly.faces.size(); ++i) {
        const Polyhedron::Face& f = poly.faces[i];
        const float num = f.offset - dot(f.normal, o);
        const float den = dot(f.normal, d);
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
    normal = poly.faces[index].normal;
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
    }
    if (hit) { normal = xf.q.rotate(normal); }
    return hit;
}

}
