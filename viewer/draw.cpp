#include "draw.h"

#include <algorithm>
#include <cmath>

#include "rlgl.h"

using namespace twofold;

namespace {

constexpr float kPi = 3.14159265358979f;

const Vec3 kLight = normalize(Vec3{ 0.4f, 1.0f, 0.3f });
const Color kSelected{ 255, 220, 60, 255 };

Vector3 rl(const Vec3& v) { return { v.x, v.y, v.z }; }

Color scale(Color c, float k) {
    auto ch = [k](unsigned char v) { return (unsigned char)std::min(255.0f, (float)v * k); };
    return { ch(c.r), ch(c.g), ch(c.b), c.a };
}

Color shade(Color c, const Vec3& n) { return scale(c, 0.45f + 0.6f * std::fmax(0.0f, dot(n, kLight))); }

void color(Color c) { rlColor4ub(c.r, c.g, c.b, c.a); }
void vertex(const Vec3& p) { rlVertex3f(p.x, p.y, p.z); }

void triangle(const Vec3& a, const Vec3& b, const Vec3& c, Color col) {
    color(col);
    vertex(a);
    vertex(b);
    vertex(c);
}

std::vector<Vec2> counterClockwise(std::vector<Vec2> poly) {
    float area = 0.0f;
    for (size_t i = 0; i < poly.size(); ++i) { area += cross(poly[i], poly[(i + 1) % poly.size()]); }
    if (area < 0.0f) { std::reverse(poly.begin(), poly.end()); }
    return poly;
}

std::vector<Vec2> circle(float radius, int sides) {
    std::vector<Vec2> p;
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        p.push_back({ radius * std::cos(a), radius * std::sin(a) });
    }
    return p;
}

void drawPrism(const twofold::Transform& xf, const std::vector<Vec2>& input, float z0, float z1, const Color* fill,
               Color edge, bool verticals) {
    const std::vector<Vec2> poly = counterClockwise(input);
    const size_t n = poly.size();
    if (n < 3) { return; }
    std::vector<Vec3> lo(n), hi(n);
    for (size_t i = 0; i < n; ++i) {
        lo[i] = xf.apply({ poly[i].x, poly[i].y, z0 });
        hi[i] = xf.apply({ poly[i].x, poly[i].y, z1 });
    }

    if (fill) {
        rlCheckRenderBatchLimit((int)(n * 12));
        rlBegin(RL_TRIANGLES);
        for (size_t i = 0; i < n; ++i) {
            const size_t j = (i + 1) % n;
            const Vec2 e = poly[j] - poly[i];
            const Color c = shade(*fill, xf.q.rotate(normalize(Vec3{ e.y, -e.x, 0.0f })));
            triangle(lo[i], lo[j], hi[j], c);
            triangle(lo[i], hi[j], hi[i], c);
        }
        const Color top = shade(*fill, xf.q.rotate({ 0, 0, 1 }));
        const Color bottom = shade(*fill, xf.q.rotate({ 0, 0, -1 }));
        for (size_t i = 1; i + 1 < n; ++i) {
            triangle(hi[0], hi[i], hi[i + 1], top);
            triangle(lo[0], lo[i + 1], lo[i], bottom);
        }
        rlEnd();
    }

    rlCheckRenderBatchLimit((int)(n * 6));
    rlBegin(RL_LINES);
    color(edge);
    for (size_t i = 0; i < n; ++i) {
        const size_t j = (i + 1) % n;
        vertex(lo[i]);
        vertex(lo[j]);
        vertex(hi[i]);
        vertex(hi[j]);
        if (verticals) {
            vertex(lo[i]);
            vertex(hi[i]);
        }
    }
    rlEnd();
}

void drawSphere(const twofold::Transform& xf, float r, const Color* fill, Color edge) {
    const int rings = 12, slices = 24;
    if (fill) {
        auto point = [&](int i, int j, Vec3& normal) {
            const float th = kPi * (float)i / (float)rings;
            const float ph = 2.0f * kPi * (float)j / (float)slices;
            const Vec3 u{ std::sin(th) * std::cos(ph), std::cos(th), std::sin(th) * std::sin(ph) };
            normal = xf.q.rotate(u);
            return xf.apply(u * r);
        };
        rlCheckRenderBatchLimit(rings * slices * 6);
        rlBegin(RL_TRIANGLES);
        for (int i = 0; i < rings; ++i) {
            for (int j = 0; j < slices; ++j) {
                const Color c = ((j / 3) % 2 == 0) ? *fill : scale(*fill, 0.75f);
                Vec3 na, nb, nc, nd;
                const Vec3 a = point(i, j, na), b = point(i + 1, j, nb);
                const Vec3 cc = point(i + 1, j + 1, nc), d = point(i, j + 1, nd);
                color(shade(c, na)); vertex(a);
                color(shade(c, nc)); vertex(cc);
                color(shade(c, nb)); vertex(b);
                color(shade(c, na)); vertex(a);
                color(shade(c, nd)); vertex(d);
                color(shade(c, nc)); vertex(cc);
            }
        }
        rlEnd();
        return;
    }
    rlCheckRenderBatchLimit(3 * slices * 2);
    rlBegin(RL_LINES);
    color(edge);
    for (int axis = 0; axis < 3; ++axis) {
        for (int j = 0; j < slices; ++j) {
            for (int k = j; k <= j + 1; ++k) {
                const float a = 2.0f * kPi * (float)k / (float)slices;
                const float c = r * std::cos(a), s = r * std::sin(a);
                const Vec3 u = axis == 0 ? Vec3{ 0, c, s } : axis == 1 ? Vec3{ c, 0, s } : Vec3{ c, s, 0 };
                vertex(xf.apply(u));
            }
        }
    }
    rlEnd();
}

void drawHalfSpace(const twofold::Transform& xf, const Vec3& focus, bool filled) {
    const int half = 30;
    const Vec3 local = xf.applyInverse(focus);
    const float cx = std::floor(local.x), cy = std::floor(local.y);
    const Vec3 n = xf.q.rotate({ 0, 0, 1 });
    if (!filled) {
        rlCheckRenderBatchLimit(half * 8);
        rlBegin(RL_LINES);
        color({ 90, 96, 108, 255 });
        for (int i = -half; i <= half; ++i) {
            vertex(xf.apply({ cx + (float)i, cy - half, 0 }));
            vertex(xf.apply({ cx + (float)i, cy + half, 0 }));
            vertex(xf.apply({ cx - half, cy + (float)i, 0 }));
            vertex(xf.apply({ cx + half, cy + (float)i, 0 }));
        }
        rlEnd();
        return;
    }
    for (int i = -half; i < half; ++i) {
        rlCheckRenderBatchLimit(half * 2 * 6);
        rlBegin(RL_TRIANGLES);
        for (int j = -half; j < half; ++j) {
            const float x0 = cx + (float)i, y0 = cy + (float)j;
            const bool dark = ((int)(x0 + 10000) + (int)(y0 + 10000)) % 2 == 0;
            const Color c = shade(dark ? Color{ 58, 62, 70, 255 } : Color{ 72, 77, 86, 255 }, n);
            const Vec3 a = xf.apply({ x0, y0, 0 }), b = xf.apply({ x0 + 1, y0, 0 });
            const Vec3 d = xf.apply({ x0 + 1, y0 + 1, 0 }), e = xf.apply({ x0, y0 + 1, 0 });
            triangle(a, b, d, c);
            triangle(a, d, e, c);
        }
        rlEnd();
    }
}

void drawHull(const twofold::Transform& xf, const HullMesh& hull, const Color* fill, Color edge) {
    std::vector<Vec3> p;
    for (const Vec3& v : hull.vertices) { p.push_back(xf.apply(v)); }
    int corners = 0;
    for (const std::vector<int>& f : hull.faces) { corners += (int)f.size(); }
    if (fill) {
        rlCheckRenderBatchLimit(corners * 3);
        rlBegin(RL_TRIANGLES);
        for (const std::vector<int>& f : hull.faces) {
            const Color c = shade(*fill, normalize(cross(p[f[1]] - p[f[0]], p[f[2]] - p[f[0]])));
            for (size_t i = 1; i + 1 < f.size(); ++i) { triangle(p[f[0]], p[f[i]], p[f[i + 1]], c); }
        }
        rlEnd();
    }
    rlCheckRenderBatchLimit(corners * 2);
    rlBegin(RL_LINES);
    color(edge);
    for (const std::vector<int>& f : hull.faces) {
        for (size_t i = 0; i < f.size(); ++i) {
            vertex(p[f[i]]);
            vertex(p[f[(i + 1) % f.size()]]);
        }
    }
    rlEnd();
}

void drawGeometry(const twofold::Transform& xf, const Visual& v, const Color* fill, Color edge, const Vec3& focus) {
    const Geometry& g = v.geometry;
    if (const Box* b = std::get_if<Box>(&g)) {
        const Vec3 h = b->halfExtents;
        drawPrism(xf, { { -h.x, -h.y }, { h.x, -h.y }, { h.x, h.y }, { -h.x, h.y } }, -h.z, h.z, fill, edge, true);
    } else if (const Sphere* sp = std::get_if<Sphere>(&g)) {
        drawSphere(xf, sp->radius, fill, edge);
    } else if (const Cylinder* cy = std::get_if<Cylinder>(&g)) {
        drawPrism(xf, circle(cy->radius, 32), -cy->halfHeight, cy->halfHeight, fill, edge, false);
    } else if (const Prism* pr = std::get_if<Prism>(&g)) {
        drawPrism(xf, pr->polygon, pr->zMin, pr->zMax, fill, edge, true);
    } else if (std::holds_alternative<Hull>(g)) {
        drawHull(xf, v.hull, fill, edge);
    } else {
        drawHalfSpace(xf, focus, fill != nullptr);
    }
}

Color hue(float degrees, float saturation, float value) { return ColorFromHSV(std::fmod(degrees, 360.0f), saturation, value); }

Color bodyColor(const World& w, const Visual& v, ColorMode mode) {
    const BodyId b = v.body;
    const MotionType type = w.motionType(b);
    const bool asleep = type == MotionType::Dynamic && !w.isAwake(b);
    Color c = v.color;
    switch (mode) {
        case ColorMode::Material: break;
        case ColorMode::Motion:
            if (type == MotionType::Static) {
                c = { 150, 150, 158, 255 };
            } else if (type == MotionType::Kinematic) {
                c = { 90, 200, 110, 255 };
            } else {
                c = asleep ? Color{ 90, 100, 130, 255 } : Color{ 235, 120, 60, 255 };
            }
            return c;
        case ColorMode::Island:
            if (type != MotionType::Dynamic) { return { 150, 150, 158, 255 }; }
            c = hue((float)w.island(b) * 137.508f, 0.65f, 0.95f);
            break;
        case ColorMode::Speed: {
            if (type == MotionType::Static) { return { 150, 150, 158, 255 }; }
            const float speed = length(w.linearVelocity(b));
            c = hue(220.0f * (1.0f - std::fmin(speed / 10.0f, 1.0f)), 0.75f, 0.95f);
            break;
        }
    }
    return asleep ? scale(c, 0.55f) : c;
}

bool include(const DrawOptions& o, BodyId body, BodyId selected) { return !o.selectedOnly || body == selected; }

void drawContacts(const Scene& s, const DrawOptions& o, BodyId selected) {
    const World& w = s.world;
    const float toForce = 1.0f / s.dt;
    for (const ContactInfo& c : w.contacts(o.speculative)) {
        if (o.selectedOnly && c.bodyA != selected && c.bodyB != selected) { continue; }
        for (int i = 0; i < c.pointCount; ++i) {
            const Vec3 p = c.points[i];
            const bool touching = c.separations[i] < 0.0f || c.normalImpulses[i] > 0.0f;
            DrawSphereEx(rl(p), 0.04f, 4, 6, touching ? RED : Color{ 150, 150, 160, 255 });
            const float len = o.forces ? c.normalImpulses[i] * toForce * o.forceScale : 0.3f;
            drawLine(p, p + c.normal * len, touching ? YELLOW : Color{ 120, 120, 130, 255 });
            if (o.friction) { drawLine(p, p + c.tangentImpulses[i] * (toForce * o.forceScale), SKYBLUE); }
        }
    }
}

void drawPlanes(const Scene& s, const Vec3& focus) {
    const World& w = s.world;
    struct Plane {
        Vec3 n;
        float d;
        Quat frame;
    };
    std::vector<Plane> planes;
    for (BodyId b : s.bodies) {
        if (!w.isValid(b) || (w.locks(b) & LockPlanar) != LockPlanar) { continue; }
        const Quat frame = w.lockFrame(b);
        const Vec3 n = frame.rotate({ 0, 0, 1 });
        const float d = dot(n, w.centerOfMass(b));
        const bool known = std::any_of(planes.begin(), planes.end(), [&](const Plane& p) {
            return dot(p.n, n) > 0.9999f && std::fabs(p.d - d) < 1e-3f;
        });
        if (!known) { planes.push_back({ n, d, frame }); }
    }
    const int half = 12;
    for (const Plane& p : planes) {
        const Vec3 u = p.frame.rotate({ 1, 0, 0 }), v = p.frame.rotate({ 0, 1, 0 });
        const Vec3 center = focus - p.n * (dot(p.n, focus) - p.d);
        const float cu = std::floor(dot(center, u)), cv = std::floor(dot(center, v));
        const Vec3 origin = p.n * p.d;
        for (int i = -half; i <= half; ++i) {
            const Color c{ 245, 150, 60, (unsigned char)(i == 0 ? 160 : 70) };
            drawLine(origin + u * (cu + (float)i) + v * (cv - half), origin + u * (cu + (float)i) + v * (cv + half), c);
            drawLine(origin + v * (cv + (float)i) + u * (cu - half), origin + v * (cv + (float)i) + u * (cu + half), c);
        }
    }
}

}

void drawLine(const Vec3& a, const Vec3& b, Color color) { DrawLine3D(rl(a), rl(b), color); }

void drawScene(const Scene& s, const DrawOptions& o, const Vec3& focus, BodyId selected) {
    const World& w = s.world;
    for (const Visual& v : s.visuals) {
        if (!w.isValid(v.body)) { continue; }
        const twofold::Transform xf = w.transform(v.body) * v.local;
        const Color c = bodyColor(w, v, o.colorMode);
        if (o.wireframe) {
            drawGeometry(xf, v, nullptr, c, focus);
        } else {
            drawGeometry(xf, v, &c, scale(c, 0.45f), focus);
        }
    }

    if (o.planes) { drawPlanes(s, focus); }

    for (const Visual& v : s.visuals) {
        if (!w.isValid(v.body) || !include(o, v.body, selected)) { continue; }
        if (o.bounds && !std::holds_alternative<HalfSpace>(v.geometry)) {
            const Aabb box = w.bounds(v.shape);
            DrawBoundingBox({ rl(box.lower), rl(box.upper) }, { 120, 220, 120, 255 });
        }
    }

    for (BodyId b : s.bodies) {
        if (!w.isValid(b) || !include(o, b, selected) || w.motionType(b) == MotionType::Static) { continue; }
        const Vec3 com = w.centerOfMass(b);
        if (o.centers) {
            const Quat q = w.rotation(b);
            DrawSphereEx(rl(com), 0.05f, 4, 6, WHITE);
            drawLine(com, com + q.rotate({ 0.5f, 0, 0 }), RED);
            drawLine(com, com + q.rotate({ 0, 0.5f, 0 }), GREEN);
            drawLine(com, com + q.rotate({ 0, 0, 0.5f }), BLUE);
        }
        if (o.velocities) {
            drawLine(com, com + w.linearVelocity(b) * 0.2f, MAGENTA);
            drawLine(com, com + w.angularVelocity(b) * 0.2f, ORANGE);
        }
    }

    if (o.contacts) { drawContacts(s, o, selected); }

    if (w.isValid(selected)) {
        rlDrawRenderBatchActive();
        rlDisableDepthTest();
        for (const Visual& v : s.visuals) {
            if (v.body != selected || std::holds_alternative<HalfSpace>(v.geometry)) { continue; }
            drawGeometry(w.transform(v.body) * v.local, v, nullptr, kSelected, focus);
        }
        rlDrawRenderBatchActive();
        rlEnableDepthTest();
    }
}
