#include <cmath>
#include <cstdint>
#include <vector>

#include "scene.h"

using namespace twofold;

namespace {

constexpr float kPi = 3.14159265358979f;

const Quat kUp = Quat::fromTo({ 0, 0, 1 }, { 0, 1, 0 });

const Color kStatic{ 150, 150, 158, 255 };
const Color kKinematic{ 90, 200, 110, 255 };
const Color kSolid[] = { { 66, 135, 245, 255 }, { 86, 180, 233, 255 }, { 120, 110, 230, 255 }, { 70, 190, 170, 255 } };
const Color kFlat[] = { { 245, 150, 60, 255 }, { 240, 105, 80, 255 }, { 250, 190, 70, 255 } };

struct Rng {
    uint32_t s = 0x9e3779b9u;

    float next() {
        s ^= s << 13;
        s ^= s >> 17;
        s ^= s << 5;
        return (float)(s & 0xffffff) / 16777216.0f;
    }
    float range(float lo, float hi) { return lo + (hi - lo) * next(); }
};

Color solid(int i) { return kSolid[i % 4]; }
Color flat(int i) { return kFlat[i % 3]; }

ShapeDesc shape(const Geometry& g, float friction = 0.6f, float restitution = 0.0f) {
    ShapeDesc sd;
    sd.geometry = g;
    sd.friction = friction;
    sd.restitution = restitution;
    return sd;
}

BodyDesc body(const Vec3& p, MotionType type = MotionType::Dynamic, const Quat& q = {}) {
    BodyDesc bd;
    bd.type = type;
    bd.position = p;
    bd.rotation = q;
    return bd;
}

BodyDesc planar(const Vec3& p, const Quat& q = {}) {
    BodyDesc bd = body(p, MotionType::Dynamic, q);
    bd.locks = LockPlanar;
    return bd;
}

void ground(Scene& s, float friction = 0.6f) {
    s.add(body({}, MotionType::Static, kUp), shape(HalfSpace{}, friction), kStatic);
}

std::vector<Vec2> regular(int sides, float radius) {
    std::vector<Vec2> p;
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        p.push_back({ radius * std::cos(a), radius * std::sin(a) });
    }
    return p;
}

void tower(Scene& s) {
    ground(s);
    for (int i = 0; i < 14; ++i) { s.add(body({ 0, 0.5f + (float)i, 0 }), shape(Box{}), solid(i)); }
}

void pyramid(Scene& s) {
    ground(s);
    const int base = 14;
    int k = 0;
    for (int row = 0; row < base; ++row) {
        const int n = base - row;
        for (int i = 0; i < n; ++i) {
            const float x = ((float)i - 0.5f * (float)(n - 1)) * 1.02f;
            s.add(body({ x, 0.5f + (float)row, 0 }), shape(Box{}), solid(k++));
        }
    }
}

void rain(Scene& s) {
    ground(s);
    s.add(body({ -2, 1, 0 }, MotionType::Static, Quat::axisAngle({ 0, 0, 1 }, 0.4f)), shape(Box{ { 3, 0.3f, 2 } }), kStatic);
    s.add(body({ 3, 1, 1 }, MotionType::Static, kUp), shape(Cylinder{ 1.2f, 1.0f }), kStatic);
    Rng rng;
    for (int i = 0; i < 180; ++i) {
        const int cell = i % 12, level = i / 12;
        const Vec3 p{ -3.0f + 2.0f * (float)(cell % 4) + (float)(level % 2) + rng.range(-0.1f, 0.1f), 5.0f + 1.8f * (float)level,
                      -2.0f + 2.0f * (float)(cell / 4) + rng.range(-0.1f, 0.1f) };
        const Quat q = Quat::axisAngle({ rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1) }, rng.range(0, kPi));
        Geometry g;
        switch (i % 4) {
            case 0: g = Box{ { rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f), rng.range(0.2f, 0.5f) } }; break;
            case 1: g = Sphere{ rng.range(0.2f, 0.45f) }; break;
            case 2: g = Cylinder{ rng.range(0.2f, 0.4f), rng.range(0.15f, 0.4f) }; break;
            default: g = Prism{ regular(3 + i % 5, rng.range(0.25f, 0.45f)), -0.2f, 0.2f }; break;
        }
        s.add(body(p, MotionType::Dynamic, q), shape(g), solid(i));
    }
}

std::vector<Vec3> rock(Rng& rng, const Vec3& size, int count) {
    std::vector<Vec3> p;
    for (int i = 0; i < count; ++i) {
        const Vec3 d = normalize(Vec3{ rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1) });
        const float r = rng.range(0.8f, 1.0f);
        p.push_back({ d.x * size.x * r, d.y * size.y * r, d.z * size.z * r });
    }
    return p;
}

std::vector<Vec3> gem(int sides, float radius, float crown, float pavilion) {
    std::vector<Vec3> p{ { 0, -pavilion, 0 } };
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        p.push_back({ radius * std::cos(a), 0, radius * std::sin(a) });
        p.push_back({ 0.6f * radius * std::cos(a + kPi / (float)sides), crown, 0.6f * radius * std::sin(a + kPi / (float)sides) });
    }
    return p;
}

std::vector<Vec3> cone(int sides, float radius, float height) {
    std::vector<Vec3> p{ { 0, height, 0 } };
    for (int i = 0; i < sides; ++i) {
        const float a = 2.0f * kPi * (float)i / (float)sides;
        p.push_back({ radius * std::cos(a), 0, radius * std::sin(a) });
    }
    return p;
}

void hulls(Scene& s) {
    ground(s);
    Rng rng;
    s.add(body({ 0, 0, 0 }, MotionType::Static), shape(Hull{ rock(rng, { 3, 1.6f, 2.5f }, 200) }), kStatic);
    for (int i = 0; i < 96; ++i) {
        const int cell = i % 16, level = i / 16;
        const Vec3 p{ -3.0f + 2.0f * (float)(cell % 4) + (float)(level % 2) + rng.range(-0.1f, 0.1f), 4.0f + 1.6f * (float)level,
                      -3.0f + 2.0f * (float)(cell / 4) + rng.range(-0.1f, 0.1f) };
        const Quat q = Quat::axisAngle({ rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1) }, rng.range(0, kPi));
        Hull h;
        switch (i % 3) {
            case 0: h.points = rock(rng, { rng.range(0.3f, 0.6f), rng.range(0.2f, 0.4f), rng.range(0.3f, 0.5f) }, 120); break;
            case 1: h.points = gem(8, rng.range(0.3f, 0.45f), 0.15f, 0.35f); break;
            default: h.points = cone(16, rng.range(0.25f, 0.35f), rng.range(0.5f, 0.8f)); break;
        }
        h.maxVertices = 32;
        s.add(body(p, MotionType::Dynamic, q), shape(h), solid(i));
    }
}

// Rasterizes a signed distance in pixels, centered with y up, then cuts the traced outline into prisms.
template <typename F>
std::vector<Prism> cutout(int width, int height, float scale, float depth, F distance) {
    std::vector<uint8_t> alpha((size_t)width * height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float d = distance(Vec2{ (float)x + 0.5f - 0.5f * (float)width, 0.5f * (float)height - (float)y - 0.5f });
            alpha[(size_t)y * width + x] = (uint8_t)std::lround(std::fmin(std::fmax(0.5f - d, 0.0f), 1.0f) * 255.0f);
        }
    }
    std::vector<Prism> prisms;
    for (const std::vector<Vec2>& piece : alphaPieces({ alpha.data(), width, height, width, 1 })) {
        Prism pr{ {}, -depth, depth };
        for (const Vec2& q : piece) { pr.polygon.push_back({ (q.x - 0.5f * (float)width) * scale, (0.5f * (float)height - q.y) * scale }); }
        prisms.push_back(std::move(pr));
    }
    return prisms;
}

BodyId addCutout(Scene& s, const BodyDesc& bd, const std::vector<Prism>& prisms, Color color) {
    const BodyId b = s.add(bd, shape(prisms[0], 0.5f), color);
    for (size_t i = 1; i < prisms.size(); ++i) { s.attach(b, shape(prisms[i], 0.5f), color); }
    return b;
}

void cutouts(Scene& s) {
    ground(s);
    auto ring = [](Vec2 p) { return std::fmax(length(p) - 40.0f, 24.0f - length(p)); };
    auto gear = [](Vec2 p) {
        const float r = length(p), teeth = std::cos(8.0f * std::atan2(p.y, p.x)) > 0.0f ? 44.0f : 36.0f;
        return std::fmax(r - teeth, 12.0f - r);
    };
    auto star = [](Vec2 p) {
        const float a = std::atan2(p.y, p.x);
        return length(p) - (18.0f + 26.0f * std::pow(std::fabs(std::cos(2.5f * a + 0.5f * kPi)), 4.0f));
    };
    auto crescent = [](Vec2 p) { return std::fmax(length(p) - 42.0f, 32.0f - length(p - Vec2{ 16, 6 })); };
    auto hook = [](Vec2 p) {
        const float r = length(p);
        return std::fabs(std::atan2(p.y, p.x) + 0.3f) < 0.75f ? 1.0f : std::fmax(r - 40.0f, 26.0f - r);
    };
    auto tee = [](Vec2 p) {
        const float bar = std::fmax(std::fabs(p.x) - 42.0f, std::fabs(p.y - 30.0f) - 10.0f);
        const float stem = std::fmax(std::fabs(p.x) - 10.0f, std::fabs(p.y + 6.0f) - 36.0f);
        return std::fmin(bar, stem);
    };
    const float scale = 0.0125f, depth = 0.05f;
    const std::vector<std::vector<Prism>> shapes{ cutout(96, 96, scale, depth, ring),     cutout(96, 96, scale, depth, gear),
                                                  cutout(96, 96, scale, depth, star),     cutout(96, 96, scale, depth, crescent),
                                                  cutout(96, 96, scale, depth, hook),     cutout(96, 96, scale, depth, tee) };

    const std::vector<Prism> bowl = cutout(400, 140, 0.025f, 0.6f, [](Vec2 p) {
        const float outer = std::fmax(std::fabs(p.x) - 190.0f, std::fabs(p.y) - 60.0f);
        const float hollow = std::fmax(std::fabs(p.x) - 96.0f - 0.6f * (p.y + 30.0f), -30.0f - p.y);
        return std::fmax(outer, -hollow);
    });
    addCutout(s, body({ 0, 1.5f, 0 }, MotionType::Static), bowl, kStatic);
    for (int i = 0; i < 4; ++i) {
        s.add(body({ -3.3f + 2.2f * (float)i, 6.0f, 0 }, MotionType::Static), shape(Box{ { 0.08f, 0.08f, 0.6f } }), kStatic);
    }

    Rng rng;
    for (int i = 0; i < 36; ++i) {
        const int cell = i % 6, level = i / 6;
        const Vec3 p{ -4.0f + 1.6f * (float)cell + 0.8f * (float)(level % 2) + rng.range(-0.1f, 0.1f), 8.0f + 1.6f * (float)level, 0 };
        addCutout(s, planar(p, Quat::axisAngle({ 0, 0, 1 }, rng.range(0, 2.0f * kPi))), shapes[(size_t)(i % 6)], flat(i));
    }
}

void dominoes(Scene& s) {
    ground(s);
    const float radius = 8.0f;
    const int count = 44;
    const float step = 0.85f / radius;
    BodyId first;
    for (int i = 0; i < count; ++i) {
        const float a = (float)i * step;
        const Quat q = Quat::axisAngle({ 0, 1, 0 }, -(a + 0.5f * kPi));
        const BodyId b = s.add(body({ radius * std::cos(a), 1.0f, radius * std::sin(a) }, MotionType::Dynamic, q),
                               shape(Box{ { 0.1f, 1.0f, 0.5f } }), solid(i));
        if (i == 0) { first = b; }
    }
    const Vec3 tangent{ 0, 0, 1 };
    s.world.addImpulseAtPoint(first, tangent * 0.6f, s.world.position(first) + Vec3{ 0, 0.9f, 0 });
}

void flatLayers(Scene& s) {
    ground(s);
    for (int i = 0; i < 5; ++i) {
        s.add(body({ -3.0f + 1.5f * (float)i, 0.25f + 0.5f * (float)i, 0.6f }, MotionType::Static),
              shape(Box{ { 0.75f, 0.25f + 0.5f * (float)i, 1.6f } }), kStatic);
    }
    s.add(body({ 5.0f, 2.0f, 0.6f }, MotionType::Static, Quat::axisAngle({ 0, 0, 1 }, -0.5f)),
          shape(Box{ { 2.0f, 0.15f, 1.6f } }), kStatic);

    Rng rng;
    for (int i = 0; i < 60; ++i) {
        const float z = (i % 3 == 0) ? 1.2f : 0.0f;
        const int cell = i % 4, level = i / 4;
        const Vec3 p{ -4.0f + 2.5f * (float)cell + 1.25f * (float)(level % 2) + rng.range(-0.3f, 0.3f), 5.0f + 0.9f * (float)level, z };
        const Quat q = Quat::axisAngle({ 0, 0, 1 }, rng.range(0, kPi));
        Geometry g;
        switch (i % 3) {
            case 0: g = Box{ { rng.range(0.2f, 0.45f), rng.range(0.15f, 0.35f), 0.04f } }; break;
            case 1: g = Prism{ regular(5 + i % 3, rng.range(0.25f, 0.4f)), -0.04f, 0.04f }; break;
            default: g = Cylinder{ rng.range(0.2f, 0.35f), 0.04f }; break;
        }
        s.add(planar(p, q), shape(g, 0.5f), flat(i));
    }
    for (int i = 0; i < 4; ++i) {
        BodyDesc ball = body({ -9.0f - 2.0f * (float)i, 4.5f, 0.0f });
        ball.linearVelocity = { 7, 0, 0 };
        ShapeDesc sd = shape(Sphere{ 0.35f });
        sd.density = 3.0f;
        s.add(ball, sd, solid(i));
    }
}

void sweeperUpdate(Scene& s, float dt) {
    const float speed = 0.8f;
    s.world.setKinematicTarget(s.driven[0], { 0, 0.5f, 0 }, Quat::axisAngle({ 0, 1, 0 }, (s.time + dt) * speed), dt);
}

void sweeper(Scene& s) {
    ground(s, 0.4f);
    BodyDesc paddle = body({ 0, 0.5f, 0 }, MotionType::Kinematic);
    const BodyId p = s.world.createBody(paddle);
    s.bodies.push_back(p);
    ShapeDesc arm = shape(Box{ { 3.5f, 0.45f, 0.15f } });
    arm.position = { 4.0f, 0, 0 };
    s.attach(p, arm, kKinematic);
    s.driven.push_back(p);
    s.update = sweeperUpdate;

    Rng rng;
    std::vector<Vec3> placed;
    for (int i = 0; i < 90; ++i) {
        Vec3 pos;
        bool clear = false;
        while (!clear) {
            const float a = rng.range(0, 2.0f * kPi), r = rng.range(1.5f, 7.0f);
            pos = { r * std::cos(a), 0.3f + 0.6f * (float)(i % 2), r * std::sin(a) };
            clear = pos.x < 0.0f || std::fabs(pos.z) > 0.6f;
            for (const Vec3& q : placed) { clear = clear && lengthSquared(q - pos) >= 0.81f; }
        }
        placed.push_back(pos);
        Geometry g = (i % 3 == 0) ? Geometry{ Sphere{ 0.3f } } : Geometry{ Box{ { 0.3f, 0.3f, 0.3f } } };
        s.add(body(pos), shape(g, 0.4f), solid(i));
    }
}

void slopes(Scene& s) {
    ground(s);
    const float frictions[] = { 0.02f, 0.25f, 0.9f };
    const Quat tilt = Quat::axisAngle({ 0, 0, 1 }, 0.35f);
    for (int i = 0; i < 3; ++i) {
        const Vec3 center{ 0, 2.0f, -3.0f + 3.0f * (float)i };
        s.add(body(center, MotionType::Static, tilt), shape(Box{ { 5, 0.25f, 1.2f } }, frictions[i]), kStatic);
        s.add(body(center + tilt.rotate({ 3.5f, 0.65f, 0 }), MotionType::Dynamic, tilt),
              shape(Box{ { 0.4f, 0.4f, 0.4f } }, frictions[i]), solid(i));
    }
    const float bounces[] = { 0.0f, 0.5f, 0.9f };
    for (int i = 0; i < 3; ++i) {
        s.add(body({ 9.0f, 5.0f, -3.0f + 3.0f * (float)i }), shape(Sphere{ 0.5f }, 0.6f, bounces[i]), solid(i + 1));
    }
}

void dormant(Scene& s) {
    ground(s);
    for (int y = 0; y < 3; ++y) {
        for (int x = 0; x < 6; ++x) {
            for (int z = 0; z < 3; ++z) {
                BodyDesc bd = body({ -3.0f + 1.2f * (float)x, 3.0f + 1.2f * (float)y, -1.2f + 1.2f * (float)z });
                bd.awake = false;
                s.add(bd, shape(Box{ { 0.45f, 0.45f, 0.45f } }), solid(x + y + z));
            }
        }
    }
}

}

BodyId Scene::add(const BodyDesc& desc, const ShapeDesc& shapeDesc, Color color) {
    const BodyId b = world.createBody(desc);
    bodies.push_back(b);
    attach(b, shapeDesc, color);
    return b;
}

void Scene::attach(BodyId b, const ShapeDesc& shapeDesc, Color color) {
    const ShapeId shape = world.addShape(b, shapeDesc);
    Visual v{ b, shape, shapeDesc.geometry, { shapeDesc.position, shapeDesc.rotation }, color, {} };
    if (const Hull* h = std::get_if<Hull>(&shapeDesc.geometry)) { v.hull = computeHull(h->points, h->maxVertices); }
    visuals.push_back(std::move(v));
}

const std::vector<SceneInfo>& sceneList() {
    static const std::vector<SceneInfo> scenes = {
        { "Tower", "Fourteen boxes stacked on a ground plane.", tower, { 0, 6, 0 }, 22, 0.6f, 0.25f },
        { "Pyramid", "A pyramid of 105 boxes.", pyramid, { 0, 6, 0 }, 26, 0.3f, 0.2f },
        { "Rain", "Boxes, spheres, cylinders and prisms falling onto obstacles.", rain, { 0, 3, 0 }, 22, 0.7f, 0.45f },
        { "Dominoes", "A curved row of dominoes.", dominoes, { 2, 1, 4 }, 24, 0.4f, 0.7f },
        { "Flat layers", "Planar bodies in two planes, with 3D balls rolling through the front plane.", flatLayers,
          { 1, 2.5f, 0 }, 18, 0.35f, 0.25f },
        { "Sweeper", "A kinematic arm sweeping loose bodies.", sweeper, { 0, 1, 0 }, 22, 0.5f, 0.6f },
        { "Slopes", "Friction 0.02, 0.25 and 0.9. Restitution 0, 0.5 and 0.9.", slopes, { 3, 2, 0 }, 20, 0.9f, 0.3f },
        { "Sleeping", "Bodies created asleep until something touches them.", dormant, { 0, 3, 0 }, 18, 0.6f, 0.3f },
        { "Hulls", "Rocks, gems and cones as convex hulls of point clouds, falling onto a static rock.", hulls, { 0, 2, 0 }, 20,
          0.7f, 0.4f },
        { "Cutouts", "Traced 2D layers with holes and concave outlines, catching on pegs and piling into a traced bowl.",
          cutouts, { 0, 4.5f, 0 }, 17, 0.25f, 0.15f },
    };
    return scenes;
}
