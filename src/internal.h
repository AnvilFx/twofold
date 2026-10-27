#pragma once

#include <cstdint>
#include <vector>

#include "twofold/world.h"

namespace twofold {

struct Polyhedron {
    struct Face {
        Vec3 normal;
        float offset = 0.0f;
        int first = 0;
        int count = 0;
    };
    struct Edge {
        int v0 = 0, v1 = 0;
        int f0 = 0, f1 = 0;
    };

    std::vector<Vec3> vertices;
    std::vector<int> indices;
    std::vector<Face> faces;
    std::vector<Edge> edges;
    Vec3 centroid;
};

// Faces are vertex loops, counter-clockwise seen from outside.
Polyhedron makePolyhedron(std::vector<Vec3> vertices, const std::vector<std::vector<int>>& faces);
Polyhedron makePrism(const std::vector<Vec2>& polygon, float zMin, float zMax);
Polyhedron makeBox(const Vec3& halfExtents);
MassProperties polyhedronMass(const Polyhedron& poly, float density);

enum class ShapeKind : uint8_t { Sphere, Convex, HalfSpace };

struct Shape {
    uint32_t generation = 0;
    bool alive = false;
    uint32_t body = 0;
    Geometry geometry;
    Transform local;
    float density = 1.0f;
    float friction = 0.5f;
    float restitution = 0.0f;
    uint32_t category = 1;
    uint32_t mask = 0xffffffffu;
    void* userData = nullptr;

    ShapeKind kind = ShapeKind::Convex;
    float radius = 0.0f;
    Polyhedron poly;
    // Distance from the shape origin to its farthest point.
    float extent = 0.0f;
    Vec3 centroid;
    // Distance from the centroid to the nearest face.
    float minExtent = 0.0f;
};

void buildShape(Shape& shape);
MassProperties shapeMass(const Shape& shape);
void shapeBounds(const Shape& shape, const Transform& xf, Vec3& lo, Vec3& hi);

// Points in shape space, grown by radius. A half-space has no points.
struct Convex {
    const Vec3* points = nullptr;
    int count = 0;
    float radius = 0.0f;
    bool halfSpace = false;
};

Convex convexOf(const Shape& shape);

// Normal points from A to B. Overlapping shapes have a distance of zero or less.
struct DistanceResult {
    Vec3 pointA;
    Vec3 pointB;
    Vec3 normal;
    float distance = 0.0f;
};

DistanceResult distance(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb);

// Center of mass and rotation of a body moving from its pose at t = 0 to its pose at t = 1.
struct Sweep {
    Vec3 localCenter;
    Vec3 c0, c1;
    Quat q0;
    Vec3 axis;
    float angle = 0.0f;

    Transform at(float t) const;
};

Sweep makeSweep(const Vec3& localCenter, const Vec3& c0, const Quat& q0, const Vec3& c1, const Quat& q1);

enum class Impact : uint8_t { Miss, Hit, Touching };

// First t at which two swept shapes come within target of each other. Reach bounds the distance from
// the center of mass to any point of the shape. Shapes already that close at t = 0 are Touching.
Impact timeOfImpact(const Convex& a, const Transform& localA, const Sweep& sweepA, float reachA, const Convex& b,
                    const Transform& localB, const Sweep& sweepB, float reachB, float target, float& t);

// Casts origin + t * d for t in [0, maxFraction]. A ray that starts inside the shape misses.
bool rayShape(const Shape& shape, const Transform& xf, const Vec3& origin, const Vec3& d, float maxFraction,
              float& fraction, Vec3& normal);

struct Body {
    uint32_t generation = 0;
    bool alive = false;
    MotionType type = MotionType::Dynamic;

    Transform xf;
    Vec3 com;
    Vec3 localCenter;
    Vec3 v;
    Vec3 w;
    Vec3 force;
    Vec3 torque;

    float mass = 0.0f;
    float invMass = 0.0f;
    Mat3 inertia = Mat3::zero();
    bool massOverride = false;
    MassProperties overrideMass;

    float linearDamping = 0.0f;
    float angularDamping = 0.0f;
    float gravityScale = 1.0f;
    uint8_t locks = 0;
    Quat lockFrame;

    bool awake = true;
    bool canSleep = true;
    float sleepTime = 0.0f;
    uint32_t island = UINT32_MAX;

    bool hasTarget = false;
    Transform target;

    std::vector<uint32_t> shapes;
    void* userData = nullptr;

    Mat3 invMassWorld = Mat3::zero();
    Mat3 invInertiaWorld = Mat3::zero();
    Mat3 linearProj;
    Mat3 angularProj;
    // Motion since qBase, the orientation contact anchors are measured in.
    Vec3 dp;
    Quat dq;
    Quat qBase;
    Quat q0;
    Vec3 com0;
    // How far any point of the body could move at the velocities the step started with.
    float expectedTravel = 0.0f;

    bool planar() const { return (locks & LockPlanar) == LockPlanar; }
    Vec3 planeNormal() const { return lockFrame.rotate({ 0, 0, 1 }); }
};

struct ManifoldPoint {
    Vec3 position;
    float separation = 0.0f;
    float predicted = 0.0f;
};

// Normal points from shape A to shape B.
struct Manifold {
    Vec3 normal;
    ManifoldPoint points[4];
    int count = 0;
};

struct Tolerances {
    float linearSlop = 0.005f;
    float speculative = 0.02f;
};

// Velocities of the two bodies over the coming step.
struct Motion {
    Vec3 linearA;
    Vec3 angularA;
    Vec3 centerA;
    Vec3 linearB;
    Vec3 angularB;
    Vec3 centerB;
    float dt = 0.0f;

    // Separation at p after the step, with n pointing from A to B.
    float predict(const Vec3& p, const Vec3& n, float separation) const {
        const Vec3 va = linearA + cross(angularA, p - centerA);
        const Vec3 vb = linearB + cross(angularB, p - centerB);
        return separation + dt * dot(n, vb - va);
    }

    Motion swapped() const { return { linearB, angularB, centerB, linearA, angularA, centerA, dt }; }
};

void collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
             float margin, const Tolerances& tol, const Motion& motion, Manifold& m);

// Collides two shapes that both extrude along n, as 2D profiles in the plane normal to n.
// Returns false when either shape does not extrude along n.
bool collidePlanar(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
                   const Vec3& n, float margin, const Tolerances& tol, Manifold& m);

void reduceManifold(std::vector<ManifoldPoint>& points, const Vec3& normal, float band, Manifold& m);

struct ContactPoint {
    Vec3 localA;
    Vec3 anchorA;
    Vec3 anchorB;
    float separation = 0.0f;
    float adjustedSeparation = 0.0f;
    float normalImpulse = 0.0f;
    Vec3 tangentImpulse;
    float maxNormalImpulse = 0.0f;
    float normalMass = 0.0f;
    float tangentMass1 = 0.0f;
    float tangentMass2 = 0.0f;
    float relativeVelocity = 0.0f;
};

struct Contact {
    uint64_t key = 0;
    uint32_t shapeA = 0, shapeB = 0;
    uint32_t bodyA = 0, bodyB = 0;
    Vec3 normal;
    Vec3 tangent1, tangent2;
    float friction = 0.0f;
    float restitution = 0.0f;
    int count = 0;
    ContactPoint points[4];

    bool touching() const {
        for (int i = 0; i < count; ++i) {
            if (points[i].separation < 0.0f || points[i].normalImpulse > 0.0f) { return true; }
        }
        return false;
    }
};

struct Softness {
    float biasRate = 0.0f;
    float massScale = 1.0f;
    float impulseScale = 0.0f;
};

Softness makeSoft(float hertz, float zeta, float h);

struct SolverContext {
    float h = 0.0f;
    float invH = 0.0f;
    Softness contactSoft;
    Softness staticSoft;
    float maxPushout = 0.0f;
    float restitutionThreshold = 0.0f;
};

void prepareContact(const Body& a, const Body& b, Contact& c);
void prepareContacts(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active);
// Re-expresses anchors as the bodies' current lever arms, before their dp and dq are reset.
void rebaseContact(const Body& a, const Body& b, Contact& c);
void warmStart(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active);
void solveContacts(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active,
                   const SolverContext& ctx, bool useBias);
void applyRestitution(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active,
                      const SolverContext& ctx);

}
