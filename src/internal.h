#pragma once

#include <cstdint>
#include <vector>

#include "tree.h"
#include "twofold/mesh.h"
#include "twofold/world.h"

#if !defined(TWOFOLD_NO_SIMD) && (defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#define TWOFOLD_SSE2 1
#elif !defined(TWOFOLD_NO_SIMD) && (defined(__aarch64__) || defined(_M_ARM64))
#define TWOFOLD_NEON 1
#endif

namespace twofold {

// Contacts solved side by side.
#if defined(TWOFOLD_SSE2) || defined(TWOFOLD_NEON)
constexpr int kLanes = 4;
#else
constexpr int kLanes = 1;
#endif
constexpr int kColors = 32;

// std::fmax and std::fmin, NaN and signed zeros included, but inline. MSVC makes them library calls.
// Ordered unequal values take one select. Only equal values and NaN take the branch.
inline float maxf(float a, float b) {
    if (a == b || a != a || b != b) {
        if (a != a) { return b; }
        if (b != b) { return a; }
        // Equal. Only zeros can differ, and the sum of two zeros is negative only when both are.
        return a == 0.0f ? a + b : a;
    }
    return a > b ? a : b;
}

inline float minf(float a, float b) {
    if (a == b || a != a || b != b) {
        if (a != a) { return b; }
        if (b != b) { return a; }
        return a == 0.0f ? -(-a + -b) : a;
    }
    return a < b ? a : b;
}

struct TriangleMesh::Data {
    // A leaf holds count triangles from first. An inner node has a count of zero and its children at first
    // and first + 1.
    struct Node {
        Vec3 lo, hi;
        int first = 0;
        int count = 0;
    };

    std::vector<Vec3> vertices;
    std::vector<uint32_t> indices;
    std::vector<Vec3> normals;
    // Bit k is set when edge k, from corner k to corner k + 1, may give a contact its own normal on the
    // front side. Bit 3 + k is the same for the back side.
    std::vector<uint8_t> edges;
    std::vector<Node> nodes;

    void corners(int t, Vec3 out[3]) const {
        for (int k = 0; k < 3; ++k) { out[k] = vertices[indices[3 * t + k]]; }
    }

    template <class Visit>
    void query(const Vec3& lo, const Vec3& hi, Visit&& visit) const {
        if (nodes.empty()) { return; }
        int stack[64];
        int top = 0;
        stack[top++] = 0;
        while (top > 0) {
            const Node& n = nodes[stack[--top]];
            if (n.lo.x > hi.x || n.hi.x < lo.x || n.lo.y > hi.y || n.hi.y < lo.y || n.lo.z > hi.z || n.hi.z < lo.z) {
                continue;
            }
            if (n.count > 0) {
                for (int t = n.first; t < n.first + n.count; ++t) { visit(t); }
            } else {
                stack[top++] = n.first + 1;
                stack[top++] = n.first;
            }
        }
    }
};

struct Polyhedron {
    // A face's vertex loop in indices.
    struct Face {
        int first = 0;
        int count = 0;
    };
    struct Edge {
        int v0 = 0, v1 = 0;
        int f0 = 0, f1 = 0;
        float length = 0.0f;
    };

    std::vector<Vec3> vertices;
    std::vector<int> indices;
    std::vector<Face> faces;
    std::vector<Vec3> normals;
    std::vector<float> offsets;
    std::vector<Edge> edges;
    // Cross of each edge's two face normals: the edge's arc on the Gauss map.
    std::vector<Vec3> arcs;
    Vec3 centroid;
    Vec3 lo, hi;
};

// Faces are vertex loops, counter-clockwise seen from outside.
Polyhedron makePolyhedron(std::vector<Vec3> vertices, const std::vector<std::vector<int>>& faces);
Polyhedron makePrism(const std::vector<Vec2>& polygon, float zMin, float zMax);
Polyhedron makeBox(const Vec3& halfExtents);
MassProperties polyhedronMass(const Polyhedron& poly, float density);

enum class ShapeKind : uint8_t { Sphere, Convex, HalfSpace, Mesh, Cylinder };

enum class ProfileKind : uint8_t { None, Polygon, Circle, HalfPlane };

// A shape seen along a plane normal, in its body's frame, with 2D coordinates along u and v and depth along
// key, the body-frame normal it was made for. None when the shape doesn't extrude along key.
struct PlanarProfile {
    Vec3 key, u, v;
    ProfileKind kind = ProfileKind::None;
    std::vector<Vec2> points, normals;
    Vec2 center;
    float radius = 0.0f;
    // Half the size of the box around the polygon or circle, centered on center.
    Vec2 extent;
    // The half-plane holds the points p at depth z with dot(normal, p) <= offset - tilt * z.
    Vec2 normal;
    float offset = 0.0f, tilt = 0.0f;
    float zLo = 0.0f, zHi = 0.0f;
};

// A body's pose in a plane: the plane normal it was placed for, that normal in the body's frame with the
// body's 2D axes around it, and the body's position, turn and depth in the plane's own 2D frame.
struct PlanarPose {
    Vec3 normal;
    Vec3 key, u, v;
    Vec3 axisU, axisV;
    Vec2 position, turn;
    float depth = 0.0f;
    uint32_t stamp = 0;
};

void placeInPlane(PlanarPose& pose, const Transform& xf, const Vec3& n);

struct Shape {
    uint32_t generation = 0;
    bool alive = false;
    uint32_t body = 0;
    Geometry geometry;
    Transform local;
    float density = 1.0f;
    float friction = 0.5f;
    float restitution = 0.0f;
    float rollingResistance = 0.0f;
    uint32_t category = 1;
    uint32_t mask = 0xffffffffu;
    bool sensor = false;
    void* userData = nullptr;

    ShapeKind kind = ShapeKind::Convex;
    // A cylinder runs along local z.
    float radius = 0.0f;
    float halfHeight = 0.0f;
    Polyhedron poly;
    // Distance from the shape origin to its farthest point.
    float extent = 0.0f;
    Vec3 centroid;
    // Distance from the centroid to the nearest face.
    float minExtent = 0.0f;

    // Bounds at the body's current pose, and the shape's leaf in its body's own tree when it has one.
    Vec3 lo, hi;
    int part = -1;
    // Made on first use by the planar path, for its body's pose key, and again when that key changes.
    mutable PlanarProfile profile;
    // Distance from the centroid to the nearest face along the plane a planar body moves in, for that
    // plane's normal in the shape's frame. Negative until found.
    mutable Vec3 innerKey;
    mutable float inner = -1.0f;
};

void buildShape(Shape& shape);
MassProperties shapeMass(const Shape& shape);
void shapeBounds(const Shape& shape, const Transform& xf, Vec3& lo, Vec3& hi);
// Half the size of the bounds of a cylinder along each axis.
Vec3 cylinderReach(const Vec3& axis, float radius, float halfHeight);
// The point of a cylinder along local z farthest along d.
Vec3 cylinderSupport(float radius, float halfHeight, const Vec3& d);
// Null unless the shape is a mesh with at least one triangle.
const TriangleMesh::Data* meshOf(const Shape& shape);
// Bounds in the mesh's frame of the box lo, hi in the world.
void boundsInMesh(const Transform& xm, const Vec3& lo, const Vec3& hi, Vec3& meshLo, Vec3& meshHi);

// Points in shape space, grown by radius. A half-space has no points, and neither has a cylinder along z
// with a rim radius.
struct Convex {
    const Vec3* points = nullptr;
    int count = 0;
    float radius = 0.0f;
    bool halfSpace = false;
    bool cylinder = false;
    float rim = 0.0f;
    float halfHeight = 0.0f;
};

Convex convexOf(const Shape& shape);

// Normal points from A to B. Overlapping shapes have a distance of zero or less.
struct DistanceResult {
    Vec3 pointA;
    Vec3 pointB;
    Vec3 normal;
    float distance = 0.0f;
};

// A nonzero start is a guess at the normal, which only changes how fast the search converges.
DistanceResult distance(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb,
                        const Vec3& start = {});
// Whether distance(a, xa, b, xb).distance > limit, decided as soon as the search can tell.
bool apart(const Convex& a, const Transform& xa, const Convex& b, const Transform& xb, float limit);

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

// Casts o + t * d in the mesh's frame. The normal faces the ray.
bool rayMesh(const TriangleMesh::Data& mesh, const Vec3& o, const Vec3& d, float maxFraction, float& fraction,
             Vec3& normal);

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
    // How far the shapes reach from the center of mass.
    float extent = 0.0f;
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
    // Bounds of the shapes other than half-spaces, and the body's leaf in the world's tree.
    Vec3 lo, hi;
    int proxy = -1;
    // A body with many shapes keeps their leaves in a tree of its own, so finding the ones near a box doesn't
    // test them all.
    BoundsTree parts;
    bool hasParts = false;
    // Placed by the planar path at most once per round of collisions.
    mutable PlanarPose planePose;

    Mat3 invMassWorld = Mat3::zero();
    Mat3 invInertiaWorld = Mat3::zero();
    // For a body locked to a plane, its inverse mass and its inverse inertia about the plane normal. Zero while
    // it can't move.
    float planeInvMass = 0.0f;
    float planeInvInertia = 0.0f;
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

// The faces or edges of two polyhedra that a manifold was clipped from.
struct Feature {
    enum Kind : uint8_t { None, FaceA, FaceB, Edges };
    Kind kind = None;
    int a = -1, b = -1;
};

// Normal points from shape A to shape B.
struct Manifold {
    Vec3 normal;
    ManifoldPoint points[4];
    int count = 0;
    Feature feature;
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
    // The pair travels further this step than the speculative distance.
    bool fast = false;

    // Separation at p after the step, with n pointing from A to B.
    float predict(const Vec3& p, const Vec3& n, float separation) const {
        const Vec3 va = linearA + cross(angularA, p - centerA);
        const Vec3 vb = linearB + cross(angularB, p - centerB);
        return separation + dt * dot(n, vb - va);
    }

    Motion swapped() const { return { linearB, angularB, centerB, linearA, angularA, centerA, dt, fast }; }
};

// With reuse, two polyhedra test the edge pair it names, while its arcs still cross, in place of every edge pair.
// Two boxes test every axis in closed form and need no reuse.
void collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb,
             float margin, const Tolerances& tol, const Motion& motion, Manifold& m, const Feature* reuse = nullptr);

// Collides two shapes that both extrude along n, as 2D profiles in the plane normal to n. Takes the poses of
// the shapes' bodies. Returns false when either shape does not extrude along n.
bool collidePlanar(const Shape& a, const PlanarPose& bodyA, const Shape& b, const PlanarPose& bodyB, float margin,
                   const Tolerances& tol, Manifold& m);

// A mesh triangle found no nearer to the other shape than distance.
struct TriangleBound {
    int triangle = 0;
    float distance = 0.0f;
};

// Triangles ruled out by distance are added to found. Those in known, sorted by triangle, are skipped while
// their distance less moved, how far the shapes have moved since, stays past the margin.
struct TriangleBounds {
    std::vector<TriangleBound>* found = nullptr;
    const TriangleBound* known = nullptr;
    int knownCount = 0;
    float moved = 0.0f;
};

// Collides a sphere or polyhedron with a mesh. Each group of triangles facing the same way gives one
// manifold, with its normal pointing from the mesh.
void collideMesh(const Shape& mesh, const Transform& xm, const Shape& other, const Transform& xo, float margin,
                 const Tolerances& tol, const Motion& motion, std::vector<Manifold>& out,
                 TriangleBounds* bounds = nullptr);

// The first point kept is the first one within tie of the deepest.
void reduceManifold(std::vector<ManifoldPoint>& points, const Vec3& normal, float band, Manifold& m, float tie = 0.0f,
                    bool fast = false);

// Whether origin + t * d for some t in [0, maxFraction] lies in the box lo, hi.
bool segmentHitsBox(const Vec3& lo, const Vec3& hi, const Vec3& origin, const Vec3& d, float maxFraction);

// Whether the point is within tolerance of the shape. A mesh only holds points near its triangles.
bool shapeContains(const Shape& shape, const Transform& xf, const Vec3& point, float tolerance);

// Whether the convex shape q at xq comes within tolerance of the shape.
bool shapesOverlap(const Shape& shape, const Transform& xf, const Shape& q, const Transform& xq, float tolerance);

// The point of the core of c farthest along d. A half-space has none.
Vec3 convexSupport(const Convex& c, const Transform& xf, const Vec3& d);

// Sweeps the convex shape q from xq along d without turning, as far as fraction, and leaves the fraction where
// it first comes within touch of the shape there. Point and normal are on the shape, the normal facing q. A start
// within touch is Touching, with the deepest point of the overlap and the normal that pushes q out.
Impact castShape(const Shape& shape, const Transform& xf, const Shape& q, const Transform& xq, const Vec3& d,
                 float touch, const Tolerances& tol, float& fraction, Vec3& point, Vec3& normal);

struct ContactPoint {
    // Where the point touches each shape, in its body's frame.
    Vec3 surfaceA, surfaceB;
    Vec3 localA;
    Vec3 anchorA;
    Vec3 anchorB;
    float separation = 0.0f;
    float adjustedSeparation = 0.0f;
    float normalImpulse = 0.0f;
    Vec3 tangentImpulse;
    float maxNormalImpulse = 0.0f;
    float relativeVelocity = 0.0f;
};

// One direction at a contact point: its angular Jacobians and the velocity change a unit impulse
// along it gives each body.
struct Axis {
    Vec3 armA, armB;
    Vec3 linearA, angularA;
    Vec3 linearB, angularB;
    float mass = 0.0f;
};

// Built once per step. Contact ci uses entries 4 * ci to 4 * ci + 3.
struct PointAxes {
    Axis normal, tangent1, tangent2;
};

// The axes of a point of a planar contact: for the normal and the tangent in the plane, the turn about the
// plane normal a unit impulse gives each body per unit of its inverse inertia, and the effective mass.
struct PlanarAxes {
    float normalA = 0.0f, normalB = 0.0f, normalMass = 0.0f;
    float tangentA = 0.0f, tangentB = 0.0f, tangentMass = 0.0f;
};

struct Contact {
    uint64_t key = 0;
    // Orders the contacts of one pair. Only a mesh pair has more than one.
    uint32_t part = 0;
    uint32_t shapeA = 0, shapeB = 0;
    uint32_t bodyA = 0, bodyB = 0;
    Vec3 normal;
    Vec3 tangent1, tangent2;
    float friction = 0.0f;
    float restitution = 0.0f;
    float rollingResistance = 0.0f;
    // Mean height of the turning body's center of mass above the contact points.
    float rollingRadius = 0.0f;
    // Angular impulse on body B.
    Vec3 rollingImpulse;
    // Solved in the plane with normal plane, since every body it moves is locked to that plane and the normal
    // lies in it. Its tangent1 is then the tangent in the plane.
    bool planar = false;
    Vec3 plane;
    int count = 0;
    ContactPoint points[4];
    // The normal and the pose of body B in body A's frame when the points were collided.
    Vec3 localNormal;
    Vec3 baseOffset;
    Quat baseTurn;
    Feature feature;
    bool reusable = false;
    // A contact removed between steps stays in place until the next step drops it.
    bool removed = false;
    // The pair was touching after the last step, as contact events report it. Every contact of a pair agrees.
    bool touched = false;
    // Solved in a lane batch this step, which holds its axes in place of SolverData::axes.
    bool laned = false;

    bool touching() const {
        for (int i = 0; i < count; ++i) {
            if (points[i].separation < 0.0f || points[i].normalImpulse > 0.0f) { return true; }
        }
        return false;
    }
};

// Joint ends that use this index are fixed in the world.
constexpr uint32_t kWorldBody = UINT32_MAX;

struct Joint {
    uint32_t generation = 0;
    bool alive = false;
    JointType type = JointType::Point;
    uint32_t bodyA = kWorldBody, bodyB = 0;
    bool collideConnected = false;
    void* userData = nullptr;

    // Anchors, axis and a reference perpendicular to the axis, in each body's frame.
    Vec3 localA, localB;
    Vec3 axisA, axisB;
    Vec3 refA, refB;
    // Body B's rotation in body A's frame when joined.
    Quat restTurn;
    float length = 0.0f;
    bool limited = false;
    float lower = 0.0f, upper = 0.0f;
    float motorSpeed = 0.0f, maxMotor = 0.0f;
    float hertz = 0.0f, dampingRatio = 0.0f;
    float maxFriction = 0.0f;
    // A swing-twist joint's frame in each body, x along the axis and y along the normal, and the sines of half
    // its swing limits.
    Quat frameA, frameB;
    float swingY = 1.0f, swingZ = 1.0f;
    float breakForce = 0.0f, breakTorque = 0.0f;
    bool broken = false;

    Vec3 linearImpulse;
    Vec3 angularImpulse;
    float axialImpulse = 0.0f;
    float motorImpulse = 0.0f;
    float lowerImpulse = 0.0f, upperImpulse = 0.0f;
    float frictionImpulse = 0.0f;
    float swingImpulse = 0.0f;
    Vec3 frictionTurn;
    // Force and torque about anchor B per unit impulse of the axial or twist row and of the swing row, as last
    // solved.
    Vec3 rowForce[2], rowTorque[2];
    // On body B, summed over the substeps of the last step the joint was active, and then divided by its length.
    Vec3 force, torque;
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
    Softness jointSoft;
    float maxPushout = 0.0f;
    float restitutionThreshold = 0.0f;
};

struct SolverData {
    std::vector<Body>& bodies;
    std::vector<Contact>& contacts;
    std::vector<PointAxes>& axes;
    std::vector<PlanarAxes>& planarAxes;
    std::vector<Joint>& joints;
    // Stands in for the world at joint ends fixed in it.
    Body& world;
};

// One direction at each lane's point.
struct LaneAxis {
    float armA[3][kLanes], armB[3][kLanes];
    float linearA[3][kLanes], angularA[3][kLanes];
    float linearB[3][kLanes], angularB[3][kLanes];
    float mass[kLanes];
};

// Up to kLanes contacts that move no body in common, solved side by side. Masks have every bit set or none.
struct ContactLanes {
    // How many lanes hold a contact. The rest repeat the first contact and move nothing.
    int size = 0;
    // The most points any contact has, and whether every contact has that many.
    int points = 0;
    bool even = true;
    uint32_t contact[kLanes];
    uint32_t bodyA[kLanes], bodyB[kLanes];
    uint32_t moveA[kLanes], moveB[kLanes];
    uint32_t count[kLanes];
    uint32_t used[4][kLanes];
    float biasRate[kLanes], massScale[kLanes], impulseScale[kLanes];
    float normal[3][kLanes], tangent1[3][kLanes], tangent2[3][kLanes];
    float friction[kLanes];
    float anchorA[4][3][kLanes], anchorB[4][3][kLanes];
    float adjustedSeparation[4][kLanes];
    LaneAxis normalAxis[4];
    LaneAxis tangentAxis[4][2];
};

// Contacts grouped by color: no two contacts of one color move the same body, so a color's contacts can be
// solved in any order, and kLanes at a time. Colors are solved in order, then the rest one by one. Within a
// color, contacts with more points come first so that lanes in a batch tend to have as many points.
struct ContactBatches {
    std::vector<ContactLanes> lanes;
    std::vector<uint32_t> rest;
    std::vector<uint32_t> byColor;
};

struct ContactColors {
    // By body index, a bit for each color one of its contacts has.
    std::vector<uint32_t> bodies;
    // By contact index. kColors is no color.
    std::vector<uint8_t> of;
};

Axis makeAxis(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& d);
float speedAlong(const Body& a, const Body& b, const Vec3& d, const Axis& x);
void push(Body& a, Body& b, const Axis& x, float impulse);

void prepareContact(const Body& a, const Body& b, Contact& c, PlanarAxes* planar);
void prepareContacts(SolverData& s, const std::vector<uint32_t>& active);
// The point axes of a contact that is solved on its own, once it is known not to be laned.
void makeContactAxes(const Body& a, const Body& b, const Contact& c, PointAxes* axes);
// Re-expresses anchors as the bodies' current lever arms, before their dp and dq are reset.
void rebaseContact(const Body& a, const Body& b, Contact& c, PointAxes* axes, PlanarAxes* planar);
void warmStart(SolverData& s, const std::vector<uint32_t>& active);
void solveContacts(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx, bool useBias);
// Colors the contacts and marks those that will be laned. A planar or rolling contact gets no color.
void colorContacts(SolverData& s, const std::vector<uint32_t>& list, ContactColors& colors);
// Groups the listed contacts by their colors, then fills the lanes as prepareLanes does.
void batchContacts(SolverData& s, const std::vector<uint32_t>& list, const ContactColors& colors, const SolverContext& ctx,
                   ContactBatches& out);
// Takes the contacts' current points into their lanes and works out the axes there.
void prepareLanes(SolverData& s, const SolverContext& ctx, ContactBatches& batches);
void warmStartBatches(SolverData& s, const ContactBatches& batches);
void solveBatches(SolverData& s, const ContactBatches& batches, const SolverContext& ctx, bool useBias);
void restituteBatches(SolverData& s, const ContactBatches& batches, const SolverContext& ctx);
void applyRestitution(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx);

void warmStartJoints(SolverData& s, const std::vector<uint32_t>& active);
// The relax pass adds each joint's impulses on body B to its force and torque.
void solveJoints(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx, bool useBias);
// Moves joined bodies back together without changing their velocities.
void settleJoints(SolverData& s, const std::vector<uint32_t>& active, int iterations);
// The hinge angle, slider translation, distance or twist.
float jointValue(const Joint& j, const Body& a, const Body& b);
// Adds the impulse and the angular impulse about anchor B the joint gave body B in the substep just solved.
void addReaction(Joint& j);

}
