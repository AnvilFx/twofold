#pragma once

#include <cstdint>
#include <memory>
#include <variant>
#include <vector>

#include "math.h"

namespace twofold {

struct BodyId {
    uint32_t index = UINT32_MAX;
    uint32_t generation = 0;

    bool operator==(const BodyId& b) const { return index == b.index && generation == b.generation; }
    bool operator!=(const BodyId& b) const { return !(*this == b); }
};

struct ShapeId {
    uint32_t index = UINT32_MAX;
    uint32_t generation = 0;

    bool operator==(const ShapeId& b) const { return index == b.index && generation == b.generation; }
    bool operator!=(const ShapeId& b) const { return !(*this == b); }
};

enum class MotionType : uint8_t { Static, Kinematic, Dynamic };

// Degrees of freedom a dynamic body may not use, measured along the axes of BodyDesc::lockFrame.
enum LockFlags : uint8_t {
    LockTranslationX = 1 << 0,
    LockTranslationY = 1 << 1,
    LockTranslationZ = 1 << 2,
    LockRotationX = 1 << 3,
    LockRotationY = 1 << 4,
    LockRotationZ = 1 << 5,
    // Moves in the lock frame's xy plane and turns about its z axis.
    LockPlanar = LockTranslationZ | LockRotationX | LockRotationY,
};

struct Sphere {
    float radius = 0.5f;
};

struct Box {
    Vec3 halfExtents{ 0.5f, 0.5f, 0.5f };
};

// Axis along local z.
struct Cylinder {
    float radius = 0.5f;
    float halfHeight = 0.5f;
};

// A convex polygon in the local xy plane, extruded along z. Either winding is accepted.
struct Prism {
    std::vector<Vec2> polygon;
    float zMin = -0.5f;
    float zMax = 0.5f;
};

// The convex hull of a point cloud, as computeHull builds it.
struct Hull {
    std::vector<Vec3> points;
    int maxVertices = 64;
};

// Everything below the local xy plane. Static bodies only.
struct HalfSpace {};

using Geometry = std::variant<Sphere, Box, Cylinder, Prism, Hull, HalfSpace>;

struct ShapeDesc {
    Geometry geometry = Box{};
    Vec3 position;
    Quat rotation;
    float density = 1.0f;
    float friction = 0.5f;
    float restitution = 0.0f;
    // Two shapes collide when each one's category overlaps the other's mask.
    uint32_t category = 1;
    uint32_t mask = 0xffffffffu;
    void* userData = nullptr;
};

struct BodyDesc {
    MotionType type = MotionType::Dynamic;
    Vec3 position;
    Quat rotation;
    Vec3 linearVelocity;
    Vec3 angularVelocity;
    float linearDamping = 0.0f;
    float angularDamping = 0.0f;
    float gravityScale = 1.0f;
    uint8_t locks = 0;
    Quat lockFrame;
    // A body created asleep stays put until something touches it or it is woken.
    bool awake = true;
    bool canSleep = true;
    void* userData = nullptr;
};

// Center and inertia are in the body's local frame, inertia about the center.
struct MassProperties {
    float mass = 0.0f;
    Vec3 center;
    Mat3 inertia = Mat3::zero();
};

struct WorldSettings {
    Vec3 gravity{ 0.0f, -9.81f, 0.0f };
    // World units per meter. Every tolerance and speed below is in meters and scales by this.
    float lengthUnit = 1.0f;
    int substeps = 4;
    float contactHertz = 30.0f;
    float contactDampingRatio = 10.0f;
    float maxPushoutSpeed = 3.0f;
    float restitutionThreshold = 1.0f;
    float sleepThreshold = 0.05f;
    float timeToSleep = 0.5f;
};

struct MaterialInput {
    float value;
    MotionType type;
};

using MixFunction = float (*)(const MaterialInput& a, const MaterialInput& b);

struct ContactInfo {
    BodyId bodyA;
    BodyId bodyB;
    ShapeId shapeA;
    ShapeId shapeB;
    Vec3 normal;
    float friction = 0.0f;
    float restitution = 0.0f;
    int pointCount = 0;
    Vec3 points[4];
    float separations[4] = {};
    float normalImpulses[4] = {};
    Vec3 tangentImpulses[4];
};

struct Aabb {
    Vec3 lower;
    Vec3 upper;
};

struct RayHit {
    BodyId body;
    ShapeId shape;
    Vec3 point;
    Vec3 normal;
    // Distance to the hit as a fraction of the ray translation.
    float fraction = 1.0f;
};

// Timings of the last step in milliseconds, and what it processed.
struct StepProfile {
    float total = 0.0f;
    float broadphase = 0.0f;
    float narrowphase = 0.0f;
    float solver = 0.0f;
    float sleep = 0.0f;
    int candidatePairs = 0;
    int manifolds = 0;
    int touchingContacts = 0;
    int awakeBodies = 0;
    int awakeIslands = 0;
};

class World {
public:
    explicit World(const WorldSettings& settings = {});
    ~World();
    // A copy steps exactly like the original, so it can serve as a snapshot.
    World(const World& other);
    World& operator=(const World& other);
    World(World&&) noexcept;
    World& operator=(World&&) noexcept;

    const WorldSettings& settings() const;
    void setSettings(const WorldSettings& settings);
    void setGravity(const Vec3& gravity);
    void setFrictionMix(MixFunction mix);
    void setRestitutionMix(MixFunction mix);

    void step(float dt);
    const StepProfile& profile() const;

    BodyId createBody(const BodyDesc& desc);
    void destroyBody(BodyId body);
    bool isValid(BodyId body) const;
    std::vector<BodyId> bodies() const;

    ShapeId addShape(BodyId body, const ShapeDesc& desc);
    void removeShape(ShapeId shape);
    bool isValid(ShapeId shape) const;
    BodyId shapeBody(ShapeId shape) const;
    std::vector<ShapeId> shapes(BodyId body) const;
    Geometry geometry(ShapeId shape) const;
    Transform localTransform(ShapeId shape) const;
    float friction(ShapeId shape) const;
    float restitution(ShapeId shape) const;
    Aabb bounds(ShapeId shape) const;
    void setGeometry(ShapeId shape, const Geometry& geometry);
    void setFriction(ShapeId shape, float friction);
    void setRestitution(ShapeId shape, float restitution);
    void setFilter(ShapeId shape, uint32_t category, uint32_t mask);

    Transform transform(BodyId body) const;
    Vec3 position(BodyId body) const;
    Quat rotation(BodyId body) const;
    Vec3 centerOfMass(BodyId body) const;
    void setTransform(BodyId body, const Vec3& position, const Quat& rotation);

    Vec3 linearVelocity(BodyId body) const;
    Vec3 angularVelocity(BodyId body) const;
    void setLinearVelocity(BodyId body, const Vec3& v);
    void setAngularVelocity(BodyId body, const Vec3& w);
    // Sets the velocities that carry a kinematic body to the target pose over the next step of length dt.
    void setKinematicTarget(BodyId body, const Vec3& position, const Quat& rotation, float dt);

    MotionType motionType(BodyId body) const;
    void setMotionType(BodyId body, MotionType type);

    // Forces and torques last for one step. A sleeping body ignores them unless woken.
    void addForce(BodyId body, const Vec3& force, bool wake = true);
    void addForceAtPoint(BodyId body, const Vec3& force, const Vec3& point, bool wake = true);
    void addTorque(BodyId body, const Vec3& torque, bool wake = true);
    void addImpulse(BodyId body, const Vec3& impulse, bool wake = true);
    void addImpulseAtPoint(BodyId body, const Vec3& impulse, const Vec3& point, bool wake = true);

    MassProperties massProperties(BodyId body) const;
    // Replaces the mass computed from the shapes until resetMassProperties is called.
    void setMassProperties(BodyId body, const MassProperties& mass);
    void resetMassProperties(BodyId body);

    uint8_t locks(BodyId body) const;
    Quat lockFrame(BodyId body) const;
    void setLocks(BodyId body, uint8_t locks, const Quat& frame = {});
    float linearDamping(BodyId body) const;
    float angularDamping(BodyId body) const;
    void setDamping(BodyId body, float linear, float angular);
    float gravityScale(BodyId body) const;
    void setGravityScale(BodyId body, float scale);

    bool isAwake(BodyId body) const;
    void wake(BodyId body);
    void sleep(BodyId body);
    bool canSleep(BodyId body) const;
    void setCanSleep(BodyId body, bool canSleep);
    // Bodies in one island sleep and wake together. The island is named by its lowest body index
    // and is updated each step while the body is awake. UINT32_MAX for static and kinematic bodies.
    uint32_t island(BodyId body) const;

    void* userData(BodyId body) const;
    void* userData(ShapeId shape) const;

    // Contacts after the last step. Speculative contacts have points that are not yet touching.
    std::vector<ContactInfo> contacts(bool includeSpeculative = false) const;

    // Closest shape hit by the segment from origin to origin + translation, ignoring shapes whose
    // category does not overlap mask and shapes that contain the origin.
    bool castRay(const Vec3& origin, const Vec3& translation, RayHit& hit, uint32_t mask = 0xffffffffu) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}
