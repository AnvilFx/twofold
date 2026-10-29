#include <algorithm>
#include <cfloat>
#include <chrono>
#include <iterator>

#include "internal.h"
#include "tree.h"

namespace twofold {

namespace {

constexpr float kAligned = 0.9999f;

constexpr int kSettleIterations = 2;

constexpr float kPi = 3.14159265358979f;

// Sine of half a swing limit. Limits are kept above zero so the cone has a direction to push along.
float coneSine(float swing) {
    if (swing >= kPi) { return 1.0f; }
    return std::sin(0.5f * maxf(swing, 1e-3f));
}

// How far past its bounds a shape's leaf in the tree reaches, in meters.
constexpr float kTreeMargin = 0.1f;
// Bodies with more shapes than this keep them in a tree of their own.
constexpr size_t kPartTree = 64;

using Clock = std::chrono::steady_clock;

float millisecondsSince(Clock::time_point start) {
    return std::chrono::duration<float, std::milli>(Clock::now() - start).count();
}

float mixFriction(const MaterialInput& a, const MaterialInput& b) { return std::sqrt(a.value * b.value); }
float mixRestitution(const MaterialInput& a, const MaterialInput& b) { return maxf(a.value, b.value); }

uint64_t pairKey(uint32_t a, uint32_t b) { return ((uint64_t)a << 32) | b; }

Mat3 invertSubspace(const Mat3& m, const Vec3* axes, int k) {
    if (k == 0) { return Mat3::zero(); }
    if (k == 3) { return inverse(m); }
    if (k == 1) {
        const float d = dot(axes[0], m * axes[0]);
        return d > 0.0f ? Mat3::outer(axes[0], axes[0]) * (1.0f / d) : Mat3::zero();
    }
    const float a = dot(axes[0], m * axes[0]);
    const float b = dot(axes[0], m * axes[1]);
    const float c = dot(axes[1], m * axes[0]);
    const float d = dot(axes[1], m * axes[1]);
    float det = a * d - b * c;
    if (det == 0.0f) { return Mat3::zero(); }
    det = 1.0f / det;
    return Mat3::outer(axes[0], axes[0]) * (d * det) - Mat3::outer(axes[0], axes[1]) * (b * det)
         - Mat3::outer(axes[1], axes[0]) * (c * det) + Mat3::outer(axes[1], axes[1]) * (a * det);
}

void updateInverseMass(Body& b) {
    const Mat3 lf = b.lockFrame.toMat3();
    Vec3 ang[3];
    int na = 0;
    b.linearProj = Mat3::zero();
    b.angularProj = Mat3::zero();
    for (int i = 0; i < 3; ++i) {
        const Vec3 axis = lf.column(i);
        if (!(b.locks & (LockTranslationX << i))) { b.linearProj = b.linearProj + Mat3::outer(axis, axis); }
        if (!(b.locks & (LockRotationX << i))) {
            ang[na++] = axis;
            b.angularProj = b.angularProj + Mat3::outer(axis, axis);
        }
    }
    if (b.type != MotionType::Dynamic || !b.awake) {
        b.invMassWorld = Mat3::zero();
        b.invInertiaWorld = Mat3::zero();
        b.planeInvMass = 0.0f;
        b.planeInvInertia = 0.0f;
        return;
    }
    b.invMassWorld = b.linearProj * b.invMass;
    const Mat3 r = b.xf.q.toMat3();
    const Mat3 iw = r * b.inertia * transpose(r);
    b.invInertiaWorld = invertSubspace(iw, ang, na);
    if (b.planar()) {
        const Vec3 n = b.planeNormal();
        b.planeInvMass = b.invMass;
        b.planeInvInertia = dot(n, b.invInertiaWorld * n);
    }
}

bool movingAwake(const Body& b) {
    if (b.type == MotionType::Dynamic) { return b.awake; }
    if (b.type == MotionType::Kinematic) { return lengthSquared(b.v) > 0.0f || lengthSquared(b.w) > 0.0f; }
    return false;
}

// Each new point continues the nearest old point in reach.
// An old point continued by several new points splits its impulse between them.
void carryImpulses(const Contact& old, Contact& c, float matchSq) {
    int from[4], claims[4] = {};
    for (int i = 0; i < c.count; ++i) {
        from[i] = -1;
        float bestSq = matchSq;
        for (int k = 0; k < old.count; ++k) {
            const float d = lengthSquared(old.points[k].localA - c.points[i].localA);
            if (d <= bestSq) {
                from[i] = k;
                bestSq = d;
            }
        }
        if (from[i] >= 0) { ++claims[from[i]]; }
    }
    for (int i = 0; i < c.count; ++i) {
        if (from[i] < 0) { continue; }
        const ContactPoint& o = old.points[from[i]];
        const float share = 1.0f / (float)claims[from[i]];
        c.points[i].normalImpulse = o.normalImpulse * share;
        c.points[i].tangentImpulse = (o.tangentImpulse - c.normal * dot(c.normal, o.tangentImpulse)) * share;
    }
    c.rollingImpulse = old.rollingImpulse;
}

int nearestPoint(const Contact& old, const Vec3& localA) {
    int best = -1;
    float bestSq = FLT_MAX;
    for (int k = 0; k < old.count; ++k) {
        const float d = lengthSquared(old.points[k].localA - localA);
        if (d < bestSq) {
            best = k;
            bestSq = d;
        }
    }
    return best;
}

bool boxesOverlap(const Vec3& loA, const Vec3& hiA, const Vec3& loB, const Vec3& hiB) {
    return loA.x <= hiB.x && loB.x <= hiA.x && loA.y <= hiB.y && loB.y <= hiA.y && loA.z <= hiB.z && loB.z <= hiA.z;
}

}

struct World::Impl {
    WorldSettings settings;
    Tolerances tol;
    MixFunction frictionMix = mixFriction;
    MixFunction restitutionMix = mixRestitution;

    std::vector<Body> bodies;
    std::vector<uint32_t> freeBodies;
    std::vector<Shape> shapes;
    std::vector<uint32_t> freeShapes;
    std::vector<Contact> contacts;
    // Kept between steps so the rebuilt contact list doesn't allocate.
    std::vector<Contact> nextContacts;
    // Lists every step rebuilds, kept for the same reason.
    std::vector<float> growth;
    std::vector<uint64_t> freshPairs, candidatePairs;
    std::vector<uint32_t> activeContacts, fastContacts, fastSolidContacts;
    std::vector<PointAxes> axes;
    std::vector<PlanarAxes> planarAxes;
    // Mesh triangles a fast pair ruled out at the start of the step, so collisions within the step can skip
    // them. Each pair's run is sorted by triangle, and the runs are in key order.
    struct BoundRun {
        uint64_t key;
        uint32_t begin, end;
    };
    std::vector<TriangleBound> triangleBounds;
    std::vector<BoundRun> boundRuns;
    // One leaf per body, over its shapes but half-spaces, whose bounds are infinite and which are kept apart
    // in planes.
    BoundsTree tree;
    std::vector<uint32_t> planes;
    // Sorted keys of body pairs whose leaves overlap, and the bodies whose leaves changed since the pairs
    // were last found.
    std::vector<uint64_t> leafPairs;
    std::vector<uint32_t> refitted;
    std::vector<Joint> joints;
    std::vector<uint32_t> freeJoints;
    std::vector<uint32_t> brokenJoints;
    // Body pairs held by a joint that keeps them from colliding, sorted.
    std::vector<uint64_t> joinedPairs;
    // What links each body to others for waking: links[linkStart[i]] up to links[linkStart[i + 1]]. Built
    // when first needed after a step adds contacts or a joint is made. Contacts and joints removed since stay
    // listed and are skipped.
    struct Link {
        uint64_t key;
        uint32_t other;
        uint32_t joint;
    };
    std::vector<uint32_t> linkStart;
    std::vector<Link> links;
    bool linksValid = false;
    Body worldBody;
    StepProfile profile;
    // Counts rounds of collisions, so planar poses are placed again once bodies have moved.
    uint32_t poseStamp = 1;

    // Sensor shapes in index order, and the sorted keys of sensor and visitor pairs overlapping after the
    // last step, with the sensor in the high half.
    std::vector<uint32_t> sensors;
    std::vector<uint64_t> sensorPairs;
    // A call moved, reshaped or refiltered a shape, so the next step tests still pairs again too.
    bool sensorsStale = false;
    std::vector<ContactEvent> contactEvents;
    std::vector<HitEvent> hitEvents;
    std::vector<SensorEvent> sensorEvents;
    // Events from calls between steps, reported with the next step.
    std::vector<ContactEvent> queuedContactEvents;
    std::vector<SensorEvent> queuedSensorEvents;

    explicit Impl(const WorldSettings& s) {
        worldBody.alive = true;
        worldBody.type = MotionType::Static;
        apply(s);
    }

    void apply(const WorldSettings& s) {
        settings = s;
        tol.linearSlop = 0.005f * s.lengthUnit;
        tol.speculative = 4.0f * tol.linearSlop;
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            if (bodies[bi].alive) { place(bi); }
        }
    }

    // Every leaf covers its body's bounds grown by the speculative distance.
    void fit(uint32_t bi, const Vec3& lo, const Vec3& hi) {
        const Body& b = bodies[bi];
        if (b.proxy < 0) { return; }
        const Vec3 g{ tol.speculative, tol.speculative, tol.speculative };
        if (tree.update(b.proxy, lo - g, hi + g, kTreeMargin * settings.lengthUnit)) { refitted.push_back(bi); }
    }

    // Finds the bounds of a body's shapes, and returns false when it has none but half-spaces.
    bool bound(Body& b) {
        bool any = false;
        const Vec3 g{ tol.speculative, tol.speculative, tol.speculative };
        for (uint32_t si : b.shapes) {
            Shape& s = shapes[si];
            shapeBounds(s, b.xf * s.local, s.lo, s.hi);
            if (s.kind == ShapeKind::HalfSpace) { continue; }
            if (s.part >= 0) { b.parts.update(s.part, s.lo - g, s.hi + g, kTreeMargin * settings.lengthUnit); }
            b.lo = any ? minVec(b.lo, s.lo) : s.lo;
            b.hi = any ? maxVec(b.hi, s.hi) : s.hi;
            any = true;
        }
        return any;
    }

    // A new leaf covers its bounds grown by the speculative distance and the margin.
    int insertLeaf(BoundsTree& t, const Vec3& lo, const Vec3& hi, uint32_t item) const {
        const float g = tol.speculative + kTreeMargin * settings.lengthUnit;
        return t.insert(lo - Vec3{ g, g, g }, hi + Vec3{ g, g, g }, item);
    }

    // Makes the body's leaf and its tree of shapes again after its shapes changed.
    void retrack(uint32_t bi) {
        Body& b = bodies[bi];
        if (b.proxy >= 0) {
            tree.remove(b.proxy);
            b.proxy = -1;
        }
        b.parts = BoundsTree{};
        b.hasParts = false;
        for (uint32_t si : b.shapes) { shapes[si].part = -1; }
        if (!b.alive || !bound(b)) { return; }
        b.proxy = insertLeaf(tree, b.lo, b.hi, bi);
        refitted.push_back(bi);
        if (b.shapes.size() <= kPartTree) { return; }
        b.hasParts = true;
        for (uint32_t si : b.shapes) {
            Shape& s = shapes[si];
            if (s.kind != ShapeKind::HalfSpace) { s.part = insertLeaf(b.parts, s.lo, s.hi, si); }
        }
    }

    void track(uint32_t si) {
        Shape& s = shapes[si];
        if (s.kind == ShapeKind::HalfSpace) { planes.push_back(si); }
        Body& b = bodies[s.body];
        if (!b.hasParts || s.kind == ShapeKind::HalfSpace) {
            retrack(s.body);
            return;
        }
        // A body that already has a tree of shapes takes one more without finding the others again.
        shapeBounds(s, b.xf * s.local, s.lo, s.hi);
        s.part = insertLeaf(b.parts, s.lo, s.hi, si);
        b.lo = minVec(b.lo, s.lo);
        b.hi = maxVec(b.hi, s.hi);
        fit(s.body, b.lo, b.hi);
        refitted.push_back(s.body);
    }

    void untrack(uint32_t si) {
        Shape& s = shapes[si];
        if (s.kind == ShapeKind::HalfSpace) {
            planes.erase(std::remove(planes.begin(), planes.end(), si), planes.end());
        }
        if (s.part >= 0) {
            bodies[s.body].parts.remove(s.part);
            s.part = -1;
        }
    }

    // Updates the bounds of a body's shapes after it moved.
    void place(uint32_t bi) {
        if (bound(bodies[bi])) { fit(bi, bodies[bi].lo, bodies[bi].hi); }
    }

    // Calls visit(shape) for the shapes of a body found in the world's tree. A body with a tree of shapes
    // visits those its search finds, any other body all of them.
    template <class Search, class Visit>
    void forEachShape(uint32_t bi, Search&& search, Visit&& visit) const {
        const Body& b = bodies[bi];
        if (b.hasParts) {
            search(b.parts, visit);
            return;
        }
        for (uint32_t si : b.shapes) {
            if (shapes[si].kind != ShapeKind::HalfSpace) { visit(si); }
        }
    }

    template <class Visit>
    void forEachShape(uint32_t bi, const Vec3& lo, const Vec3& hi, Visit&& visit) const {
        forEachShape(bi, [&](const BoundsTree& t, Visit& v) { t.query(lo, hi, v); }, visit);
    }

    bool pairable(uint32_t a, uint32_t b) const {
        return a != b && (bodies[a].type == MotionType::Dynamic || bodies[b].type == MotionType::Dynamic);
    }
    bool pairable(const Shape& a, const Shape& b) const { return !a.sensor && !b.sensor && pairable(a.body, b.body); }

    // Visits every shape whose leaf overlaps lo, hi, and every half-space. A half-space visits every shape.
    template <class Visit>
    void forEachNear(const Shape& s, const Vec3& lo, const Vec3& hi, Visit&& visit) const {
        if (s.kind == ShapeKind::HalfSpace) {
            for (uint32_t si = 0; si < shapes.size(); ++si) {
                if (shapes[si].alive) { visit(si); }
            }
            return;
        }
        tree.query(lo, hi, [&](uint32_t bi) { forEachShape(bi, lo, hi, visit); });
        for (uint32_t si : planes) { visit(si); }
    }

    Body* get(BodyId id) {
        if (id.index >= bodies.size()) { return nullptr; }
        Body& b = bodies[id.index];
        return (b.alive && b.generation == id.generation) ? &b : nullptr;
    }
    const Body* get(BodyId id) const { return const_cast<Impl*>(this)->get(id); }

    Shape* get(ShapeId id) {
        if (id.index >= shapes.size()) { return nullptr; }
        Shape& s = shapes[id.index];
        return (s.alive && s.generation == id.generation) ? &s : nullptr;
    }
    const Shape* get(ShapeId id) const { return const_cast<Impl*>(this)->get(id); }

    Joint* get(JointId id) {
        if (id.index >= joints.size()) { return nullptr; }
        Joint& j = joints[id.index];
        return (j.alive && j.generation == id.generation) ? &j : nullptr;
    }
    const Joint* get(JointId id) const { return const_cast<Impl*>(this)->get(id); }

    Body& end(uint32_t i) { return i == kWorldBody ? worldBody : bodies[i]; }
    const Body& end(uint32_t i) const { return i == kWorldBody ? worldBody : bodies[i]; }

    void updateJoinedPairs() {
        joinedPairs.clear();
        for (const Joint& j : joints) {
            if (!j.alive || j.collideConnected || j.bodyA == kWorldBody) { continue; }
            joinedPairs.push_back(pairKey(std::min(j.bodyA, j.bodyB), std::max(j.bodyA, j.bodyB)));
        }
        std::sort(joinedPairs.begin(), joinedPairs.end());
    }

    bool joined(uint32_t a, uint32_t b) const {
        return !joinedPairs.empty() &&
               std::binary_search(joinedPairs.begin(), joinedPairs.end(), pairKey(std::min(a, b), std::max(a, b)));
    }

    // Wakes the dynamic bodies of a joint, without spreading through a static or kinematic end.
    void wakeJoint(const Joint& j) {
        std::vector<uint32_t> start;
        for (uint32_t i : { j.bodyA, j.bodyB }) {
            if (i != kWorldBody && bodies[i].type == MotionType::Dynamic) { start.push_back(i); }
        }
        wakeBodies(start);
    }

    void destroyJoint(uint32_t ji) {
        Joint& j = joints[ji];
        wakeJoint(j);
        j.alive = false;
        ++j.generation;
        freeJoints.push_back(ji);
        updateJoinedPairs();
    }

    BodyId bodyId(uint32_t i) const { return { i, bodies[i].generation }; }

    JointInfo jointInfo(uint32_t i) const {
        const Joint& j = joints[i];
        const Body& a = end(j.bodyA);
        const Body& b = end(j.bodyB);
        JointInfo info;
        info.joint = { i, j.generation };
        info.type = j.type;
        if (j.bodyA != kWorldBody) { info.bodyA = bodyId(j.bodyA); }
        info.bodyB = bodyId(j.bodyB);
        info.anchorA = a.xf.apply(j.localA);
        info.anchorB = b.xf.apply(j.localB);
        info.axis = a.xf.q.rotate(j.axisA);
        info.value = jointValue(j, a, b);
        info.force = j.force;
        info.torque = j.torque;
        info.broken = j.broken;
        return info;
    }
    ShapeId shapeId(uint32_t i) const { return { i, shapes[i].generation }; }

    ContactEvent contactEvent(const Contact& c, bool begin) const {
        return { bodyId(c.bodyA), bodyId(c.bodyB), shapeId(c.shapeA), shapeId(c.shapeB), begin };
    }

    SensorEvent sensorEvent(uint64_t key, bool begin) const {
        const uint32_t s = (uint32_t)(key >> 32), v = (uint32_t)key;
        return { bodyId(shapes[s].body), bodyId(shapes[v].body), shapeId(s), shapeId(v), begin };
    }

    void updateMass(Body& b) {
        MassProperties mp;
        if (b.massOverride) {
            mp = b.overrideMass;
        } else {
            Vec3 weighted;
            std::vector<MassProperties> parts;
            for (uint32_t si : b.shapes) {
                if (shapes[si].sensor) { continue; }
                const MassProperties sm = shapeMass(shapes[si]);
                if (sm.mass <= 0.0f) { continue; }
                mp.mass += sm.mass;
                weighted += sm.center * sm.mass;
                parts.push_back(sm);
            }
            if (mp.mass > 0.0f) {
                mp.center = weighted / mp.mass;
                for (const MassProperties& sm : parts) {
                    const Vec3 d = sm.center - mp.center;
                    const float dd = dot(d, d);
                    mp.inertia = mp.inertia + sm.inertia + (Mat3::diagonal({ dd, dd, dd }) - Mat3::outer(d, d)) * sm.mass;
                }
            }
        }
        b.localCenter = mp.center;
        b.mass = mp.mass;
        b.inertia = mp.inertia;
        if (b.type == MotionType::Dynamic && b.mass <= 0.0f) {
            b.mass = 1.0f;
            b.inertia = Mat3{};
        }
        b.invMass = (b.type == MotionType::Dynamic && b.mass > 0.0f) ? 1.0f / b.mass : 0.0f;
        b.com = b.xf.apply(b.localCenter);
        b.extent = 0.0f;
        for (uint32_t si : b.shapes) {
            const Shape& s = shapes[si];
            if (s.kind == ShapeKind::HalfSpace || s.sensor) { continue; }
            b.extent = maxf(b.extent, s.extent + length(s.local.p - b.localCenter));
        }
        updateInverseMass(b);
    }

    float bodyExtent(const Body& b) const { return b.extent; }

    // How far any point of the body has moved since the step started.
    float movedThisStep(const Body& b) const {
        const Quat turn = b.q0.conjugate() * b.xf.q;
        const float s = minf(std::sqrt(turn.x * turn.x + turn.y * turn.y + turn.z * turn.z), 1.0f);
        return length(b.com - b.com0) + 2.0f * std::asin(s) * bodyExtent(b);
    }

    float relativeTravel(const Body& a, const Body& b, float dt) const {
        return (length(b.v - a.v) + length(a.w) * bodyExtent(a) + length(b.w) * bodyExtent(b)) * dt;
    }

    // A body's pose in the plane, placed once per round of collisions.
    const PlanarPose& poseOf(const Body& b, const Vec3& n) const {
        PlanarPose& p = b.planePose;
        if (p.stamp != poseStamp || !(p.normal.x == n.x && p.normal.y == n.y && p.normal.z == n.z)) {
            placeInPlane(p, b.xf, n);
            p.stamp = poseStamp;
        }
        return p;
    }

    // One manifold, or one for each way the triangles of a mesh face.
    void collidePair(const Shape& sa, const Shape& sb, float dt, std::vector<Manifold>& out,
                     TriangleBounds* bounds = nullptr) const {
        const Body& ba = bodies[sa.body];
        const Body& bb = bodies[sb.body];
        const float travel = relativeTravel(ba, bb, dt);
        const float margin = tol.speculative + travel;
        out.clear();
        Manifold m;
        Vec3 n;
        if (sa.kind != ShapeKind::Mesh && sb.kind != ShapeKind::Mesh && usePlanar(ba, bb, n) &&
            collidePlanar(sa, poseOf(ba, n), sb, poseOf(bb, n), margin, tol, m)) {
            if (m.count > 0) { out.push_back(m); }
            return;
        }
        const Transform xa = ba.xf * sa.local;
        const Transform xb = bb.xf * sb.local;
        const Motion mo{ ba.v, ba.w, ba.com, bb.v, bb.w, bb.com, dt, travel > tol.speculative };
        if (sa.kind == ShapeKind::Mesh) {
            collideMesh(sa, xa, sb, xb, margin, tol, mo, out, bounds);
            return;
        }
        if (sb.kind == ShapeKind::Mesh) {
            collideMesh(sb, xb, sa, xa, margin, tol, mo.swapped(), out, bounds);
            for (Manifold& x : out) { x.normal = -x.normal; }
            return;
        }
        collide(sa, xa, sb, xb, margin, tol, mo, m);
        if (m.count > 0) { out.push_back(m); }
    }

    // Points further apart than the speculative distance anchor on each body's own surface, so friction on a
    // fast body acts where it lands. Planar pairs keep the midpoint.
    Vec3 lift(const Body& a, const Body& b, const Vec3& normal, float separation) const {
        if (a.planar() || b.planar()) { return Vec3{}; }
        return separation > tol.speculative ? normal * (0.5f * separation) : Vec3{};
    }

    // The old contact of a pair whose normal is closest to n.
    static const Contact* closestNormal(const Contact* first, const Contact* last, const Vec3& n) {
        const Contact* best = nullptr;
        for (const Contact* c = first; c != last; ++c) {
            if (c->count > 0 && (!best || dot(c->normal, n) > dot(best->normal, n))) { best = c; }
        }
        return best;
    }

    bool meshTouches(const Shape& s, const Transform& xf, const Vec3& lo, const Vec3& hi) const {
        const TriangleMesh::Data* m = meshOf(s);
        if (!m) { return false; }
        Vec3 meshLo, meshHi;
        boundsInMesh(xf, lo, hi, meshLo, meshHi);
        bool any = false;
        m->query(meshLo, meshHi, [&](int t) {
            Vec3 c[3];
            m->corners(t, c);
            any = any || boxesOverlap(minVec(minVec(c[0], c[1]), c[2]), maxVec(maxVec(c[0], c[1]), c[2]), meshLo, meshHi);
        });
        return any;
    }

    // Marks the body's contacts that pick accepts as removed, for the next step to drop, and returns the
    // dynamic bodies they touched. When their shapes are going away, touching pairs end now.
    template <class Pick>
    std::vector<uint32_t> dropContacts(uint32_t body, bool ending, Pick&& pick) {
        if (!linksValid) { buildLinks(); }
        std::vector<uint32_t> touched;
        if (body + 1 >= linkStart.size()) { return touched; }
        for (uint32_t k = linkStart[body]; k < linkStart[body + 1]; ++k) {
            const Link& l = links[k];
            if (l.joint != UINT32_MAX || !pick(l)) { continue; }
            bool hit = false;
            auto it = std::lower_bound(contacts.begin(), contacts.end(), l.key,
                                       [](const Contact& x, uint64_t key) { return x.key < key; });
            if (ending && it != contacts.end() && it->key == l.key && it->touched) {
                queuedContactEvents.push_back(contactEvent(*it, false));
            }
            for (; it != contacts.end() && it->key == l.key; ++it) {
                if (ending) { it->touched = false; }
                if (it->removed) { continue; }
                it->removed = true;
                it->reusable = false;
                it->count = 0;
                hit = true;
            }
            if (hit && bodies[l.other].type == MotionType::Dynamic) { touched.push_back(l.other); }
        }
        return touched;
    }

    // All of a body's contacts, or those of one of its shapes.
    void removeContacts(uint32_t body, int shape, bool ending) {
        wakeBodies(dropContacts(body, ending, [&](const Link& l) {
            return shape < 0 || (uint32_t)(l.key >> 32) == (uint32_t)shape || (uint32_t)l.key == (uint32_t)shape;
        }));
    }

    // Contacts of a body that stopped being dynamic, with bodies that aren't dynamic either.
    void removeAnchoredContacts(uint32_t body) {
        dropContacts(body, false, [&](const Link& l) { return bodies[l.other].type != MotionType::Dynamic; });
    }

    // Ends the sensor overlaps of shapes that pick accepts, which are going away.
    template <class Pick>
    void dropSensorPairs(Pick&& pick) {
        size_t kept = 0;
        for (uint64_t key : sensorPairs) {
            if (pick((uint32_t)(key >> 32)) || pick((uint32_t)key)) {
                queuedSensorEvents.push_back(sensorEvent(key, false));
                continue;
            }
            sensorPairs[kept++] = key;
        }
        sensorPairs.resize(kept);
    }

    // Sleeping bodies an anchor was moved into have no contact with it yet.
    void pushOverlapping(uint32_t anchor, std::vector<uint32_t>& stack) const {
        const Vec3 grow{ tol.speculative, tol.speculative, tol.speculative };
        std::vector<uint32_t> found;
        for (uint32_t ai : bodies[anchor].shapes) {
            const Shape& sa = shapes[ai];
            if (sa.sensor) { continue; }
            const Transform xa = bodies[anchor].xf * sa.local;
            const Vec3 n = xa.q.rotate({ 0, 0, 1 });
            found.clear();
            forEachNear(sa, sa.lo, sa.hi, [&](uint32_t si) {
                const Shape& sb = shapes[si];
                const Body& b = bodies[sb.body];
                if (sb.sensor || b.type != MotionType::Dynamic || b.awake) { return; }
                if (!(sa.category & sb.mask) || !(sb.category & sa.mask)) { return; }
                const Vec3 blo = sb.lo - grow, bhi = sb.hi + grow;
                if (!boxesOverlap(sa.lo, sa.hi, blo, bhi)) { return; }
                if (sa.kind == ShapeKind::HalfSpace) {
                    const Vec3 c = (blo + bhi) * 0.5f, e = (bhi - blo) * 0.5f;
                    const float reach = e.x * std::fabs(n.x) + e.y * std::fabs(n.y) + e.z * std::fabs(n.z);
                    if (dot(c - xa.p, n) > reach) { return; }
                }
                if (sa.kind == ShapeKind::Mesh && !meshTouches(sa, xa, blo, bhi)) { return; }
                found.push_back(sb.body);
            });
            std::sort(found.begin(), found.end());
            found.erase(std::unique(found.begin(), found.end()), found.end());
            stack.insert(stack.end(), found.begin(), found.end());
        }
    }

    void buildLinks() {
        linkStart.assign(bodies.size() + 1, 0);
        links.clear();
        // The contacts of a mesh pair share a key and need one link.
        auto distinct = [&](size_t ci) { return ci == 0 || contacts[ci - 1].key != contacts[ci].key; };
        for (size_t ci = 0; ci < contacts.size(); ++ci) {
            if (!distinct(ci)) { continue; }
            ++linkStart[contacts[ci].bodyA + 1];
            ++linkStart[contacts[ci].bodyB + 1];
        }
        for (const Joint& j : joints) {
            if (!j.alive || j.bodyA == kWorldBody) { continue; }
            ++linkStart[j.bodyA + 1];
            ++linkStart[j.bodyB + 1];
        }
        for (size_t i = 1; i < linkStart.size(); ++i) { linkStart[i] += linkStart[i - 1]; }
        links.resize(linkStart.back());
        std::vector<uint32_t> fill(linkStart.begin(), linkStart.end() - 1);
        for (size_t ci = 0; ci < contacts.size(); ++ci) {
            if (!distinct(ci)) { continue; }
            const Contact& c = contacts[ci];
            links[fill[c.bodyA]++] = { c.key, c.bodyB, UINT32_MAX };
            links[fill[c.bodyB]++] = { c.key, c.bodyA, UINT32_MAX };
        }
        for (uint32_t ji = 0; ji < joints.size(); ++ji) {
            const Joint& j = joints[ji];
            if (!j.alive || j.bodyA == kWorldBody) { continue; }
            links[fill[j.bodyA]++] = { 0, j.bodyB, ji };
            links[fill[j.bodyB]++] = { 0, j.bodyA, ji };
        }
        linksValid = true;
    }

    // A contact links its bodies while it has points.
    bool linked(const Link& l) const {
        if (l.joint != UINT32_MAX) { return joints[l.joint].alive && !joints[l.joint].broken; }
        auto it = std::lower_bound(contacts.begin(), contacts.end(), l.key,
                                   [](const Contact& x, uint64_t k) { return x.key < k; });
        for (; it != contacts.end() && it->key == l.key; ++it) {
            if (it->count > 0) { return true; }
        }
        return false;
    }

    void pushLinked(uint32_t i, std::vector<uint32_t>& stack) const {
        if (i + 1 >= linkStart.size()) { return; }
        for (uint32_t k = linkStart[i]; k < linkStart[i + 1]; ++k) {
            const Link& l = links[k];
            const Body& o = bodies[l.other];
            if (o.type == MotionType::Dynamic && !o.awake && linked(l)) { stack.push_back(l.other); }
        }
    }

    // A static or kinematic start wakes the dynamic bodies touching it, without spreading past them.
    void wakeBodies(const std::vector<uint32_t>& start) {
        std::vector<uint32_t> stack;
        bool anchors = false;
        for (uint32_t i : start) {
            const Body& b = bodies[i];
            if (!b.alive) { continue; }
            if (b.type != MotionType::Dynamic) {
                anchors = true;
            } else if (!b.awake) {
                stack.push_back(i);
            }
        }
        if (stack.empty() && !anchors) { return; }

        if (!linksValid) { buildLinks(); }
        for (uint32_t i : start) {
            if (!bodies[i].alive || bodies[i].type == MotionType::Dynamic) { continue; }
            pushLinked(i, stack);
            pushOverlapping(i, stack);
        }
        while (!stack.empty()) {
            const uint32_t i = stack.back();
            stack.pop_back();
            Body& b = bodies[i];
            if (b.awake) { continue; }
            b.awake = true;
            b.sleepTime = 0.0f;
            updateInverseMass(b);
            pushLinked(i, stack);
        }
    }

    bool usePlanar(const Body& a, const Body& b, Vec3& n) const {
        const Body* first = nullptr;
        for (const Body* x : { &a, &b }) {
            if (x->type != MotionType::Dynamic) { continue; }
            if (!x->planar()) { return false; }
            if (!first) {
                first = x;
                n = x->planeNormal();
            } else if (std::fabs(dot(n, x->planeNormal())) < kAligned) {
                return false;
            }
        }
        return first != nullptr;
    }

    void updateContacts(float dt) {
        const Clock::time_point start = Clock::now();
        ++poseStamp;
        // Each shape's bounds grow by the speculative distance and how far its body can travel this step.
        std::vector<float>& grow = growth;
        grow.assign(bodies.size(), tol.speculative);
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            const Body& b = bodies[bi];
            if (!b.alive || b.type == MotionType::Static) { continue; }
            if (lengthSquared(b.v) == 0.0f && lengthSquared(b.w) == 0.0f) { continue; }
            grow[bi] += (length(b.v) + length(b.w) * bodyExtent(b)) * dt;
            const Vec3 g{ grow[bi], grow[bi], grow[bi] };
            fit(bi, b.lo - g, b.hi + g);
        }

        std::sort(refitted.begin(), refitted.end());
        refitted.erase(std::unique(refitted.begin(), refitted.end()), refitted.end());
        std::vector<uint64_t>& fresh = freshPairs;
        fresh.clear();
        for (uint32_t bi : refitted) {
            const Body& b = bodies[bi];
            if (!b.alive || b.proxy < 0) { continue; }
            tree.query(b.proxy, [&](uint32_t oi) {
                if (pairable(bi, oi)) { fresh.push_back(pairKey(std::min(bi, oi), std::max(bi, oi))); }
            });
        }
        refitted.clear();
        if (!fresh.empty()) {
            std::sort(fresh.begin(), fresh.end());
            fresh.erase(std::unique(fresh.begin(), fresh.end()), fresh.end());
            std::vector<uint64_t> merged;
            merged.reserve(leafPairs.size() + fresh.size());
            std::set_union(leafPairs.begin(), leafPairs.end(), fresh.begin(), fresh.end(), std::back_inserter(merged));
            leafPairs.swap(merged);
        }

        auto test = [&](uint32_t ia, uint32_t ib, std::vector<uint64_t>& out) {
            const Shape& sa = shapes[ia];
            const Shape& sb = shapes[ib];
            if (sa.sensor || sb.sensor) { return; }
            const Vec3 ga{ grow[sa.body], grow[sa.body], grow[sa.body] };
            const Vec3 gb{ grow[sb.body], grow[sb.body], grow[sb.body] };
            if (!boxesOverlap(sa.lo - ga, sa.hi + ga, sb.lo - gb, sb.hi + gb)) { return; }
            if (!(sa.category & sb.mask) || !(sb.category & sa.mask)) { return; }
            if (joined(sa.body, sb.body)) { return; }
            out.push_back(pairKey(std::min(ia, ib), std::max(ia, ib)));
        };

        // The shapes of each body near the other's bounds are paired, by their grown bounds.
        struct Near {
            uint32_t shape;
            Vec3 lo, hi;
        };
        std::vector<Near> nearA, nearB;
        auto near = [&](uint32_t bi, float g, const Vec3& lo, const Vec3& hi, std::vector<Near>& out) {
            out.clear();
            const Vec3 grown{ g, g, g };
            forEachShape(bi, lo, hi, [&](uint32_t si) {
                const Shape& s = shapes[si];
                if (s.sensor || !boxesOverlap(s.lo, s.hi, lo, hi)) { return; }
                out.push_back({ si, s.lo - grown, s.hi + grown });
            });
        };

        std::vector<uint64_t>& candidates = candidatePairs;
        candidates.clear();
        size_t kept = 0;
        for (uint64_t key : leafPairs) {
            const uint32_t ia = (uint32_t)(key >> 32), ib = (uint32_t)(key & 0xffffffffu);
            const Body& ba = bodies[ia];
            const Body& bb = bodies[ib];
            if (!ba.alive || !bb.alive || ba.proxy < 0 || bb.proxy < 0 || !pairable(ia, ib)) { continue; }
            if (!tree.overlap(ba.proxy, bb.proxy)) { continue; }
            leafPairs[kept++] = key;
            if (!movingAwake(ba) && !movingAwake(bb)) { continue; }
            const float g = grow[ia] + grow[ib];
            const Vec3 reach{ g, g, g };
            if (!boxesOverlap(ba.lo - reach, ba.hi + reach, bb.lo, bb.hi) || joined(ia, ib)) { continue; }
            near(ia, grow[ia], bb.lo - reach, bb.hi + reach, nearA);
            near(ib, grow[ib], ba.lo - reach, ba.hi + reach, nearB);
            for (const Near& x : nearA) {
                const Shape& sa = shapes[x.shape];
                for (const Near& y : nearB) {
                    if (!boxesOverlap(x.lo, x.hi, y.lo, y.hi)) { continue; }
                    const Shape& sb = shapes[y.shape];
                    if (!(sa.category & sb.mask) || !(sb.category & sa.mask)) { continue; }
                    candidates.push_back(pairKey(std::min(x.shape, y.shape), std::max(x.shape, y.shape)));
                }
            }
        }
        leafPairs.resize(kept);

        for (uint32_t pi : planes) {
            const Shape& sp = shapes[pi];
            const bool planeMoving = movingAwake(bodies[sp.body]);
            const float gp = grow[sp.body];
            const Transform xp = bodies[sp.body].xf * sp.local;
            const Vec3 n = xp.q.rotate({ 0, 0, 1 });
            const float offset = dot(n, xp.p);
            // Whether bounds grown by g reach below the plane.
            auto reaches = [&](const Vec3& lo, const Vec3& hi, float g) {
                const Vec3 c = (lo + hi) * 0.5f, e = (hi - lo) * 0.5f + Vec3{ g, g, g };
                return dot(c, n) - offset <= e.x * std::fabs(n.x) + e.y * std::fabs(n.y) + e.z * std::fabs(n.z);
            };
            for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
                const Body& b = bodies[bi];
                if (!b.alive || bi == sp.body || (!planeMoving && !movingAwake(b))) { continue; }
                const float g = gp + grow[bi];
                if (b.proxy >= 0 && !reaches(b.lo, b.hi, g)) { continue; }
                for (uint32_t si : b.shapes) {
                    const Shape& s = shapes[si];
                    if (s.kind != ShapeKind::HalfSpace && !reaches(s.lo, s.hi, g)) { continue; }
                    if (pairable(sp, s)) { test(pi, si, candidates); }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
        profile.candidatePairs = (int)candidates.size();
        profile.broadphase = millisecondsSince(start);
        const Clock::time_point narrowStart = Clock::now();

        std::vector<Contact>& next = nextContacts;
        next.clear();
        triangleBounds.clear();
        boundRuns.clear();
        std::vector<uint32_t> toWake;
        std::vector<Manifold> found;
        const float matchSq = tol.speculative * tol.speculative;
        // Both lists are in key order.
        auto old = contacts.begin();
        for (uint64_t key : candidates) {
            const uint32_t ia = (uint32_t)(key >> 32), ib = (uint32_t)(key & 0xffffffffu);
            const Shape& sa = shapes[ia];
            const Shape& sb = shapes[ib];
            const Body& ba = bodies[sa.body];
            const Body& bb = bodies[sb.body];
            while (old != contacts.end() && old->key < key) { ++old; }
            auto oldEnd = old;
            while (oldEnd != contacts.end() && oldEnd->key == key) { ++oldEnd; }
            const bool hasOld = old != oldEnd;
            bool keep = hasOld;
            for (auto it = old; it != oldEnd && keep; ++it) { keep = stillFits(*it, ba, bb); }
            const bool touched = hasOld && old->touched;

            const size_t first = next.size();
            if (keep) {
                for (auto it = old; it != oldEnd; ++it) {
                    next.push_back(*it);
                    measureAgain(ba, bb, next.back());
                }
            } else {
                // A fast pair is collided again within the step.
                const bool record = (sa.kind == ShapeKind::Mesh || sb.kind == ShapeKind::Mesh) &&
                                    relativeTravel(ba, bb, dt) > tol.speculative;
                const uint32_t begin = (uint32_t)triangleBounds.size();
                TriangleBounds bounds;
                bounds.found = &triangleBounds;
                collidePair(sa, sb, dt, found, record ? &bounds : nullptr);
                if (triangleBounds.size() > begin) {
                    std::sort(triangleBounds.begin() + begin, triangleBounds.end(),
                              [](const TriangleBound& x, const TriangleBound& y) { return x.triangle < y.triangle; });
                    boundRuns.push_back({ key, begin, (uint32_t)triangleBounds.size() });
                }
                if (found.empty()) { continue; }
            }
            for (size_t mi = 0; !keep && mi < found.size(); ++mi) {
                const Manifold& m = found[mi];
                Contact c;
                c.key = key;
                c.part = (uint32_t)mi;
                c.shapeA = ia;
                c.shapeB = ib;
                c.bodyA = sa.body;
                c.bodyB = sb.body;
                c.touched = touched;
                c.normal = m.normal;
                c.tangent1 = perpendicular(m.normal);
                c.tangent2 = cross(m.normal, c.tangent1);
                c.count = m.count;
                for (int i = 0; i < m.count; ++i) {
                    ContactPoint& p = c.points[i];
                    const Vec3 half = m.normal * (0.5f * m.points[i].separation);
                    p.surfaceA = ba.xf.applyInverse(m.points[i].position - half);
                    p.surfaceB = bb.xf.applyInverse(m.points[i].position + half);
                    p.localA = ba.xf.applyInverse(m.points[i].position);
                    p.separation = m.points[i].separation;
                    const Vec3 up = lift(ba, bb, m.normal, p.separation);
                    p.anchorA = m.points[i].position - up - ba.com;
                    p.anchorB = m.points[i].position + up - bb.com;
                }
                remember(ba, bb, c);
                const Contact* olds = contacts.data();
                if (const Contact* from = closestNormal(olds + (old - contacts.begin()), olds + (oldEnd - contacts.begin()), c.normal)) {
                    carryImpulses(*from, c, matchSq);
                }
                next.push_back(c);
            }
            for (size_t i = first; i < next.size(); ++i) {
                next[i].friction = frictionMix({ sa.friction, ba.type }, { sb.friction, bb.type });
                next[i].restitution = restitutionMix({ sa.restitution, ba.type }, { sb.restitution, bb.type });
                next[i].rollingResistance = maxf(sa.rollingResistance, sb.rollingResistance);
            }

            if (ba.type == MotionType::Dynamic && !ba.awake && movingAwake(bb)) { toWake.push_back(sa.body); }
            if (bb.type == MotionType::Dynamic && !bb.awake && movingAwake(ba)) { toWake.push_back(sb.body); }
        }

        const size_t collided = next.size();
        // A touching pair that wasn't collided again or carried over has ended. Both lists are in key order.
        size_t probe = 0;
        uint64_t ended = 0;
        auto lose = [&](const Contact& c) {
            if (!c.touched || c.key == ended) { return; }
            while (probe < collided && next[probe].key < c.key) { ++probe; }
            if (probe < collided && next[probe].key == c.key) { return; }
            contactEvents.push_back(contactEvent(c, false));
            ended = c.key;
        };
        for (const Contact& c : contacts) {
            const Body& a = bodies[c.bodyA];
            const Body& b = bodies[c.bodyB];
            if (c.removed || movingAwake(a) || movingAwake(b)) {
                lose(c);
                continue;
            }
            if (a.type != MotionType::Dynamic && b.type != MotionType::Dynamic) {
                lose(c);
                continue;
            }
            next.push_back(c);
        }
        std::inplace_merge(next.begin(), next.begin() + (std::ptrdiff_t)collided, next.end(), [](const Contact& a, const Contact& b) {
            return a.key != b.key ? a.key < b.key : a.part < b.part;
        });
        contacts.swap(next);
        linksValid = false;
        wakeBodies(toWake);
        profile.manifolds = (int)contacts.size();
        profile.narrowphase = millisecondsSince(narrowStart);
    }

    void remember(const Body& a, const Body& b, Contact& c) const {
        c.localNormal = a.xf.q.unrotate(c.normal);
        c.baseOffset = a.xf.q.unrotate(b.com - a.com);
        c.baseTurn = a.xf.q.conjugate() * b.xf.q;
        c.reusable = true;
    }

    bool stillFits(const Contact& c, const Body& a, const Body& b) const {
        if (!c.reusable) { return false; }
        const float reach = settings.contactReuseDistance * settings.lengthUnit;
        if (lengthSquared(a.xf.q.unrotate(b.com - a.com) - c.baseOffset) > reach * reach) { return false; }
        const Quat turn = a.xf.q.conjugate() * b.xf.q;
        return std::fabs(dot(turn, c.baseTurn)) >= std::cos(0.5f * settings.contactReuseAngle);
    }

    // Kept points stay where they were on each body. The normal turns with body A.
    void measureAgain(const Body& a, const Body& b, Contact& c) const {
        c.normal = a.xf.q.rotate(c.localNormal);
        c.tangent1 = perpendicular(c.normal);
        c.tangent2 = cross(c.normal, c.tangent1);
        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            const Vec3 pa = a.xf.apply(p.surfaceA), pb = b.xf.apply(p.surfaceB);
            const Vec3 mid = (pa + pb) * 0.5f;
            p.separation = dot(pb - pa, c.normal);
            const Vec3 up = lift(a, b, c.normal, p.separation);
            p.anchorA = mid - up - a.com;
            p.anchorB = mid + up - b.com;
            p.tangentImpulse = p.tangentImpulse - c.normal * dot(c.normal, p.tangentImpulse);
        }
    }

    // Contacts of every body in a pair closing faster than limit.
    // They get a second pass, so an impact reaches past the first contact.
    std::vector<uint32_t> impactContacts(const std::vector<uint32_t>& active, float limit) const {
        std::vector<uint8_t> struck(bodies.size(), 0);
        bool any = false;
        for (uint32_t ci : active) {
            const Contact& c = contacts[ci];
            for (int i = 0; i < c.count; ++i) {
                if (c.points[i].relativeVelocity >= -limit) { continue; }
                struck[c.bodyA] = bodies[c.bodyA].type == MotionType::Dynamic;
                struck[c.bodyB] = bodies[c.bodyB].type == MotionType::Dynamic;
                any = true;
                break;
            }
        }
        std::vector<uint32_t> out;
        if (!any) { return out; }
        for (uint32_t ci : active) {
            const Contact& c = contacts[ci];
            if (struck[c.bodyA] || struck[c.bodyB]) { out.push_back(ci); }
        }
        return out;
    }

    // A spinning body carries its contact points away from the overlap, so fast pairs are collided again within the step.
    void refreshContacts(const std::vector<uint32_t>& active, const std::vector<uint32_t>& fast, float remaining) {
        // Both lists are in contact order. Fast contacts are collided again below.
        size_t next = 0;
        for (uint32_t ci : active) {
            if (next < fast.size() && fast[next] == ci) {
                ++next;
                continue;
            }
            Contact& c = contacts[ci];
            rebaseContact(bodies[c.bodyA], bodies[c.bodyB], c, &axes[4 * ci], &planarAxes[4 * ci]);
        }
        for (Body& b : bodies) {
            b.dp = Vec3{};
            b.dq = Quat{};
            b.qBase = b.xf.q;
        }

        const float matchSq = tol.speculative * tol.speculative;
        ++poseStamp;
        std::vector<Manifold> found;
        std::vector<uint8_t> taken;
        for (size_t fi = 0; fi < fast.size(); ++fi) {
            const uint32_t ci = fast[fi];
            Contact& c = contacts[ci];
            const Body& a = bodies[c.bodyA];
            const Body& b = bodies[c.bodyB];
            // The contacts of a pair sit next to each other. The pair is collided once, and each contact takes
            // the new manifold whose normal is closest to its own.
            if (fi == 0 || contacts[fast[fi - 1]].key != c.key) {
                auto run = std::lower_bound(boundRuns.begin(), boundRuns.end(), c.key,
                                            [](const BoundRun& r, uint64_t key) { return r.key < key; });
                TriangleBounds bounds;
                if (run != boundRuns.end() && run->key == c.key) {
                    bounds.known = triangleBounds.data() + run->begin;
                    bounds.knownCount = (int)(run->end - run->begin);
                    bounds.moved = movedThisStep(a) + movedThisStep(b);
                }
                collidePair(shapes[c.shapeA], shapes[c.shapeB], remaining, found, &bounds);
                taken.assign(found.size(), 0);
            }
            int pick = -1;
            for (size_t mi = 0; mi < found.size(); ++mi) {
                if (taken[mi]) { continue; }
                if (pick < 0 || dot(found[mi].normal, c.normal) > dot(found[(size_t)pick].normal, c.normal)) { pick = (int)mi; }
            }
            if (pick < 0) {
                c.count = 0;
                continue;
            }
            taken[(size_t)pick] = 1;
            const Manifold& m = found[(size_t)pick];
            const Contact old = c;
            c.reusable = false;
            c.normal = m.normal;
            c.tangent1 = perpendicular(m.normal);
            c.tangent2 = cross(m.normal, c.tangent1);
            c.count = m.count;
            for (int i = 0; i < m.count; ++i) {
                ContactPoint& p = c.points[i];
                p = ContactPoint{};
                p.localA = a.xf.applyInverse(m.points[i].position);
                p.separation = m.points[i].separation;
                const Vec3 up = lift(a, b, m.normal, p.separation);
                p.anchorA = m.points[i].position - up - a.com;
                p.anchorB = m.points[i].position + up - b.com;
            }
            prepareContact(a, b, c, &axes[4 * ci], &planarAxes[4 * ci]);
            // Restitution needs the speed the pair approached at, from the nearest point the step started with.
            for (int i = 0; i < c.count; ++i) {
                ContactPoint& p = c.points[i];
                const int k = nearestPoint(old, p.localA);
                if (k < 0) { continue; }
                const ContactPoint& o = old.points[k];
                p.relativeVelocity = o.relativeVelocity;
                p.maxNormalImpulse = o.maxNormalImpulse;
            }
            carryImpulses(old, c, matchSq);
        }
    }

    void integrateVelocities(float h) {
        for (Body& b : bodies) {
            if (!b.alive || b.type != MotionType::Dynamic || !b.awake) { continue; }
            b.v = b.v + (settings.gravity * b.gravityScale + b.force * b.invMass) * h;
            b.v = b.linearProj * (b.v * (1.0f / (1.0f + h * b.linearDamping)));
            b.w = b.w + b.invInertiaWorld * b.torque * h;
            b.w = b.angularProj * (b.w * (1.0f / (1.0f + h * b.angularDamping)));
        }
    }

    void integratePositions(float h) {
        for (Body& b : bodies) {
            if (!b.alive || b.type == MotionType::Static) { continue; }
            if (b.type == MotionType::Dynamic && !b.awake) { continue; }
            b.com += b.v * h;
            b.dp += b.v * h;
            b.xf.q = integrate(b.xf.q, b.w, h);
            b.dq = b.xf.q * b.qBase.conjugate();
            b.xf.p = b.com - b.xf.q.rotate(b.localCenter);
        }
    }

    void updateSleep(float dt) {
        const float threshold = settings.sleepThreshold * settings.lengthUnit;
        for (Body& b : bodies) {
            if (!b.alive || b.type != MotionType::Dynamic || !b.awake) { continue; }
            const float speed = maxf(length(b.v), length(b.w) * bodyExtent(b));
            if (!b.canSleep || speed > threshold) {
                b.sleepTime = 0.0f;
            } else {
                b.sleepTime += dt;
            }
        }

        std::vector<uint32_t> parent(bodies.size());
        for (uint32_t i = 0; i < parent.size(); ++i) { parent[i] = i; }
        auto find = [&](uint32_t x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        for (const Contact& c : contacts) {
            if (c.count == 0) { continue; }
            Body& a = bodies[c.bodyA];
            Body& b = bodies[c.bodyB];
            const bool da = a.type == MotionType::Dynamic && a.awake;
            const bool db = b.type == MotionType::Dynamic && b.awake;
            if (da && db) {
                const uint32_t ra = find(c.bodyA), rb = find(c.bodyB);
                if (ra != rb) { parent[std::max(ra, rb)] = std::min(ra, rb); }
            } else if (da && b.type == MotionType::Kinematic && movingAwake(b)) {
                a.sleepTime = 0.0f;
            } else if (db && a.type == MotionType::Kinematic && movingAwake(a)) {
                b.sleepTime = 0.0f;
            }
        }
        for (const Joint& j : joints) {
            if (!j.alive || j.broken) { continue; }
            Body& a = end(j.bodyA);
            Body& b = end(j.bodyB);
            const bool da = a.type == MotionType::Dynamic && a.awake;
            const bool db = b.type == MotionType::Dynamic && b.awake;
            if (da && db) {
                const uint32_t ra = find(j.bodyA), rb = find(j.bodyB);
                if (ra != rb) { parent[std::max(ra, rb)] = std::min(ra, rb); }
            } else if (da && b.type == MotionType::Kinematic && movingAwake(b)) {
                a.sleepTime = 0.0f;
            } else if (db && a.type == MotionType::Kinematic && movingAwake(a)) {
                b.sleepTime = 0.0f;
            }
        }

        std::vector<float> islandTime(bodies.size(), FLT_MAX);
        profile.awakeBodies = 0;
        profile.awakeIslands = 0;
        for (uint32_t i = 0; i < bodies.size(); ++i) {
            Body& b = bodies[i];
            if (!b.alive || b.type != MotionType::Dynamic || !b.awake) { continue; }
            const uint32_t r = find(i);
            b.island = r;
            ++profile.awakeBodies;
            if (r == i) { ++profile.awakeIslands; }
            islandTime[r] = minf(islandTime[r], b.sleepTime);
        }
        for (uint32_t i = 0; i < bodies.size(); ++i) {
            Body& b = bodies[i];
            if (!b.alive || b.type != MotionType::Dynamic || !b.awake) { continue; }
            if (islandTime[find(i)] < settings.timeToSleep) { continue; }
            b.awake = false;
            b.v = Vec3{};
            b.w = Vec3{};
            updateInverseMass(b);
        }
    }

    float reach(const Shape& s, const Body& b) const { return length(s.local.p - b.localCenter) + s.extent; }

    // Distance from the centroid to the nearest face the body can move against.
    float innerExtent(const Shape& s, const Body& b) const {
        if ((s.kind != ShapeKind::Convex && s.kind != ShapeKind::Cylinder) || !b.planar()) { return s.minExtent; }
        const Vec3 n = s.local.q.unrotate(b.xf.q.unrotate(b.planeNormal()));
        if (s.kind == ShapeKind::Cylinder) { return std::fabs(n.z) > 0.9f ? s.radius : s.minExtent; }
        // The normal only drifts in the shape's frame, so the faces picked stay the same.
        if (s.inner >= 0.0f && lengthSquared(n - s.innerKey) <= 1e-10f) { return s.inner; }
        float inner = FLT_MAX;
        for (size_t i = 0; i < s.poly.normals.size(); ++i) {
            if (std::fabs(dot(s.poly.normals[i], n)) > 0.9f) { continue; }
            inner = minf(inner, s.poly.offsets[i] - dot(s.poly.normals[i], s.centroid));
        }
        s.innerKey = n;
        s.inner = inner < FLT_MAX ? inner : s.minExtent;
        return s.inner;
    }

    // Speculative contacts cover motion at the velocities the step started with. A shape that moved well
    // past that sweeps against everything near its path. Only a sweep that brings the shape's core into
    // the other shape moves the body back, to where the whole shape first touched it.
    void solveContinuous() {
        profile.sweptPairs = 0;
        profile.impacts = 0;
        if (!settings.continuous) { return; }

        std::vector<Sweep> sweeps(bodies.size());
        std::vector<uint8_t> fast(shapes.size(), 0);
        bool any = false;
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            const Body& b = bodies[bi];
            if (!b.alive) { continue; }
            sweeps[bi] = makeSweep(b.localCenter, b.com0, b.q0, b.com, b.xf.q);
            if (b.type != MotionType::Dynamic || !b.awake) { continue; }
            const float moved = length(b.com - b.com0);
            for (uint32_t si : b.shapes) {
                const Shape& s = shapes[si];
                if (s.kind == ShapeKind::Mesh || s.sensor) { continue; }
                const float travel = moved + sweeps[bi].angle * reach(s, b);
                if (travel > b.expectedTravel && travel > b.expectedTravel + 0.5f * innerExtent(s, b)) {
                    fast[si] = 1;
                    any = true;
                }
            }
        }
        if (!any) { return; }

        // Bounds over the whole step, found as needed.
        std::vector<uint8_t> swept(shapes.size(), 0);
        std::vector<Vec3> sweptLo(shapes.size()), sweptHi(shapes.size());
        auto sweptBounds = [&](uint32_t si) {
            if (swept[si]) { return; }
            swept[si] = 1;
            const Shape& s = shapes[si];
            sweptLo[si] = s.lo;
            sweptHi[si] = s.hi;
            if (bodies[s.body].type != MotionType::Static) {
                Vec3 lo, hi;
                shapeBounds(s, sweeps[s.body].at(0.0f) * s.local, lo, hi);
                sweptLo[si] = minVec(sweptLo[si], lo);
                sweptHi[si] = maxVec(sweptHi[si], hi);
            }
        };

        std::vector<std::pair<uint32_t, uint32_t>> pairs;
        for (uint32_t fi = 0; fi < shapes.size(); ++fi) {
            if (!fast[fi]) { continue; }
            sweptBounds(fi);
            const Shape& sf = shapes[fi];
            forEachNear(sf, sweptLo[fi], sweptHi[fi], [&](uint32_t oi) {
                const Shape& so = shapes[oi];
                if (so.sensor || so.body == sf.body || joined(so.body, sf.body)) { return; }
                sweptBounds(oi);
                if (!boxesOverlap(sweptLo[fi], sweptHi[fi], sweptLo[oi], sweptHi[oi])) { return; }
                if (!(sf.category & so.mask) || !(so.category & sf.mask)) { return; }
                pairs.push_back({ fi, oi });
            });
        }
        profile.sweptPairs = (int)pairs.size();

        std::vector<float> fraction(bodies.size(), 1.0f);
        std::vector<uint32_t> carrier(bodies.size(), UINT32_MAX);
        std::vector<uint32_t> group(bodies.size());
        for (uint32_t i = 0; i < group.size(); ++i) { group[i] = i; }
        auto find = [&](uint32_t x) {
            while (group[x] != x) {
                group[x] = group[group[x]];
                x = group[x];
            }
            return x;
        };
        const float target = tol.linearSlop;
        for (const auto& [fi, oi] : pairs) {
            const Shape& f = shapes[fi];
            const Shape& o = shapes[oi];
            const Body& bf = bodies[f.body];
            const Body& bo = bodies[o.body];
            const Convex core{ &f.centroid, 1, 0.5f * innerExtent(f, bf), false };
            const float coreReach = length(f.local.apply(f.centroid) - bf.localCenter) + core.radius;
            auto sweep = [&](const Convex& other, float otherReach, float& t) {
                float coreT = 1.0f;
                if (timeOfImpact(core, f.local, sweeps[f.body], coreReach, other, o.local, sweeps[o.body], otherReach,
                                 target, coreT) != Impact::Hit) {
                    return false;
                }
                if (timeOfImpact(convexOf(f), f.local, sweeps[f.body], reach(f, bf), other, o.local, sweeps[o.body],
                                 otherReach, target, t) != Impact::Hit) {
                    t = coreT;
                }
                return true;
            };
            float t = 1.0f;
            if (o.kind == ShapeKind::Mesh) {
                const TriangleMesh::Data* m = meshOf(o);
                if (!m) { continue; }
                Vec3 lo0, hi0, lo1, hi1;
                boundsInMesh(sweeps[o.body].at(0.0f) * o.local, sweptLo[fi], sweptHi[fi], lo0, hi0);
                boundsInMesh(bo.xf * o.local, sweptLo[fi], sweptHi[fi], lo1, hi1);
                bool hit = false;
                m->query(minVec(lo0, lo1), maxVec(hi0, hi1), [&](int tri) {
                    Vec3 c[3];
                    m->corners(tri, c);
                    float triReach = 0.0f;
                    for (const Vec3& v : c) { triReach = maxf(triReach, length(o.local.apply(v) - bo.localCenter)); }
                    float tt = 1.0f;
                    if (sweep({ c, 3, 0.0f, false }, triReach, tt)) {
                        t = minf(t, tt);
                        hit = true;
                    }
                });
                if (!hit) { continue; }
            } else if (!sweep(convexOf(o), o.kind == ShapeKind::HalfSpace ? 0.0f : reach(o, bo), t)) {
                continue;
            }
            if (t < fraction[f.body]) {
                fraction[f.body] = t;
                carrier[f.body] = bo.type == MotionType::Kinematic ? o.body : UINT32_MAX;
            }
            if (bodies[o.body].type == MotionType::Dynamic) {
                const uint32_t ra = find(f.body), rb = find(o.body);
                if (ra != rb) { group[std::max(ra, rb)] = std::min(ra, rb); }
            }
        }

        // Dynamic bodies that hit each other move back together, so none is left inside another.
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            const uint32_t r = find(bi);
            if (fraction[bi] < fraction[r]) {
                fraction[r] = fraction[bi];
                carrier[r] = carrier[bi];
            }
        }
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            const uint32_t r = find(bi);
            const float t = fraction[r];
            Body& b = bodies[bi];
            if (t >= 1.0f || b.type != MotionType::Dynamic || !b.awake) { continue; }
            b.xf = sweeps[bi].at(t);
            // A kinematic body keeps moving after the hit, so what hit it rides along to the end of the step.
            if (carrier[r] != UINT32_MAX) {
                const Transform k = sweeps[carrier[r]].at(t);
                const Transform ride{ k.applyInverse(b.xf.p), normalize(k.q.conjugate() * b.xf.q) };
                b.xf = bodies[carrier[r]].xf * ride;
            }
            b.com = b.xf.apply(b.localCenter);
            updateInverseMass(b);
            place(bi);
            ++profile.impacts;
        }
    }

    // Counts touching contacts and reports the pairs that started or stopped touching, and the hits. A pair
    // starts when the solver pushed on it or it overlaps, and stops once neither holds and every point is
    // further apart than the linear slop. Contacts of sleeping pairs keep their last solve, so they stay put.
    void reportTouches() {
        // Only contacts solved this step have fresh approach speeds.
        auto solved = [&](const Contact& c) {
            const Body& a = bodies[c.bodyA];
            const Body& b = bodies[c.bodyB];
            return (a.type == MotionType::Dynamic && a.awake) || (b.type == MotionType::Dynamic && b.awake);
        };
        profile.touchingContacts = 0;
        size_t last = 0;
        for (size_t first = 0; first < contacts.size(); first = last) {
            const uint64_t key = contacts[first].key;
            const bool was = contacts[first].touched;
            const float apart = was ? tol.linearSlop : 0.0f;
            bool touching = false;
            float fastest = settings.hitThreshold * settings.lengthUnit;
            const Contact* struck = nullptr;
            int point = 0;
            for (last = first; last < contacts.size() && contacts[last].key == key; ++last) {
                const Contact& c = contacts[last];
                if (c.touching()) { ++profile.touchingContacts; }
                for (int i = 0; i < c.count; ++i) {
                    const ContactPoint& p = c.points[i];
                    if (p.maxNormalImpulse > 0.0f) {
                        touching = true;
                        if (-p.relativeVelocity > fastest && solved(c)) {
                            fastest = -p.relativeVelocity;
                            struck = &c;
                            point = i;
                        }
                    } else if (p.separation < apart) {
                        touching = true;
                    }
                }
            }
            if (touching != was) {
                contactEvents.push_back(contactEvent(contacts[first], touching));
                for (size_t k = first; k < last; ++k) { contacts[k].touched = touching; }
            }
            if (struck) {
                const Vec3 p = bodies[struck->bodyA].xf.apply(struck->points[point].localA);
                hitEvents.push_back({ bodyId(struck->bodyA), bodyId(struck->bodyB), shapeId(struck->shapeA),
                                      shapeId(struck->shapeB), p, struck->normal, fastest });
            }
        }
    }

    bool meshOverlaps(const Shape& mesh, const Shape& other) const {
        const TriangleMesh::Data* m = meshOf(mesh);
        if (!m || other.kind == ShapeKind::HalfSpace) { return false; }
        const Transform xm = bodies[mesh.body].xf * mesh.local;
        const Transform xo = bodies[other.body].xf * other.local;
        const Convex co = convexOf(other);
        Vec3 lo, hi;
        boundsInMesh(xm, other.lo, other.hi, lo, hi);
        bool any = false;
        m->query(lo, hi, [&](int t) {
            if (any) { return; }
            Vec3 c[3];
            m->corners(t, c);
            any = distance({ c, 3, 0.0f, false }, xm, co, xo).distance <= 0.0f;
        });
        return any;
    }

    // Two meshes, a mesh and a half-space, or two half-spaces never overlap, as they never collide.
    bool overlapping(const Shape& a, const Shape& b) const {
        if (a.kind == ShapeKind::Mesh) { return b.kind != ShapeKind::Mesh && meshOverlaps(a, b); }
        if (b.kind == ShapeKind::Mesh) { return meshOverlaps(b, a); }
        return distance(convexOf(a), bodies[a.body].xf * a.local, convexOf(b), bodies[b.body].xf * b.local).distance <=
               0.0f;
    }

    // Finds what overlaps each sensor at the end of the step, and reports what changed. A pair of shapes
    // that both stayed still keeps its last result.
    void updateSensors() {
        if (sensors.empty()) { return; }
        std::vector<uint64_t> found;
        std::vector<uint32_t> near;
        size_t old = 0;
        for (uint32_t si : sensors) {
            const Shape& s = shapes[si];
            const bool moving = sensorsStale || movingAwake(bodies[s.body]);
            while (old < sensorPairs.size() && (uint32_t)(sensorPairs[old] >> 32) < si) { ++old; }
            size_t oldEnd = old;
            while (oldEnd < sensorPairs.size() && (uint32_t)(sensorPairs[oldEnd] >> 32) == si) { ++oldEnd; }
            near.clear();
            forEachNear(s, s.lo, s.hi, [&](uint32_t vi) { near.push_back(vi); });
            std::sort(near.begin(), near.end());
            for (uint32_t vi : near) {
                const Shape& v = shapes[vi];
                if (v.sensor || v.body == s.body || !(s.category & v.mask) || !(v.category & s.mask)) { continue; }
                if (!boxesOverlap(s.lo, s.hi, v.lo, v.hi)) { continue; }
                const uint64_t key = pairKey(si, vi);
                bool inside;
                if (moving || movingAwake(bodies[v.body])) {
                    inside = overlapping(s, v);
                } else {
                    inside = std::binary_search(sensorPairs.begin() + (std::ptrdiff_t)old,
                                                sensorPairs.begin() + (std::ptrdiff_t)oldEnd, key);
                }
                if (inside) { found.push_back(key); }
            }
            old = oldEnd;
        }
        sensorsStale = false;

        size_t i = 0, j = 0;
        while (i < sensorPairs.size() || j < found.size()) {
            if (j == found.size() || (i < sensorPairs.size() && sensorPairs[i] < found[j])) {
                sensorEvents.push_back(sensorEvent(sensorPairs[i++], false));
            } else if (i == sensorPairs.size() || found[j] < sensorPairs[i]) {
                sensorEvents.push_back(sensorEvent(found[j++], true));
            } else {
                ++i;
                ++j;
            }
        }
        sensorPairs.swap(found);
    }

    // Turns the impulses each active joint applied over the step into a force and torque, and breaks the joints
    // past their limits, which then leave active.
    void breakJoints(std::vector<uint32_t>& active, float dt) {
        size_t kept = 0;
        for (uint32_t ji : active) {
            Joint& j = joints[ji];
            j.force = j.force / dt;
            j.torque = j.torque / dt;
            if (length(j.force) > j.breakForce || length(j.torque) > j.breakTorque) {
                j.broken = true;
                j.linearImpulse = j.angularImpulse = j.frictionTurn = Vec3{};
                j.axialImpulse = j.motorImpulse = j.lowerImpulse = j.upperImpulse = 0.0f;
                j.frictionImpulse = j.swingImpulse = 0.0f;
                brokenJoints.push_back(ji);
                continue;
            }
            active[kept++] = ji;
        }
        active.resize(kept);
    }

    void step(float dt) {
        if (dt <= 0.0f) { return; }
        const Clock::time_point start = Clock::now();
        const int n = std::max(1, settings.substeps);
        const float h = dt / (float)n;

        contactEvents.swap(queuedContactEvents);
        queuedContactEvents.clear();
        sensorEvents.swap(queuedSensorEvents);
        queuedSensorEvents.clear();
        hitEvents.clear();
        brokenJoints.clear();
        const size_t queued = contactEvents.size();

        std::vector<uint32_t> pulled;
        for (const Joint& j : joints) {
            if (!j.alive || j.broken) { continue; }
            const Body& a = end(j.bodyA);
            const Body& b = end(j.bodyB);
            if (a.type == MotionType::Dynamic && !a.awake && movingAwake(b)) { pulled.push_back(j.bodyA); }
            if (b.type == MotionType::Dynamic && !b.awake && movingAwake(a)) { pulled.push_back(j.bodyB); }
        }
        wakeBodies(pulled);

        updateContacts(dt);
        const size_t lost = contactEvents.size();

        for (Body& b : bodies) {
            if (!b.alive) { continue; }
            b.q0 = b.xf.q;
            b.qBase = b.xf.q;
            b.com0 = b.com;
            b.expectedTravel = (length(b.v) + length(b.w) * bodyExtent(b)) * dt;
            b.dp = Vec3{};
            b.dq = Quat{};
            // The world inertia stays fixed through the substeps to match the effective masses built from it.
            updateInverseMass(b);
        }

        // Fast pairs are collided again mid-step, and those without a planar body before every substep.
        std::vector<uint32_t>& active = activeContacts;
        std::vector<uint32_t>& fast = fastContacts;
        std::vector<uint32_t>& fastSolid = fastSolidContacts;
        active.clear();
        fast.clear();
        fastSolid.clear();
        for (uint32_t i = 0; i < contacts.size(); ++i) {
            const Contact& c = contacts[i];
            const Body& a = bodies[c.bodyA];
            const Body& b = bodies[c.bodyB];
            if ((a.type == MotionType::Dynamic && a.awake) || (b.type == MotionType::Dynamic && b.awake)) {
                active.push_back(i);
                if (relativeTravel(a, b, dt) > tol.speculative) {
                    fast.push_back(i);
                    if (!a.planar() && !b.planar()) { fastSolid.push_back(i); }
                }
            }
        }

        std::vector<uint32_t> activeJoints;
        for (uint32_t i = 0; i < joints.size(); ++i) {
            const Joint& j = joints[i];
            if (!j.alive || j.broken) { continue; }
            const Body& a = end(j.bodyA);
            const Body& b = end(j.bodyB);
            if ((a.type == MotionType::Dynamic && a.awake) || (b.type == MotionType::Dynamic && b.awake)) {
                activeJoints.push_back(i);
            }
        }

        SolverContext ctx;
        ctx.h = h;
        ctx.invH = 1.0f / h;
        const float hertz = minf(settings.contactHertz, 0.25f * ctx.invH);
        ctx.contactSoft = makeSoft(hertz, settings.contactDampingRatio, h);
        ctx.staticSoft = makeSoft(2.0f * hertz, settings.contactDampingRatio, h);
        ctx.jointSoft = makeSoft(minf(settings.jointHertz, 0.5f * ctx.invH), settings.jointDampingRatio, h);
        ctx.maxPushout = settings.maxPushoutSpeed * settings.lengthUnit;
        ctx.restitutionThreshold = settings.restitutionThreshold * settings.lengthUnit;

        const Clock::time_point solveStart = Clock::now();
        axes.resize(4 * contacts.size());
        planarAxes.resize(4 * contacts.size());
        for (uint32_t ji : activeJoints) { joints[ji].force = joints[ji].torque = Vec3{}; }
        SolverData data{ bodies, contacts, axes, planarAxes, joints, worldBody };
        prepareContacts(data, active);
        const std::vector<uint32_t> impact = impactContacts(active, ctx.restitutionThreshold);
        for (int i = 0; i < n; ++i) {
            const std::vector<uint32_t>& again = i == n / 2 ? fast : fastSolid;
            if (i > 0 && !again.empty()) { refreshContacts(active, again, (float)(n - i) * h); }
            integrateVelocities(h);
            warmStartJoints(data, activeJoints);
            warmStart(data, active);
            solveJoints(data, activeJoints, ctx, true);
            solveContacts(data, active, ctx, true);
            if (!impact.empty()) { solveContacts(data, impact, ctx, true); }
            integratePositions(h);
            solveJoints(data, activeJoints, ctx, false);
            solveContacts(data, active, ctx, false);
            if (!impact.empty()) { solveContacts(data, impact, ctx, false); }
        }
        applyRestitution(data, active, ctx);
        if (!activeJoints.empty()) {
            breakJoints(activeJoints, dt);
            settleJoints(data, activeJoints, kSettleIterations);
            for (Body& b : bodies) {
                if (b.alive && b.type == MotionType::Dynamic && b.awake) { updateInverseMass(b); }
            }
        }
        profile.solver = millisecondsSince(solveStart);
        reportTouches();
        if (lost > queued && contactEvents.size() > lost) {
            auto byPair = [](const ContactEvent& a, const ContactEvent& b) {
                return a.shapeA.index != b.shapeA.index ? a.shapeA.index < b.shapeA.index : a.shapeB.index < b.shapeB.index;
            };
            auto first = contactEvents.begin() + (std::ptrdiff_t)queued;
            std::inplace_merge(first, contactEvents.begin() + (std::ptrdiff_t)lost, contactEvents.end(), byPair);
        }

        for (Body& b : bodies) {
            if (!b.alive) { continue; }
            b.force = Vec3{};
            b.torque = Vec3{};
            if (b.hasTarget) {
                b.xf = b.target;
                b.com = b.xf.apply(b.localCenter);
                b.hasTarget = false;
            }
        }
        // Leaves of moved bodies cover the step's start and end, which the sweeps need.
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            Body& b = bodies[bi];
            if (!b.alive || b.type == MotionType::Static || (b.type == MotionType::Dynamic && !b.awake)) { continue; }
            const Vec3 lo = b.lo, hi = b.hi;
            if (bound(b)) { fit(bi, minVec(lo, b.lo), maxVec(hi, b.hi)); }
        }
        const Clock::time_point continuousStart = Clock::now();
        solveContinuous();
        profile.continuous = millisecondsSince(continuousStart);
        updateSensors();
        const Clock::time_point sleepStart = Clock::now();
        updateSleep(dt);
        profile.sleep = millisecondsSince(sleepStart);
        profile.total = millisecondsSince(start);
    }
};

World::World(const WorldSettings& settings) : impl_(std::make_unique<Impl>(settings)) {}
World::~World() = default;
World::World(const World& other) : impl_(std::make_unique<Impl>(*other.impl_)) {}

World& World::operator=(const World& other) {
    if (this != &other) { impl_ = std::make_unique<Impl>(*other.impl_); }
    return *this;
}

World::World(World&&) noexcept = default;
World& World::operator=(World&&) noexcept = default;

const WorldSettings& World::settings() const { return impl_->settings; }
void World::setSettings(const WorldSettings& settings) { impl_->apply(settings); }
void World::setGravity(const Vec3& gravity) { impl_->settings.gravity = gravity; }
void World::setFrictionMix(MixFunction mix) { impl_->frictionMix = mix ? mix : mixFriction; }
void World::setRestitutionMix(MixFunction mix) { impl_->restitutionMix = mix ? mix : mixRestitution; }

void World::step(float dt) { impl_->step(dt); }
const StepProfile& World::profile() const { return impl_->profile; }

BodyId World::createBody(const BodyDesc& d) {
    Impl& w = *impl_;
    uint32_t i;
    if (!w.freeBodies.empty()) {
        i = w.freeBodies.back();
        w.freeBodies.pop_back();
    } else {
        i = (uint32_t)w.bodies.size();
        w.bodies.emplace_back();
    }
    Body& b = w.bodies[i];
    const uint32_t gen = b.generation;
    b = Body{};
    b.generation = gen;
    b.alive = true;
    b.type = d.type;
    b.xf = { d.position, normalize(d.rotation) };
    b.v = d.type == MotionType::Static ? Vec3{} : d.linearVelocity;
    b.w = d.type == MotionType::Static ? Vec3{} : d.angularVelocity;
    b.linearDamping = d.linearDamping;
    b.angularDamping = d.angularDamping;
    b.gravityScale = d.gravityScale;
    b.locks = d.locks;
    b.lockFrame = normalize(d.lockFrame);
    b.awake = d.awake;
    b.canSleep = d.canSleep;
    b.userData = d.userData;
    w.updateMass(b);
    b.v = b.linearProj * b.v;
    b.w = b.angularProj * b.w;
    return w.bodyId(i);
}

void World::destroyBody(BodyId id) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    w.removeContacts(id.index, -1, true);
    for (uint32_t ji = 0; ji < w.joints.size(); ++ji) {
        const Joint& j = w.joints[ji];
        if (j.alive && (j.bodyA == id.index || j.bodyB == id.index)) { w.destroyJoint(ji); }
    }
    if (!w.sensorPairs.empty()) {
        w.dropSensorPairs([&](uint32_t si) { return w.shapes[si].body == id.index; });
    }
    for (uint32_t si : b->shapes) {
        if (w.shapes[si].sensor) { w.sensors.erase(std::find(w.sensors.begin(), w.sensors.end(), si)); }
        w.untrack(si);
        w.shapes[si].alive = false;
        ++w.shapes[si].generation;
        w.freeShapes.push_back(si);
    }
    b->shapes.clear();
    b->alive = false;
    ++b->generation;
    w.retrack(id.index);
    w.freeBodies.push_back(id.index);
}

bool World::isValid(BodyId id) const { return impl_->get(id) != nullptr; }

std::vector<BodyId> World::bodies() const {
    std::vector<BodyId> out;
    for (uint32_t i = 0; i < impl_->bodies.size(); ++i) {
        if (impl_->bodies[i].alive) { out.push_back(impl_->bodyId(i)); }
    }
    return out;
}

ShapeId World::addShape(BodyId id, const ShapeDesc& d) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return {}; }
    uint32_t i;
    if (!w.freeShapes.empty()) {
        i = w.freeShapes.back();
        w.freeShapes.pop_back();
    } else {
        i = (uint32_t)w.shapes.size();
        w.shapes.emplace_back();
    }
    Shape& s = w.shapes[i];
    const uint32_t gen = s.generation;
    s = Shape{};
    s.generation = gen;
    s.alive = true;
    s.body = id.index;
    s.geometry = d.geometry;
    s.local = { d.position, normalize(d.rotation) };
    s.density = d.density;
    s.friction = d.friction;
    s.restitution = d.restitution;
    s.rollingResistance = d.rollingResistance;
    s.category = d.category;
    s.mask = d.mask;
    s.sensor = d.sensor;
    s.userData = d.userData;
    buildShape(s);

    Body& body = w.bodies[id.index];
    body.shapes.push_back(i);
    w.track(i);
    w.updateMass(body);
    if (s.sensor) { w.sensors.insert(std::lower_bound(w.sensors.begin(), w.sensors.end(), i), i); }
    w.sensorsStale = true;
    return w.shapeId(i);
}

void World::removeShape(ShapeId id) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    const uint32_t bi = s->body;
    w.removeContacts(bi, (int)id.index, true);
    if (!w.sensorPairs.empty()) { w.dropSensorPairs([&](uint32_t si) { return si == id.index; }); }
    if (s->sensor) { w.sensors.erase(std::find(w.sensors.begin(), w.sensors.end(), id.index)); }
    Body& b = w.bodies[bi];
    b.shapes.erase(std::remove(b.shapes.begin(), b.shapes.end(), id.index), b.shapes.end());
    w.untrack(id.index);
    w.retrack(bi);
    s->alive = false;
    ++s->generation;
    w.freeShapes.push_back(id.index);
    w.updateMass(b);
    if (!s->sensor) { w.wakeBodies({ bi }); }
}

bool World::isValid(ShapeId id) const { return impl_->get(id) != nullptr; }

BodyId World::shapeBody(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? impl_->bodyId(s->body) : BodyId{};
}

std::vector<ShapeId> World::shapes(BodyId id) const {
    std::vector<ShapeId> out;
    if (const Body* b = impl_->get(id)) {
        for (uint32_t si : b->shapes) { out.push_back(impl_->shapeId(si)); }
    }
    return out;
}

Geometry World::geometry(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->geometry : Geometry{};
}

Transform World::localTransform(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->local : Transform{};
}

float World::friction(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->friction : 0.0f;
}

float World::restitution(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->restitution : 0.0f;
}

float World::rollingResistance(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->rollingResistance : 0.0f;
}

Aabb World::bounds(ShapeId id) const {
    const Shape* s = impl_->get(id);
    if (!s) { return {}; }
    Aabb box;
    shapeBounds(*s, impl_->bodies[s->body].xf * s->local, box.lower, box.upper);
    return box;
}

void World::setGeometry(ShapeId id, const Geometry& geometry) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    w.untrack(id.index);
    s->geometry = geometry;
    buildShape(*s);
    w.track(id.index);
    for (Contact& c : w.contacts) {
        if (c.shapeA == id.index || c.shapeB == id.index) { c.reusable = false; }
    }
    w.updateMass(w.bodies[s->body]);
    w.wakeBodies({ s->body });
    w.sensorsStale = true;
}

void World::setLocalTransform(ShapeId id, const Vec3& position, const Quat& rotation) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    w.untrack(id.index);
    s->local = { position, normalize(rotation) };
    s->profile = PlanarProfile{};
    w.track(id.index);
    for (Contact& c : w.contacts) {
        if (c.shapeA == id.index || c.shapeB == id.index) { c.reusable = false; }
    }
    w.updateMass(w.bodies[s->body]);
    w.wakeBodies({ s->body });
    w.sensorsStale = true;
}

void World::setFriction(ShapeId id, float friction) {
    if (Shape* s = impl_->get(id)) { s->friction = friction; }
}

void World::setRestitution(ShapeId id, float restitution) {
    if (Shape* s = impl_->get(id)) { s->restitution = restitution; }
}

void World::setRollingResistance(ShapeId id, float rollingResistance) {
    if (Shape* s = impl_->get(id)) { s->rollingResistance = rollingResistance; }
}

void World::setFilter(ShapeId id, uint32_t category, uint32_t mask) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    s->category = category;
    s->mask = mask;
    w.removeContacts(s->body, (int)id.index, false);
    w.sensorsStale = true;
}

Transform World::transform(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->xf : Transform{};
}

Vec3 World::position(BodyId id) const { return transform(id).p; }
Quat World::rotation(BodyId id) const { return transform(id).q; }

Vec3 World::centerOfMass(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->com : Vec3{};
}

void World::setTransform(BodyId id, const Vec3& position, const Quat& rotation) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    b->xf = { position, normalize(rotation) };
    b->com = b->xf.apply(b->localCenter);
    b->hasTarget = false;
    w.place(id.index);
    w.wakeBodies({ id.index });
    w.sensorsStale = true;
}

Vec3 World::linearVelocity(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->v : Vec3{};
}

Vec3 World::angularVelocity(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->w : Vec3{};
}

void World::setLinearVelocity(BodyId id, const Vec3& v) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type == MotionType::Static) { return; }
    b->v = b->linearProj * v;
    w.wakeBodies({ id.index });
}

void World::setAngularVelocity(BodyId id, const Vec3& av) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type == MotionType::Static) { return; }
    b->w = b->angularProj * av;
    w.wakeBodies({ id.index });
}

void World::setKinematicTarget(BodyId id, const Vec3& position, const Quat& rotation, float dt) {
    Body* b = impl_->get(id);
    if (!b || b->type != MotionType::Kinematic || dt <= 0.0f) { return; }
    const Quat q = normalize(rotation);
    const Vec3 com = position + q.rotate(b->localCenter);
    b->v = (com - b->com) / dt;
    b->w = twofold::angularVelocity(b->xf.q, q, dt);
    b->target = { position, q };
    b->hasTarget = true;
}

MotionType World::motionType(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->type : MotionType::Static;
}

void World::setMotionType(BodyId id, MotionType type) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type == type) { return; }
    b->type = type;
    b->hasTarget = false;
    if (type == MotionType::Static) {
        b->v = Vec3{};
        b->w = Vec3{};
    }
    if (type == MotionType::Dynamic) { b->awake = true; }
    if (type != MotionType::Dynamic) { w.removeAnchoredContacts(id.index); }
    w.refitted.push_back(id.index);
    w.updateMass(*b);
    w.wakeBodies({ id.index });
}

void World::addForce(BodyId id, const Vec3& force, bool wake) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    if (wake) { w.wakeBodies({ id.index }); }
    if (b->awake) { b->force += force; }
}

void World::addForceAtPoint(BodyId id, const Vec3& force, const Vec3& point, bool wake) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    if (wake) { w.wakeBodies({ id.index }); }
    if (!b->awake) { return; }
    b->force += force;
    b->torque += cross(point - b->com, force);
}

void World::addTorque(BodyId id, const Vec3& torque, bool wake) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    if (wake) { w.wakeBodies({ id.index }); }
    if (b->awake) { b->torque += torque; }
}

void World::addImpulse(BodyId id, const Vec3& impulse, bool wake) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    if (wake) { w.wakeBodies({ id.index }); }
    if (b->awake) { b->v += b->invMassWorld * impulse; }
}

void World::addImpulseAtPoint(BodyId id, const Vec3& impulse, const Vec3& point, bool wake) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    if (wake) { w.wakeBodies({ id.index }); }
    if (!b->awake) { return; }
    b->v += b->invMassWorld * impulse;
    b->w += b->invInertiaWorld * cross(point - b->com, impulse);
}

MassProperties World::massProperties(BodyId id) const {
    const Body* b = impl_->get(id);
    if (!b) { return {}; }
    return { b->mass, b->localCenter, b->inertia };
}

void World::setMassProperties(BodyId id, const MassProperties& mass) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    b->massOverride = true;
    b->overrideMass = mass;
    w.updateMass(*b);
}

void World::resetMassProperties(BodyId id) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    b->massOverride = false;
    w.updateMass(*b);
}

uint8_t World::locks(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->locks : 0;
}

Quat World::lockFrame(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->lockFrame : Quat{};
}

void World::setLocks(BodyId id, uint8_t locks, const Quat& frame) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    b->locks = locks;
    b->lockFrame = normalize(frame);
    updateInverseMass(*b);
    b->v = b->linearProj * b->v;
    b->w = b->angularProj * b->w;
    w.wakeBodies({ id.index });
}

float World::linearDamping(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->linearDamping : 0.0f;
}

float World::angularDamping(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->angularDamping : 0.0f;
}

void World::setDamping(BodyId id, float linear, float angular) {
    if (Body* b = impl_->get(id)) {
        b->linearDamping = linear;
        b->angularDamping = angular;
    }
}

float World::gravityScale(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->gravityScale : 0.0f;
}

void World::setGravityScale(BodyId id, float scale) {
    if (Body* b = impl_->get(id)) { b->gravityScale = scale; }
}

bool World::isAwake(BodyId id) const {
    const Body* b = impl_->get(id);
    return b && b->awake;
}

void World::wake(BodyId id) {
    if (impl_->get(id)) { impl_->wakeBodies({ id.index }); }
}

void World::sleep(BodyId id) {
    Body* b = impl_->get(id);
    if (!b || b->type != MotionType::Dynamic) { return; }
    b->awake = false;
    b->v = Vec3{};
    b->w = Vec3{};
    b->sleepTime = 0.0f;
    updateInverseMass(*b);
}

bool World::canSleep(BodyId id) const {
    const Body* b = impl_->get(id);
    return b && b->canSleep;
}

void World::setCanSleep(BodyId id, bool canSleep) {
    Impl& w = *impl_;
    Body* b = w.get(id);
    if (!b) { return; }
    b->canSleep = canSleep;
    if (!canSleep) { w.wakeBodies({ id.index }); }
}

uint32_t World::island(BodyId id) const {
    const Body* b = impl_->get(id);
    return (b && b->type == MotionType::Dynamic) ? b->island : UINT32_MAX;
}

void* World::userData(BodyId id) const {
    const Body* b = impl_->get(id);
    return b ? b->userData : nullptr;
}

JointId World::createJoint(const JointDesc& d) {
    Impl& w = *impl_;
    if (!w.get(d.bodyB)) { return {}; }
    uint32_t ai = kWorldBody;
    if (d.bodyA.index != UINT32_MAX) {
        if (!w.get(d.bodyA) || d.bodyA.index == d.bodyB.index) { return {}; }
        ai = d.bodyA.index;
    }
    uint32_t i;
    if (!w.freeJoints.empty()) {
        i = w.freeJoints.back();
        w.freeJoints.pop_back();
    } else {
        i = (uint32_t)w.joints.size();
        w.joints.emplace_back();
    }
    Joint& j = w.joints[i];
    const uint32_t gen = j.generation;
    j = Joint{};
    j.generation = gen;
    j.alive = true;
    j.type = (JointType)d.kind.index();
    j.bodyA = ai;
    j.bodyB = d.bodyB.index;
    j.collideConnected = d.collideConnected;
    j.breakForce = d.breakForce;
    j.breakTorque = d.breakTorque;
    j.userData = d.userData;

    const Body& a = w.end(j.bodyA);
    const Body& b = w.end(j.bodyB);
    Vec3 anchorA, anchorB, axis{ 1, 0, 0 };
    Vec3 normal;
    if (const PointJoint* pj = std::get_if<PointJoint>(&d.kind)) {
        anchorA = anchorB = pj->anchor;
    } else if (const FixedJoint* fj = std::get_if<FixedJoint>(&d.kind)) {
        anchorA = anchorB = fj->anchor;
    } else if (const HingeJoint* hj = std::get_if<HingeJoint>(&d.kind)) {
        anchorA = anchorB = hj->anchor;
        axis = normalize(hj->axis);
        j.limited = hj->limited;
        j.lower = hj->lowerAngle;
        j.upper = hj->upperAngle;
        j.motorSpeed = hj->motorSpeed;
        j.maxMotor = hj->maxMotorTorque;
        j.maxFriction = hj->frictionTorque;
    } else if (const SliderJoint* sj = std::get_if<SliderJoint>(&d.kind)) {
        anchorA = anchorB = sj->anchor;
        axis = normalize(sj->axis);
        j.limited = sj->limited;
        j.lower = sj->lowerTranslation;
        j.upper = sj->upperTranslation;
        j.motorSpeed = sj->motorSpeed;
        j.maxMotor = sj->maxMotorForce;
        j.maxFriction = sj->frictionForce;
    } else if (const DistanceJoint* dj = std::get_if<DistanceJoint>(&d.kind)) {
        anchorA = dj->anchorA;
        anchorB = dj->anchorB;
        j.length = dj->length < 0.0f ? length(anchorB - anchorA) : dj->length;
        j.hertz = dj->hertz;
        j.dampingRatio = dj->dampingRatio;
        j.limited = dj->limited;
        j.lower = dj->minDistance;
        j.upper = dj->maxDistance;
    } else if (const SwingTwistJoint* tj = std::get_if<SwingTwistJoint>(&d.kind)) {
        anchorA = anchorB = tj->anchor;
        axis = normalize(tj->axis);
        normal = tj->normal;
        j.limited = tj->upperTwist - tj->lowerTwist < 2.0f * kPi;
        j.lower = tj->lowerTwist;
        j.upper = tj->upperTwist;
        j.swingY = coneSine(tj->swingY);
        j.swingZ = coneSine(tj->swingZ);
        j.maxFriction = tj->frictionTorque;
    }
    const Vec3 ref = perpendicular(axis);
    j.localA = a.xf.applyInverse(anchorA);
    j.localB = b.xf.applyInverse(anchorB);
    j.axisA = a.xf.q.unrotate(axis);
    j.axisB = b.xf.q.unrotate(axis);
    j.refA = a.xf.q.unrotate(ref);
    j.refB = b.xf.q.unrotate(ref);
    j.restTurn = a.xf.q.conjugate() * b.xf.q;
    if (j.type == JointType::SwingTwist) {
        Vec3 y = normal - axis * dot(normal, axis);
        y = lengthSquared(y) > 1e-12f ? normalize(y) : ref;
        const Vec3 z = cross(axis, y);
        Mat3 m;
        m.r[0] = { axis.x, y.x, z.x };
        m.r[1] = { axis.y, y.y, z.y };
        m.r[2] = { axis.z, y.z, z.z };
        const Quat frame = fromMat3(m);
        j.frameA = normalize(a.xf.q.conjugate() * frame);
        j.frameB = normalize(b.xf.q.conjugate() * frame);
    }

    if (!j.collideConnected && ai != kWorldBody) {
        auto between = [&](const Contact& c) {
            return (c.bodyA == j.bodyA && c.bodyB == j.bodyB) || (c.bodyA == j.bodyB && c.bodyB == j.bodyA);
        };
        for (size_t ci = 0; ci < w.contacts.size(); ++ci) {
            const Contact& c = w.contacts[ci];
            if (c.touched && between(c) && (ci == 0 || w.contacts[ci - 1].key != c.key)) {
                w.queuedContactEvents.push_back(w.contactEvent(c, false));
            }
        }
        w.contacts.erase(std::remove_if(w.contacts.begin(), w.contacts.end(), between), w.contacts.end());
    }
    w.updateJoinedPairs();
    w.linksValid = false;
    w.wakeJoint(j);
    return { i, j.generation };
}

void World::destroyJoint(JointId id) {
    if (impl_->get(id)) { impl_->destroyJoint(id.index); }
}

bool World::isValid(JointId id) const { return impl_->get(id) != nullptr; }

void World::setJointLimits(JointId id, bool limited, float lower, float upper) {
    Impl& w = *impl_;
    Joint* j = w.get(id);
    if (!j) { return; }
    j->limited = limited && (j->type != JointType::SwingTwist || upper - lower < 2.0f * kPi);
    j->lower = lower;
    j->upper = upper;
    j->lowerImpulse = 0.0f;
    j->upperImpulse = 0.0f;
    w.wakeJoint(*j);
}

void World::setJointMotor(JointId id, float speed, float maxForce) {
    Impl& w = *impl_;
    Joint* j = w.get(id);
    if (!j) { return; }
    j->motorSpeed = speed;
    j->maxMotor = maxForce;
    w.wakeJoint(*j);
}

void World::setJointSpring(JointId id, float hertz, float dampingRatio) {
    Impl& w = *impl_;
    Joint* j = w.get(id);
    if (!j) { return; }
    j->hertz = hertz;
    j->dampingRatio = dampingRatio;
    w.wakeJoint(*j);
}

void World::setJointFriction(JointId id, float maxFriction) {
    Impl& w = *impl_;
    Joint* j = w.get(id);
    if (!j) { return; }
    j->maxFriction = maxFriction;
    if (!(maxFriction > 0.0f)) {
        j->frictionImpulse = 0.0f;
        j->frictionTurn = Vec3{};
    }
    w.wakeJoint(*j);
}

void World::setJointSwing(JointId id, float swingY, float swingZ) {
    Impl& w = *impl_;
    Joint* j = w.get(id);
    if (!j) { return; }
    j->swingY = coneSine(swingY);
    j->swingZ = coneSine(swingZ);
    j->swingImpulse = 0.0f;
    w.wakeJoint(*j);
}

void World::setJointBreakLimits(JointId id, float force, float torque) {
    Joint* j = impl_->get(id);
    if (!j) { return; }
    j->breakForce = force;
    j->breakTorque = torque;
}

std::vector<JointInfo> World::joints() const {
    const Impl& w = *impl_;
    std::vector<JointInfo> out;
    for (uint32_t i = 0; i < w.joints.size(); ++i) {
        if (w.joints[i].alive) { out.push_back(w.jointInfo(i)); }
    }
    return out;
}

JointInfo World::jointInfo(JointId id) const {
    return impl_->get(id) ? impl_->jointInfo(id.index) : JointInfo{};
}

std::vector<JointId> World::brokenJoints() const {
    const Impl& w = *impl_;
    std::vector<JointId> out;
    for (uint32_t i : w.brokenJoints) {
        if (w.joints[i].alive) { out.push_back({ i, w.joints[i].generation }); }
    }
    return out;
}

void* World::userData(JointId id) const {
    const Joint* j = impl_->get(id);
    return j ? j->userData : nullptr;
}

void* World::userData(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->userData : nullptr;
}

std::vector<ContactInfo> World::contacts(bool includeSpeculative) const {
    const Impl& w = *impl_;
    std::vector<ContactInfo> out;
    for (const Contact& c : w.contacts) {
        if (c.removed || (!includeSpeculative && !c.touching())) { continue; }
        ContactInfo info;
        info.bodyA = w.bodyId(c.bodyA);
        info.bodyB = w.bodyId(c.bodyB);
        info.shapeA = w.shapeId(c.shapeA);
        info.shapeB = w.shapeId(c.shapeB);
        info.normal = c.normal;
        info.friction = c.friction;
        info.restitution = c.restitution;
        info.rollingResistance = c.rollingResistance;
        info.rollingImpulse = c.rollingImpulse;
        info.pointCount = c.count;
        for (int i = 0; i < c.count; ++i) {
            info.points[i] = w.bodies[c.bodyA].xf.apply(c.points[i].localA);
            info.separations[i] = c.points[i].separation;
            info.normalImpulses[i] = c.points[i].normalImpulse;
            info.tangentImpulses[i] = c.points[i].tangentImpulse;
        }
        out.push_back(info);
    }
    return out;
}

const std::vector<ContactEvent>& World::contactEvents() const { return impl_->contactEvents; }
const std::vector<HitEvent>& World::hitEvents() const { return impl_->hitEvents; }
const std::vector<SensorEvent>& World::sensorEvents() const { return impl_->sensorEvents; }

bool World::isSensor(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s && s->sensor;
}

bool World::castRay(const Vec3& origin, const Vec3& translation, RayHit& hit, uint32_t mask) const {
    const Impl& w = *impl_;
    const Vec3 end = origin + translation;
    const Vec3 lo = minVec(origin, end), hi = maxVec(origin, end);
    float best = 1.0f;
    bool found = false;
    // Shapes are tried in index order, so the lowest index wins a tie.
    std::vector<uint32_t> near(w.planes);
    auto onRay = [&](const BoundsTree& t, auto& visit) { t.raycast(origin, translation, visit); };
    w.tree.raycast(origin, translation, [&](uint32_t bi) {
        w.forEachShape(bi, onRay, [&](uint32_t si) { near.push_back(si); });
    });
    std::sort(near.begin(), near.end());
    for (uint32_t si : near) {
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask)) { continue; }
        if (!boxesOverlap(s.lo, s.hi, lo, hi)) { continue; }
        const Transform xf = w.bodies[s.body].xf * s.local;
        float fraction;
        Vec3 normal;
        if (!rayShape(s, xf, origin, translation, best, fraction, normal)) { continue; }
        if (found && fraction >= best) { continue; }
        best = fraction;
        found = true;
        hit.body = w.bodyId(s.body);
        hit.shape = w.shapeId(si);
        hit.point = origin + translation * fraction;
        hit.normal = normal;
        hit.fraction = fraction;
    }
    return found;
}

namespace {

// Query shapes this close to a shape, in meters, touch it.
constexpr float kTouching = 1e-5f;

bool convexShape(const Geometry& geometry, Shape& q) {
    if (std::holds_alternative<HalfSpace>(geometry) || std::holds_alternative<Mesh>(geometry)) { return false; }
    q.geometry = geometry;
    buildShape(q);
    return true;
}

// Sorts the candidates by shape index and keeps the ones keep accepts.
template <class Keep>
bool keepSorted(std::vector<ShapeId>& found, Keep&& keep) {
    std::sort(found.begin(), found.end(), [](const ShapeId& a, const ShapeId& b) { return a.index < b.index; });
    found.erase(std::remove_if(found.begin(), found.end(), [&](const ShapeId& id) { return !keep(id.index); }),
                found.end());
    return !found.empty();
}

}

bool World::castRayAll(const Vec3& origin, const Vec3& translation, std::vector<RayHit>& hits, uint32_t mask) const {
    const Impl& w = *impl_;
    const Vec3 end = origin + translation;
    const Vec3 lo = minVec(origin, end), hi = maxVec(origin, end);
    hits.clear();
    auto add = [&](uint32_t si) {
        RayHit h;
        h.shape.index = si;
        hits.push_back(h);
    };
    for (uint32_t si : w.planes) { add(si); }
    auto onRay = [&](const BoundsTree& t, auto& visit) { t.raycast(origin, translation, visit); };
    w.tree.raycast(origin, translation, [&](uint32_t bi) { w.forEachShape(bi, onRay, add); });
    std::sort(hits.begin(), hits.end(), [](const RayHit& a, const RayHit& b) { return a.shape.index < b.shape.index; });
    size_t kept = 0;
    for (size_t i = 0; i < hits.size(); ++i) {
        const uint32_t si = hits[i].shape.index;
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask) || !boxesOverlap(s.lo, s.hi, lo, hi)) { continue; }
        const Transform xf = w.bodies[s.body].xf * s.local;
        float fraction;
        Vec3 normal;
        if (!rayShape(s, xf, origin, translation, 1.0f, fraction, normal)) { continue; }
        RayHit& h = hits[kept++];
        h.body = w.bodyId(s.body);
        h.shape = w.shapeId(si);
        h.point = origin + translation * fraction;
        h.normal = normal;
        h.fraction = fraction;
    }
    hits.resize(kept);
    std::sort(hits.begin(), hits.end(), [](const RayHit& a, const RayHit& b) {
        return a.fraction != b.fraction ? a.fraction < b.fraction : a.shape.index < b.shape.index;
    });
    return !hits.empty();
}

bool World::castGeometry(const Geometry& geometry, const Transform& pose, const Vec3& translation, CastHit& hit,
                         uint32_t mask) const {
    const Impl& w = *impl_;
    Shape q;
    if (!convexShape(geometry, q)) { return false; }
    const float touch = kTouching * w.settings.lengthUnit;
    Vec3 lo, hi;
    shapeBounds(q, pose, lo, hi);
    const Vec3 center = (lo + hi) * 0.5f;
    const Vec3 reach = (hi - lo) * 0.5f + Vec3{ touch, touch, touch };
    // Shapes are tried in index order, so the lowest index wins a tie.
    std::vector<uint32_t> near(w.planes);
    auto onPath = [&](const BoundsTree& t, auto& visit) { t.boxcast(center, translation, reach, visit); };
    w.tree.boxcast(center, translation, reach, [&](uint32_t bi) {
        w.forEachShape(bi, onPath, [&](uint32_t si) { near.push_back(si); });
    });
    std::sort(near.begin(), near.end());
    float best = 1.0f;
    bool found = false;
    for (uint32_t si : near) {
        if (found && best == 0.0f) { break; }
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask)) { continue; }
        if (s.kind != ShapeKind::HalfSpace && !segmentHitsBox(s.lo - reach, s.hi + reach, center, translation, best)) {
            continue;
        }
        const Transform xf = w.bodies[s.body].xf * s.local;
        float fraction = best;
        Vec3 point, normal;
        const Impact r = castShape(s, xf, q, pose, translation, touch, w.tol, fraction, point, normal);
        if (r == Impact::Miss || (found && fraction >= best)) { continue; }
        best = fraction;
        found = true;
        hit.body = w.bodyId(s.body);
        hit.shape = w.shapeId(si);
        hit.point = point;
        hit.normal = normal;
        hit.fraction = fraction;
        hit.initialOverlap = r == Impact::Touching;
    }
    return found;
}

bool World::overlapGeometry(const Geometry& geometry, const Transform& pose, std::vector<ShapeId>& found,
                            uint32_t mask) const {
    const Impl& w = *impl_;
    found.clear();
    Shape q;
    if (!convexShape(geometry, q)) { return false; }
    const float touch = kTouching * w.settings.lengthUnit;
    Vec3 lo, hi;
    shapeBounds(q, pose, lo, hi);
    lo -= Vec3{ touch, touch, touch };
    hi += Vec3{ touch, touch, touch };
    for (uint32_t si : w.planes) { found.push_back(w.shapeId(si)); }
    w.tree.query(lo, hi, [&](uint32_t bi) {
        w.forEachShape(bi, lo, hi, [&](uint32_t si) { found.push_back(w.shapeId(si)); });
    });
    return keepSorted(found, [&](uint32_t si) {
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask)) { return false; }
        if (s.kind != ShapeKind::HalfSpace && !boxesOverlap(s.lo, s.hi, lo, hi)) { return false; }
        return shapesOverlap(s, w.bodies[s.body].xf * s.local, q, pose, touch);
    });
}

bool World::overlapAabb(const Aabb& box, std::vector<ShapeId>& found, uint32_t mask) const {
    const Impl& w = *impl_;
    found.clear();
    for (uint32_t si : w.planes) { found.push_back(w.shapeId(si)); }
    w.tree.query(box.lower, box.upper, [&](uint32_t bi) {
        w.forEachShape(bi, box.lower, box.upper, [&](uint32_t si) { found.push_back(w.shapeId(si)); });
    });
    const Vec3 c = (box.lower + box.upper) * 0.5f, e = (box.upper - box.lower) * 0.5f;
    return keepSorted(found, [&](uint32_t si) {
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask)) { return false; }
        if (s.kind != ShapeKind::HalfSpace) { return boxesOverlap(s.lo, s.hi, box.lower, box.upper); }
        const Transform xf = w.bodies[s.body].xf * s.local;
        const Vec3 n = xf.q.rotate({ 0, 0, 1 });
        return dot(c - xf.p, n) <= e.x * std::fabs(n.x) + e.y * std::fabs(n.y) + e.z * std::fabs(n.z);
    });
}

bool World::overlapPoint(const Vec3& point, std::vector<ShapeId>& found, uint32_t mask) const {
    const Impl& w = *impl_;
    found.clear();
    const float touch = kTouching * w.settings.lengthUnit;
    const Vec3 lo = point - Vec3{ touch, touch, touch }, hi = point + Vec3{ touch, touch, touch };
    for (uint32_t si : w.planes) { found.push_back(w.shapeId(si)); }
    w.tree.query(lo, hi, [&](uint32_t bi) {
        w.forEachShape(bi, lo, hi, [&](uint32_t si) { found.push_back(w.shapeId(si)); });
    });
    return keepSorted(found, [&](uint32_t si) {
        const Shape& s = w.shapes[si];
        if (s.sensor || !(s.category & mask)) { return false; }
        if (s.kind != ShapeKind::HalfSpace && !boxesOverlap(s.lo, s.hi, lo, hi)) { return false; }
        return shapeContains(s, w.bodies[s.body].xf * s.local, point, touch);
    });
}

}
