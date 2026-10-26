#include <algorithm>
#include <cfloat>
#include <chrono>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kAligned = 0.9999f;

using Clock = std::chrono::steady_clock;

float millisecondsSince(Clock::time_point start) {
    return std::chrono::duration<float, std::milli>(Clock::now() - start).count();
}

float mixFriction(const MaterialInput& a, const MaterialInput& b) { return std::sqrt(a.value * b.value); }
float mixRestitution(const MaterialInput& a, const MaterialInput& b) { return std::fmax(a.value, b.value); }

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
        return;
    }
    b.invMassWorld = b.linearProj * b.invMass;
    const Mat3 r = b.xf.q.toMat3();
    const Mat3 iw = r * b.inertia * transpose(r);
    b.invInertiaWorld = invertSubspace(iw, ang, na);
}

bool movingAwake(const Body& b) {
    if (b.type == MotionType::Dynamic) { return b.awake; }
    if (b.type == MotionType::Kinematic) { return lengthSquared(b.v) > 0.0f || lengthSquared(b.w) > 0.0f; }
    return false;
}

struct Proxy {
    Vec3 lo, hi;
    uint32_t shape = 0;
};

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
    StepProfile profile;

    explicit Impl(const WorldSettings& s) { apply(s); }

    void apply(const WorldSettings& s) {
        settings = s;
        tol.linearSlop = 0.005f * s.lengthUnit;
        tol.speculative = 4.0f * tol.linearSlop;
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

    BodyId bodyId(uint32_t i) const { return { i, bodies[i].generation }; }
    ShapeId shapeId(uint32_t i) const { return { i, shapes[i].generation }; }

    void updateMass(Body& b) {
        MassProperties mp;
        if (b.massOverride) {
            mp = b.overrideMass;
        } else {
            Vec3 weighted;
            std::vector<MassProperties> parts;
            for (uint32_t si : b.shapes) {
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
        updateInverseMass(b);
    }

    float bodyExtent(const Body& b) const {
        float e = 0.0f;
        for (uint32_t si : b.shapes) {
            const Shape& s = shapes[si];
            if (s.kind == ShapeKind::HalfSpace) { continue; }
            e = std::fmax(e, s.extent + length(s.local.p - b.localCenter));
        }
        return e;
    }

    void removeContacts(uint32_t body, int shape) {
        std::vector<uint32_t> touched;
        contacts.erase(std::remove_if(contacts.begin(), contacts.end(), [&](const Contact& c) {
            const bool hit = shape >= 0 ? (c.shapeA == (uint32_t)shape || c.shapeB == (uint32_t)shape)
                                        : (c.bodyA == body || c.bodyB == body);
            if (hit) { touched.push_back(c.bodyA == body ? c.bodyB : c.bodyA); }
            return hit;
        }), contacts.end());
        wakeBodies(touched);
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

        std::vector<std::vector<uint32_t>> adjacency(bodies.size());
        for (const Contact& c : contacts) {
            if (c.count == 0) { continue; }
            adjacency[c.bodyA].push_back(c.bodyB);
            adjacency[c.bodyB].push_back(c.bodyA);
        }
        for (uint32_t i : start) {
            if (!bodies[i].alive || bodies[i].type == MotionType::Dynamic) { continue; }
            for (uint32_t j : adjacency[i]) {
                const Body& o = bodies[j];
                if (o.type == MotionType::Dynamic && !o.awake) { stack.push_back(j); }
            }
        }
        while (!stack.empty()) {
            const uint32_t i = stack.back();
            stack.pop_back();
            Body& b = bodies[i];
            if (b.awake) { continue; }
            b.awake = true;
            b.sleepTime = 0.0f;
            updateInverseMass(b);
            for (uint32_t j : adjacency[i]) {
                const Body& o = bodies[j];
                if (o.type == MotionType::Dynamic && !o.awake) { stack.push_back(j); }
            }
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
        std::vector<Proxy> proxies;
        for (uint32_t bi = 0; bi < bodies.size(); ++bi) {
            const Body& b = bodies[bi];
            if (!b.alive) { continue; }
            float grow = tol.speculative;
            if (b.type != MotionType::Static) { grow += (length(b.v) + length(b.w) * bodyExtent(b)) * dt; }
            for (uint32_t si : b.shapes) {
                Proxy p;
                p.shape = si;
                shapeBounds(shapes[si], b.xf * shapes[si].local, p.lo, p.hi);
                const Vec3 g{ grow, grow, grow };
                p.lo -= g;
                p.hi += g;
                proxies.push_back(p);
            }
        }
        std::sort(proxies.begin(), proxies.end(), [](const Proxy& a, const Proxy& b) {
            return a.lo.x < b.lo.x || (a.lo.x == b.lo.x && a.shape < b.shape);
        });

        std::vector<uint64_t> candidates;
        for (size_t i = 0; i < proxies.size(); ++i) {
            const Proxy& pa = proxies[i];
            for (size_t j = i + 1; j < proxies.size() && proxies[j].lo.x <= pa.hi.x; ++j) {
                const Proxy& pb = proxies[j];
                if (pa.lo.y > pb.hi.y || pb.lo.y > pa.hi.y || pa.lo.z > pb.hi.z || pb.lo.z > pa.hi.z) { continue; }
                const Shape& sa = shapes[pa.shape];
                const Shape& sb = shapes[pb.shape];
                if (sa.body == sb.body) { continue; }
                if (!(sa.category & sb.mask) || !(sb.category & sa.mask)) { continue; }
                const Body& ba = bodies[sa.body];
                const Body& bb = bodies[sb.body];
                if (ba.type != MotionType::Dynamic && bb.type != MotionType::Dynamic) { continue; }
                if (!movingAwake(ba) && !movingAwake(bb)) { continue; }
                candidates.push_back(pairKey(std::min(pa.shape, pb.shape), std::max(pa.shape, pb.shape)));
            }
        }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
        profile.candidatePairs = (int)candidates.size();
        profile.broadphase = millisecondsSince(start);
        const Clock::time_point narrowStart = Clock::now();

        std::vector<Contact> next;
        std::vector<uint32_t> toWake;
        const float matchSq = tol.speculative * tol.speculative;
        for (uint64_t key : candidates) {
            const uint32_t ia = (uint32_t)(key >> 32), ib = (uint32_t)(key & 0xffffffffu);
            const Shape& sa = shapes[ia];
            const Shape& sb = shapes[ib];
            const Body& ba = bodies[sa.body];
            const Body& bb = bodies[sb.body];
            const Transform xa = ba.xf * sa.local;
            const Transform xb = bb.xf * sb.local;

            const float motion = (length(bb.v - ba.v) + length(ba.w) * bodyExtent(ba) + length(bb.w) * bodyExtent(bb)) * dt;
            const float margin = tol.speculative + motion;

            Manifold m;
            Vec3 n;
            if (!usePlanar(ba, bb, n) || !collidePlanar(sa, xa, sb, xb, n, margin, tol, m)) {
                const Motion mo{ ba.v, ba.w, ba.com, bb.v, bb.w, bb.com, dt };
                collide(sa, xa, sb, xb, margin, tol, mo, m);
            }
            if (m.count == 0) { continue; }

            Contact c;
            c.key = key;
            c.shapeA = ia;
            c.shapeB = ib;
            c.bodyA = sa.body;
            c.bodyB = sb.body;
            c.normal = m.normal;
            c.tangent1 = perpendicular(m.normal);
            c.tangent2 = cross(m.normal, c.tangent1);
            c.friction = frictionMix({ sa.friction, ba.type }, { sb.friction, bb.type });
            c.restitution = restitutionMix({ sa.restitution, ba.type }, { sb.restitution, bb.type });
            c.count = m.count;

            auto old = std::lower_bound(contacts.begin(), contacts.end(), key,
                                        [](const Contact& x, uint64_t k) { return x.key < k; });
            const bool hasOld = old != contacts.end() && old->key == key;
            for (int i = 0; i < m.count; ++i) {
                ContactPoint& p = c.points[i];
                p.localA = ba.xf.applyInverse(m.points[i].position);
                p.separation = m.points[i].separation;
                p.anchorA = m.points[i].position - ba.com;
                p.anchorB = m.points[i].position - bb.com;
                if (!hasOld) { continue; }
                for (int k = 0; k < old->count; ++k) {
                    const ContactPoint& o = old->points[k];
                    if (lengthSquared(o.localA - p.localA) > matchSq) { continue; }
                    p.normalImpulse = o.normalImpulse;
                    p.tangentImpulse = o.tangentImpulse - c.normal * dot(c.normal, o.tangentImpulse);
                    break;
                }
            }

            if (ba.type == MotionType::Dynamic && !ba.awake && movingAwake(bb)) { toWake.push_back(sa.body); }
            if (bb.type == MotionType::Dynamic && !bb.awake && movingAwake(ba)) { toWake.push_back(sb.body); }
            next.push_back(c);
        }

        for (const Contact& c : contacts) {
            if (movingAwake(bodies[c.bodyA]) || movingAwake(bodies[c.bodyB])) { continue; }
            next.push_back(c);
        }
        std::sort(next.begin(), next.end(), [](const Contact& a, const Contact& b) { return a.key < b.key; });
        contacts.swap(next);
        wakeBodies(toWake);
        profile.manifolds = (int)contacts.size();
        profile.narrowphase = millisecondsSince(narrowStart);
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
            b.dq = b.xf.q * b.q0.conjugate();
            b.xf.p = b.com - b.xf.q.rotate(b.localCenter);
            if (b.type == MotionType::Dynamic) { updateInverseMass(b); }
        }
    }

    void updateSleep(float dt) {
        const float threshold = settings.sleepThreshold * settings.lengthUnit;
        for (Body& b : bodies) {
            if (!b.alive || b.type != MotionType::Dynamic || !b.awake) { continue; }
            const float speed = std::fmax(length(b.v), length(b.w) * bodyExtent(b));
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
            islandTime[r] = std::fmin(islandTime[r], b.sleepTime);
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

    void step(float dt) {
        if (dt <= 0.0f) { return; }
        const Clock::time_point start = Clock::now();
        const int n = std::max(1, settings.substeps);
        const float h = dt / (float)n;

        updateContacts(dt);

        for (Body& b : bodies) {
            if (!b.alive) { continue; }
            b.q0 = b.xf.q;
            b.dp = Vec3{};
            b.dq = Quat{};
            updateInverseMass(b);
        }

        std::vector<uint32_t> active;
        for (uint32_t i = 0; i < contacts.size(); ++i) {
            const Contact& c = contacts[i];
            const Body& a = bodies[c.bodyA];
            const Body& b = bodies[c.bodyB];
            if ((a.type == MotionType::Dynamic && a.awake) || (b.type == MotionType::Dynamic && b.awake)) { active.push_back(i); }
        }

        SolverContext ctx;
        ctx.h = h;
        ctx.invH = 1.0f / h;
        const float hertz = std::fmin(settings.contactHertz, 0.25f * ctx.invH);
        ctx.contactSoft = makeSoft(hertz, settings.contactDampingRatio, h);
        ctx.staticSoft = makeSoft(2.0f * hertz, settings.contactDampingRatio, h);
        ctx.maxPushout = settings.maxPushoutSpeed * settings.lengthUnit;
        ctx.restitutionThreshold = settings.restitutionThreshold * settings.lengthUnit;

        const Clock::time_point solveStart = Clock::now();
        prepareContacts(bodies, contacts, active);
        for (int i = 0; i < n; ++i) {
            integrateVelocities(h);
            warmStart(bodies, contacts, active);
            solveContacts(bodies, contacts, active, ctx, true);
            integratePositions(h);
            solveContacts(bodies, contacts, active, ctx, false);
        }
        applyRestitution(bodies, contacts, active, ctx);
        profile.solver = millisecondsSince(solveStart);
        profile.touchingContacts = 0;
        for (const Contact& c : contacts) {
            if (c.touching()) { ++profile.touchingContacts; }
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
    w.removeContacts(id.index, -1);
    for (uint32_t si : b->shapes) {
        w.shapes[si].alive = false;
        ++w.shapes[si].generation;
        w.freeShapes.push_back(si);
    }
    b->shapes.clear();
    b->alive = false;
    ++b->generation;
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
    s.category = d.category;
    s.mask = d.mask;
    s.userData = d.userData;
    buildShape(s);

    Body& body = w.bodies[id.index];
    body.shapes.push_back(i);
    w.updateMass(body);
    return w.shapeId(i);
}

void World::removeShape(ShapeId id) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    const uint32_t bi = s->body;
    w.removeContacts(bi, (int)id.index);
    Body& b = w.bodies[bi];
    b.shapes.erase(std::remove(b.shapes.begin(), b.shapes.end(), id.index), b.shapes.end());
    s->alive = false;
    ++s->generation;
    w.freeShapes.push_back(id.index);
    w.updateMass(b);
    w.wakeBodies({ bi });
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
    s->geometry = geometry;
    buildShape(*s);
    w.updateMass(w.bodies[s->body]);
    w.wakeBodies({ s->body });
}

void World::setFriction(ShapeId id, float friction) {
    if (Shape* s = impl_->get(id)) { s->friction = friction; }
}

void World::setRestitution(ShapeId id, float restitution) {
    if (Shape* s = impl_->get(id)) { s->restitution = restitution; }
}

void World::setFilter(ShapeId id, uint32_t category, uint32_t mask) {
    Impl& w = *impl_;
    Shape* s = w.get(id);
    if (!s) { return; }
    s->category = category;
    s->mask = mask;
    w.removeContacts(s->body, (int)id.index);
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
    w.wakeBodies({ id.index });
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

void* World::userData(ShapeId id) const {
    const Shape* s = impl_->get(id);
    return s ? s->userData : nullptr;
}

std::vector<ContactInfo> World::contacts(bool includeSpeculative) const {
    const Impl& w = *impl_;
    std::vector<ContactInfo> out;
    for (const Contact& c : w.contacts) {
        if (!includeSpeculative && !c.touching()) { continue; }
        ContactInfo info;
        info.bodyA = w.bodyId(c.bodyA);
        info.bodyB = w.bodyId(c.bodyB);
        info.shapeA = w.shapeId(c.shapeA);
        info.shapeB = w.shapeId(c.shapeB);
        info.normal = c.normal;
        info.friction = c.friction;
        info.restitution = c.restitution;
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

bool World::castRay(const Vec3& origin, const Vec3& translation, RayHit& hit, uint32_t mask) const {
    const Impl& w = *impl_;
    const Vec3 end = origin + translation;
    const Vec3 lo = minVec(origin, end), hi = maxVec(origin, end);
    float best = 1.0f;
    bool found = false;
    for (uint32_t si = 0; si < w.shapes.size(); ++si) {
        const Shape& s = w.shapes[si];
        if (!s.alive || !(s.category & mask)) { continue; }
        const Transform xf = w.bodies[s.body].xf * s.local;
        Vec3 slo, shi;
        shapeBounds(s, xf, slo, shi);
        if (slo.x > hi.x || shi.x < lo.x || slo.y > hi.y || shi.y < lo.y || slo.z > hi.z || shi.z < lo.z) { continue; }
        float fraction;
        Vec3 normal;
        if (!rayShape(s, xf, origin, translation, best, fraction, normal)) { continue; }
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

}
