#include "internal.h"

namespace twofold {

Axis makeAxis(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& d) {
    Axis x;
    x.armA = cross(ra, d);
    x.armB = cross(rb, d);
    x.linearA = a.invMassWorld * d;
    x.angularA = a.invInertiaWorld * x.armA;
    x.linearB = b.invMassWorld * d;
    x.angularB = b.invInertiaWorld * x.armB;
    const float k = dot(d, x.linearA) + dot(d, x.linearB) + dot(x.armA, x.angularA) + dot(x.armB, x.angularB);
    // A direction the bodies are locked against reaches k only through rounding.
    x.mass = k > 1e-6f * (a.invMass + b.invMass) ? 1.0f / k : 0.0f;
    return x;
}

float speedAlong(const Body& a, const Body& b, const Vec3& d, const Axis& x) {
    return dot(d, b.v - a.v) + dot(x.armB, b.w) - dot(x.armA, a.w);
}

void push(Body& a, Body& b, const Axis& x, float impulse) {
    a.v -= x.linearA * impulse;
    a.w -= x.angularA * impulse;
    b.v += x.linearB * impulse;
    b.w += x.angularB * impulse;
}

namespace {

void makeAxes(const Body& a, const Body& b, const Contact& c, PointAxes* axes) {
    for (int i = 0; i < c.count; ++i) {
        const ContactPoint& p = c.points[i];
        axes[i].normal = makeAxis(a, b, p.anchorA, p.anchorB, c.normal);
        axes[i].tangent1 = makeAxis(a, b, p.anchorA, p.anchorB, c.tangent1);
        axes[i].tangent2 = makeAxis(a, b, p.anchorA, p.anchorB, c.tangent2);
    }
}

bool anchored(const Body& b) { return b.type != MotionType::Dynamic || !b.awake; }

// Where the body has turned r since the substep's base pose.
Vec3 turned(const Body& b, const Vec3& r) {
    return b.dq.x == 0.0f && b.dq.y == 0.0f && b.dq.z == 0.0f ? r : b.dq.rotate(r);
}

constexpr float kAligned = 0.9999f;
// How far a contact normal may lean out of the plane and still be solved in it.
constexpr float kInPlane = 1e-3f;

// Whether every body the contact moves is a dynamic body locked to one plane, n, with the normal in it.
bool planarContact(const Body& a, const Body& b, const Vec3& normal, Vec3& n) {
    const Body* first = nullptr;
    for (const Body* x : { &a, &b }) {
        if (x->type == MotionType::Static) { continue; }
        if (x->type == MotionType::Kinematic || !x->planar()) { return false; }
        if (!first) {
            first = x;
            n = x->planeNormal();
        } else if (std::fabs(dot(n, x->planeNormal())) < kAligned) {
            return false;
        }
    }
    return first && std::fabs(dot(normal, n)) < kInPlane;
}

float planarMass(const Body& a, const Body& b, float armA, float armB) {
    const float k = a.planeInvMass + b.planeInvMass + a.planeInvInertia * armA * armA + b.planeInvInertia * armB * armB;
    return k > 1e-6f * (a.invMass + b.invMass) ? 1.0f / k : 0.0f;
}

void makePlanarAxes(const Body& a, const Body& b, const Contact& c, PlanarAxes* axes) {
    const Vec3 n = c.plane;
    for (int i = 0; i < c.count; ++i) {
        const ContactPoint& p = c.points[i];
        PlanarAxes& x = axes[i];
        x.normalA = dot(n, cross(p.anchorA, c.normal));
        x.normalB = dot(n, cross(p.anchorB, c.normal));
        x.normalMass = planarMass(a, b, x.normalA, x.normalB);
        x.tangentA = dot(n, cross(p.anchorA, c.tangent1));
        x.tangentB = dot(n, cross(p.anchorB, c.tangent1));
        x.tangentMass = planarMass(a, b, x.tangentA, x.tangentB);
    }
}

// The velocities of a planar contact's bodies, with spin about the plane normal, held while its points are solved.
struct PlanarState {
    Vec3 va, vb;
    float wa, wb;
    float ma, mb, ia, ib;

    PlanarState(const Body& a, const Body& b, const Vec3& n)
        : va(a.v), vb(b.v), wa(dot(a.w, n)), wb(dot(b.w, n)), ma(a.planeInvMass), mb(b.planeInvMass),
          ia(a.planeInvInertia), ib(b.planeInvInertia) {}

    float speed(const Vec3& d, float armA, float armB) const { return dot(d, vb - va) + armB * wb - armA * wa; }

    void push(const Vec3& d, float armA, float armB, float impulse) {
        va -= d * (ma * impulse);
        wa -= ia * armA * impulse;
        vb += d * (mb * impulse);
        wb += ib * armB * impulse;
    }

    void store(Body& a, Body& b, const Vec3& n, bool moveA, bool moveB) const {
        if (moveA) {
            a.v = va;
            a.w = n * wa;
        }
        if (moveB) {
            b.v = vb;
            b.w = n * wb;
        }
    }
};

// The cosine and sine of the turn since the substep's base pose, about n.
Vec2 planarTurn(const Body& b, const Vec3& n) {
    const Quat& q = b.dq;
    return { q.w * q.w - (q.x * q.x + q.y * q.y + q.z * q.z), 2.0f * q.w * dot(n, Vec3{ q.x, q.y, q.z }) };
}

void solvePlanar(Body& a, Body& b, Contact& c, const PlanarAxes* axes, const SolverContext& ctx, bool useBias) {
    const bool moveA = !anchored(a), moveB = !anchored(b);
    const Softness soft = (!moveA || !moveB) ? ctx.staticSoft : ctx.contactSoft;
    const float shift = dot(b.dp - a.dp, c.normal);
    const Vec2 ta = planarTurn(a, c.plane), tb = planarTurn(b, c.plane);
    PlanarState st(a, b, c.plane);
    bool changed = false;

    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        const PlanarAxes& x = axes[i];
        // An arm's tangent axis is its reach along the normal, and its normal axis its reach along the
        // tangent turned a quarter, so the turned anchors reach along the normal by these.
        const float reachA = ta.x * x.tangentA + ta.y * x.normalA;
        const float reachB = tb.x * x.tangentB + tb.y * x.normalB;
        const float sep = shift + reachB - reachA + p.adjustedSeparation;

        float bias = 0.0f, massScale = 1.0f, impulseScale = 0.0f;
        if (sep > 0.0f) {
            bias = sep * ctx.invH;
        } else if (useBias) {
            bias = maxf(soft.biasRate * sep, -ctx.maxPushout);
            massScale = soft.massScale;
            impulseScale = soft.impulseScale;
        }

        const float vn = st.speed(c.normal, x.normalA, x.normalB);
        float impulse = -x.normalMass * massScale * (vn + bias) - impulseScale * p.normalImpulse;
        const float total = maxf(p.normalImpulse + impulse, 0.0f);
        impulse = total - p.normalImpulse;
        p.normalImpulse = total;
        p.maxNormalImpulse = maxf(p.maxNormalImpulse, impulse);
        if (impulse != 0.0f) {
            st.push(c.normal, x.normalA, x.normalB, impulse);
            changed = true;
        }
    }

    for (int j = 0; j < c.count; ++j) {
        const int i = useBias ? j : c.count - 1 - j;
        ContactPoint& p = c.points[i];
        const PlanarAxes& x = axes[i];
        const float old = dot(p.tangentImpulse, c.tangent1);
        // With no normal impulse the friction limit is zero, so a point holding no friction keeps none.
        if (p.normalImpulse == 0.0f && old == 0.0f) { continue; }
        changed = true;
        const float limit = c.friction * p.normalImpulse;
        const float t = minf(maxf(old - x.tangentMass * st.speed(c.tangent1, x.tangentA, x.tangentB), -limit), limit);
        p.tangentImpulse = c.tangent1 * t;
        st.push(c.tangent1, x.tangentA, x.tangentB, t - old);
    }

    if (c.rollingResistance > 0.0f) {
        changed = true;
        float total = dot(c.rollingImpulse, c.plane);
        const float spin = st.wb - st.wa;
        const float k = st.ia + st.ib;
        if (spin != 0.0f && k > 0.0f) {
            const float lambda = -spin / k;
            total += lambda;
            st.wa -= st.ia * lambda;
            st.wb += st.ib * lambda;
        }
        float normal = 0.0f;
        for (int i = 0; i < c.count; ++i) { normal += c.points[i].normalImpulse; }
        const float limit = c.rollingResistance * c.rollingRadius * normal;
        const float held = minf(maxf(total, -limit), limit);
        st.wa -= st.ia * (held - total);
        st.wb += st.ib * (held - total);
        c.rollingImpulse = c.plane * held;
    }
    if (changed) { st.store(a, b, c.plane, moveA, moveB); }
}

// Pushes only the bodies that move, since an anchored body takes no impulse.
void pushMoving(Body& a, Body& b, const Axis& x, float impulse, bool moveA, bool moveB) {
    if (moveA) {
        a.v -= x.linearA * impulse;
        a.w -= x.angularA * impulse;
    }
    if (moveB) {
        b.v += x.linearB * impulse;
        b.w += x.angularB * impulse;
    }
}

// One step along the relative spin, so a body locked to a plane only takes impulse about the axes it can
// turn on. The total is then held to the rolling limit.
void resistTurning(Body& a, Body& b, Contact& c) {
    const Vec3 spin = b.w - a.w;
    const float speed = length(spin);
    Vec3 total = c.rollingImpulse;
    if (speed > 0.0f) {
        const Vec3 d = spin * (1.0f / speed);
        const Vec3 ia = a.invInertiaWorld * d, ib = b.invInertiaWorld * d;
        const float k = dot(d, ia) + dot(d, ib);
        if (k > 0.0f) {
            const float lambda = -speed / k;
            total += d * lambda;
            a.w -= ia * lambda;
            b.w += ib * lambda;
        }
    }
    float normal = 0.0f;
    for (int i = 0; i < c.count; ++i) { normal += c.points[i].normalImpulse; }
    const float limit = c.rollingResistance * c.rollingRadius * normal;
    const float mag = length(total);
    if (mag > limit) {
        const Vec3 fix = total * (limit / mag) - total;
        a.w -= a.invInertiaWorld * fix;
        b.w += b.invInertiaWorld * fix;
        total += fix;
    }
    c.rollingImpulse = total;
}

}

Softness makeSoft(float hertz, float zeta, float h) {
    if (hertz == 0.0f) { return { 0.0f, 1.0f, 0.0f }; }
    const float omega = 2.0f * 3.14159265358979f * hertz;
    const float a1 = 2.0f * zeta + h * omega;
    const float a2 = h * omega * a1;
    const float a3 = 1.0f / (1.0f + a2);
    return { omega / a1, a2 * a3, a3 };
}

void prepareContact(const Body& a, const Body& b, Contact& c, PlanarAxes* planar) {
    c.planar = planarContact(a, b, c.normal, c.plane);
    if (c.planar) {
        c.normal = normalize(c.normal - c.plane * dot(c.normal, c.plane));
        c.tangent1 = cross(c.plane, c.normal);
        c.tangent2 = c.plane;
        makePlanarAxes(a, b, c, planar);
    }
    const PlanarState st(a, b, c.plane);
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        p.adjustedSeparation = p.separation - dot(p.anchorB - p.anchorA, c.normal);
        p.relativeVelocity = c.planar ? st.speed(c.normal, planar[i].normalA, planar[i].normalB)
                                      : dot(c.normal, b.v - a.v) + dot(cross(p.anchorB, c.normal), b.w) -
                                            dot(cross(p.anchorA, c.normal), a.w);
        p.maxNormalImpulse = 0.0f;
    }
    if (c.rollingResistance > 0.0f && c.count > 0) {
        float ha = 0.0f, hb = 0.0f;
        for (int i = 0; i < c.count; ++i) {
            ha += std::fabs(dot(c.points[i].anchorA, c.normal));
            hb += std::fabs(dot(c.points[i].anchorB, c.normal));
        }
        if (a.type != MotionType::Dynamic) { ha = 0.0f; }
        if (b.type != MotionType::Dynamic) { hb = 0.0f; }
        c.rollingRadius = maxf(ha, hb) / (float)c.count;
    }
}

void rebaseContact(const Body& a, const Body& b, Contact& c, PointAxes* axes, PlanarAxes* planar) {
    const float shift = dot(b.dp - a.dp, c.normal);
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        p.anchorA = a.dq.rotate(p.anchorA);
        p.anchorB = b.dq.rotate(p.anchorB);
        p.adjustedSeparation += shift;
    }
    if (c.planar) {
        makePlanarAxes(a, b, c, planar);
    } else if (!c.laned) {
        makeAxes(a, b, c, axes);
    }
}

void prepareContacts(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        Contact& c = s.contacts[ci];
        prepareContact(s.bodies[c.bodyA], s.bodies[c.bodyB], c, &s.planarAxes[4 * ci]);
    }
}

void makeContactAxes(const Body& a, const Body& b, const Contact& c, PointAxes* axes) {
    if (!c.planar && !c.laned) { makeAxes(a, b, c, axes); }
}

void warmStart(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        const Contact& c = s.contacts[ci];
        if (c.laned) { continue; }
        const PointAxes* axes = &s.axes[4 * ci];
        Body& a = s.bodies[c.bodyA];
        Body& b = s.bodies[c.bodyB];
        if (c.planar) {
            bool held = c.rollingResistance > 0.0f;
            for (int i = 0; i < c.count; ++i) {
                held = held || c.points[i].normalImpulse != 0.0f || lengthSquared(c.points[i].tangentImpulse) != 0.0f;
            }
            if (!held) { continue; }
            const PlanarAxes* x = &s.planarAxes[4 * ci];
            PlanarState st(a, b, c.plane);
            for (int i = 0; i < c.count; ++i) {
                const ContactPoint& p = c.points[i];
                st.push(c.normal, x[i].normalA, x[i].normalB, p.normalImpulse);
                st.push(c.tangent1, x[i].tangentA, x[i].tangentB, dot(p.tangentImpulse, c.tangent1));
            }
            if (c.rollingResistance > 0.0f) {
                const float r = dot(c.rollingImpulse, c.plane);
                st.wa -= st.ia * r;
                st.wb += st.ib * r;
            }
            st.store(a, b, c.plane, !anchored(a), !anchored(b));
            continue;
        }
        for (int i = 0; i < c.count; ++i) {
            const ContactPoint& p = c.points[i];
            push(a, b, axes[i].normal, p.normalImpulse);
            push(a, b, axes[i].tangent1, dot(p.tangentImpulse, c.tangent1));
            push(a, b, axes[i].tangent2, dot(p.tangentImpulse, c.tangent2));
        }
        if (c.rollingResistance > 0.0f) {
            a.w -= a.invInertiaWorld * c.rollingImpulse;
            b.w += b.invInertiaWorld * c.rollingImpulse;
        }
    }
}

void solveContacts(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx, bool useBias) {
    for (uint32_t ci : active) {
        Contact& c = s.contacts[ci];
        const PointAxes* axes = &s.axes[4 * ci];
        Body& a = s.bodies[c.bodyA];
        Body& b = s.bodies[c.bodyB];
        if (c.planar) {
            solvePlanar(a, b, c, &s.planarAxes[4 * ci], ctx, useBias);
            continue;
        }
        const bool moveA = !anchored(a), moveB = !anchored(b);
        const Softness soft = (!moveA || !moveB) ? ctx.staticSoft : ctx.contactSoft;
        const Vec3 dp = b.dp - a.dp;

        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            const Vec3 d = dp + turned(b, p.anchorB) - turned(a, p.anchorA);
            const float sep = dot(d, c.normal) + p.adjustedSeparation;

            float bias = 0.0f, massScale = 1.0f, impulseScale = 0.0f;
            if (sep > 0.0f) {
                bias = sep * ctx.invH;
            } else if (useBias) {
                bias = maxf(soft.biasRate * sep, -ctx.maxPushout);
                massScale = soft.massScale;
                impulseScale = soft.impulseScale;
            }

            const Axis& x = axes[i].normal;
            const float vn = speedAlong(a, b, c.normal, x);
            float impulse = -x.mass * massScale * (vn + bias) - impulseScale * p.normalImpulse;
            const float total = maxf(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = maxf(p.maxNormalImpulse, impulse);
            pushMoving(a, b, x, impulse, moveA, moveB);
        }

        // The relax pass takes friction points in reverse. One fixed order keeps stacks wobbling.
        for (int j = 0; j < c.count; ++j) {
            const int i = useBias ? j : c.count - 1 - j;
            ContactPoint& p = c.points[i];
            const Axis& x1 = axes[i].tangent1;
            const Axis& x2 = axes[i].tangent2;
            const float old1 = dot(p.tangentImpulse, c.tangent1), old2 = dot(p.tangentImpulse, c.tangent2);
            float t1 = old1 - x1.mass * speedAlong(a, b, c.tangent1, x1);
            float t2 = old2 - x2.mass * speedAlong(a, b, c.tangent2, x2);
            const float limit = c.friction * p.normalImpulse;
            const float mag = std::sqrt(t1 * t1 + t2 * t2);
            if (mag > limit) {
                const float scale = mag > 0.0f ? limit / mag : 0.0f;
                t1 *= scale;
                t2 *= scale;
            }
            p.tangentImpulse = c.tangent1 * t1 + c.tangent2 * t2;
            pushMoving(a, b, x1, t1 - old1, moveA, moveB);
            pushMoving(a, b, x2, t2 - old2, moveA, moveB);
        }

        if (c.rollingResistance > 0.0f) { resistTurning(a, b, c); }
    }
}

void applyRestitution(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx) {
    for (uint32_t ci : active) {
        Contact& c = s.contacts[ci];
        if (c.restitution == 0.0f || c.laned) { continue; }
        const PointAxes* axes = &s.axes[4 * ci];
        Body& a = s.bodies[c.bodyA];
        Body& b = s.bodies[c.bodyB];
        if (c.planar) {
            const PlanarAxes* x = &s.planarAxes[4 * ci];
            PlanarState st(a, b, c.plane);
            for (int i = 0; i < c.count; ++i) {
                ContactPoint& p = c.points[i];
                if (p.relativeVelocity > -ctx.restitutionThreshold || p.maxNormalImpulse == 0.0f) { continue; }
                const float vn = st.speed(c.normal, x[i].normalA, x[i].normalB);
                float impulse = -x[i].normalMass * (vn + c.restitution * p.relativeVelocity);
                const float total = maxf(p.normalImpulse + impulse, 0.0f);
                impulse = total - p.normalImpulse;
                p.normalImpulse = total;
                p.maxNormalImpulse = maxf(p.maxNormalImpulse, impulse);
                st.push(c.normal, x[i].normalA, x[i].normalB, impulse);
            }
            st.store(a, b, c.plane, !anchored(a), !anchored(b));
            continue;
        }
        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            if (p.relativeVelocity > -ctx.restitutionThreshold || p.maxNormalImpulse == 0.0f) { continue; }
            const Axis& x = axes[i].normal;
            const float vn = speedAlong(a, b, c.normal, x);
            float impulse = -x.mass * (vn + c.restitution * p.relativeVelocity);
            const float total = maxf(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = maxf(p.maxNormalImpulse, impulse);
            push(a, b, x, impulse);
        }
    }
}

}
