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

void prepareContact(const Body& a, const Body& b, Contact& c, PointAxes* axes) {
    makeAxes(a, b, c, axes);
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        p.adjustedSeparation = p.separation - dot(p.anchorB - p.anchorA, c.normal);
        p.relativeVelocity = speedAlong(a, b, c.normal, axes[i].normal);
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
        c.rollingRadius = std::fmax(ha, hb) / (float)c.count;
    }
}

void rebaseContact(const Body& a, const Body& b, Contact& c, PointAxes* axes) {
    const float shift = dot(b.dp - a.dp, c.normal);
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        p.anchorA = a.dq.rotate(p.anchorA);
        p.anchorB = b.dq.rotate(p.anchorB);
        p.adjustedSeparation += shift;
    }
    makeAxes(a, b, c, axes);
}

void updateResponses(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        const Contact& c = s.contacts[ci];
        const Body& a = s.bodies[c.bodyA];
        const Body& b = s.bodies[c.bodyB];
        PointAxes* axes = &s.axes[4 * ci];
        for (int i = 0; i < c.count; ++i) {
            for (Axis* x : { &axes[i].normal, &axes[i].tangent1, &axes[i].tangent2 }) {
                x->angularA = a.invInertiaWorld * x->armA;
                x->angularB = b.invInertiaWorld * x->armB;
            }
        }
    }
}

void prepareContacts(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        Contact& c = s.contacts[ci];
        prepareContact(s.bodies[c.bodyA], s.bodies[c.bodyB], c, &s.axes[4 * ci]);
    }
}

void warmStart(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        const Contact& c = s.contacts[ci];
        const PointAxes* axes = &s.axes[4 * ci];
        Body& a = s.bodies[c.bodyA];
        Body& b = s.bodies[c.bodyB];
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
        const Softness soft = (anchored(a) || anchored(b)) ? ctx.staticSoft : ctx.contactSoft;
        const Vec3 dp = b.dp - a.dp;

        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            const Vec3 d = dp + b.dq.rotate(p.anchorB) - a.dq.rotate(p.anchorA);
            const float sep = dot(d, c.normal) + p.adjustedSeparation;

            float bias = 0.0f, massScale = 1.0f, impulseScale = 0.0f;
            if (sep > 0.0f) {
                bias = sep * ctx.invH;
            } else if (useBias) {
                bias = std::fmax(soft.biasRate * sep, -ctx.maxPushout);
                massScale = soft.massScale;
                impulseScale = soft.impulseScale;
            }

            const Axis& x = axes[i].normal;
            const float vn = speedAlong(a, b, c.normal, x);
            float impulse = -x.mass * massScale * (vn + bias) - impulseScale * p.normalImpulse;
            const float total = std::fmax(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = std::fmax(p.maxNormalImpulse, impulse);
            push(a, b, x, impulse);
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
            push(a, b, x1, t1 - old1);
            push(a, b, x2, t2 - old2);
        }

        if (c.rollingResistance > 0.0f) { resistTurning(a, b, c); }
    }
}

void applyRestitution(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx) {
    for (uint32_t ci : active) {
        Contact& c = s.contacts[ci];
        if (c.restitution == 0.0f) { continue; }
        const PointAxes* axes = &s.axes[4 * ci];
        Body& a = s.bodies[c.bodyA];
        Body& b = s.bodies[c.bodyB];
        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            if (p.relativeVelocity > -ctx.restitutionThreshold || p.maxNormalImpulse == 0.0f) { continue; }
            const Axis& x = axes[i].normal;
            const float vn = speedAlong(a, b, c.normal, x);
            float impulse = -x.mass * (vn + c.restitution * p.relativeVelocity);
            const float total = std::fmax(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = std::fmax(p.maxNormalImpulse, impulse);
            push(a, b, x, impulse);
        }
    }
}

}
