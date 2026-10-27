#include "internal.h"

namespace twofold {

namespace {

Vec3 velocityAt(const Body& b, const Vec3& r) { return b.v + cross(b.w, r); }

void applyImpulse(Body& b, const Vec3& p, const Vec3& r) {
    b.v += b.invMassWorld * p;
    b.w += b.invInertiaWorld * cross(r, p);
}

float effectiveMass(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& d) {
    const Vec3 rna = cross(ra, d), rnb = cross(rb, d);
    const float k = dot(d, a.invMassWorld * d) + dot(d, b.invMassWorld * d)
                  + dot(rna, a.invInertiaWorld * rna) + dot(rnb, b.invInertiaWorld * rnb);
    // A direction the bodies are locked against reaches k only through rounding.
    return k > 1e-6f * (a.invMass + b.invMass) ? 1.0f / k : 0.0f;
}

bool anchored(const Body& b) { return b.type != MotionType::Dynamic || !b.awake; }

}

Softness makeSoft(float hertz, float zeta, float h) {
    if (hertz == 0.0f) { return { 0.0f, 1.0f, 0.0f }; }
    const float omega = 2.0f * 3.14159265358979f * hertz;
    const float a1 = 2.0f * zeta + h * omega;
    const float a2 = h * omega * a1;
    const float a3 = 1.0f / (1.0f + a2);
    return { omega / a1, a2 * a3, a3 };
}

void prepareContact(const Body& a, const Body& b, Contact& c) {
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        const Vec3 ra = p.anchorA, rb = p.anchorB;
        p.adjustedSeparation = p.separation - dot(rb - ra, c.normal);
        p.normalMass = effectiveMass(a, b, ra, rb, c.normal);
        p.tangentMass1 = effectiveMass(a, b, ra, rb, c.tangent1);
        p.tangentMass2 = effectiveMass(a, b, ra, rb, c.tangent2);
        p.relativeVelocity = dot(c.normal, velocityAt(b, rb) - velocityAt(a, ra));
        p.maxNormalImpulse = 0.0f;
    }
}

void rebaseContact(const Body& a, const Body& b, Contact& c) {
    const float shift = dot(b.dp - a.dp, c.normal);
    for (int i = 0; i < c.count; ++i) {
        ContactPoint& p = c.points[i];
        p.anchorA = a.dq.rotate(p.anchorA);
        p.anchorB = b.dq.rotate(p.anchorB);
        p.adjustedSeparation += shift;
        p.normalMass = effectiveMass(a, b, p.anchorA, p.anchorB, c.normal);
        p.tangentMass1 = effectiveMass(a, b, p.anchorA, p.anchorB, c.tangent1);
        p.tangentMass2 = effectiveMass(a, b, p.anchorA, p.anchorB, c.tangent2);
    }
}

void prepareContacts(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        Contact& c = contacts[ci];
        prepareContact(bodies[c.bodyA], bodies[c.bodyB], c);
    }
}

void warmStart(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active) {
    for (uint32_t ci : active) {
        Contact& c = contacts[ci];
        Body& a = bodies[c.bodyA];
        Body& b = bodies[c.bodyB];
        for (int i = 0; i < c.count; ++i) {
            const ContactPoint& p = c.points[i];
            const Vec3 impulse = c.normal * p.normalImpulse + p.tangentImpulse;
            applyImpulse(a, -impulse, p.anchorA);
            applyImpulse(b, impulse, p.anchorB);
        }
    }
}

void solveContacts(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active,
                   const SolverContext& ctx, bool useBias) {
    for (uint32_t ci : active) {
        Contact& c = contacts[ci];
        Body& a = bodies[c.bodyA];
        Body& b = bodies[c.bodyB];
        const Softness soft = (anchored(a) || anchored(b)) ? ctx.staticSoft : ctx.contactSoft;
        const Vec3 dp = b.dp - a.dp;

        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            const Vec3 d = dp + b.dq.rotate(p.anchorB) - a.dq.rotate(p.anchorA);
            const float s = dot(d, c.normal) + p.adjustedSeparation;

            float bias = 0.0f, massScale = 1.0f, impulseScale = 0.0f;
            if (s > 0.0f) {
                bias = s * ctx.invH;
            } else if (useBias) {
                bias = std::fmax(soft.biasRate * s, -ctx.maxPushout);
                massScale = soft.massScale;
                impulseScale = soft.impulseScale;
            }

            const float vn = dot(c.normal, velocityAt(b, p.anchorB) - velocityAt(a, p.anchorA));
            float impulse = -p.normalMass * massScale * (vn + bias) - impulseScale * p.normalImpulse;
            const float total = std::fmax(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = std::fmax(p.maxNormalImpulse, impulse);

            const Vec3 pv = c.normal * impulse;
            applyImpulse(a, -pv, p.anchorA);
            applyImpulse(b, pv, p.anchorB);
        }

        // The relax pass takes friction points in reverse. One fixed order keeps stacks wobbling.
        for (int j = 0; j < c.count; ++j) {
            ContactPoint& p = c.points[useBias ? j : c.count - 1 - j];
            const Vec3 dv = velocityAt(b, p.anchorB) - velocityAt(a, p.anchorA);
            const Vec3 delta = c.tangent1 * (-p.tangentMass1 * dot(dv, c.tangent1))
                             + c.tangent2 * (-p.tangentMass2 * dot(dv, c.tangent2));
            Vec3 total = p.tangentImpulse + delta;
            const float limit = c.friction * p.normalImpulse;
            const float mag = length(total);
            if (mag > limit) { total = mag > 0.0f ? total * (limit / mag) : Vec3{}; }
            const Vec3 applied = total - p.tangentImpulse;
            p.tangentImpulse = total;
            applyImpulse(a, -applied, p.anchorA);
            applyImpulse(b, applied, p.anchorB);
        }
    }
}

void applyRestitution(std::vector<Body>& bodies, std::vector<Contact>& contacts, const std::vector<uint32_t>& active,
                      const SolverContext& ctx) {
    for (uint32_t ci : active) {
        Contact& c = contacts[ci];
        if (c.restitution == 0.0f) { continue; }
        Body& a = bodies[c.bodyA];
        Body& b = bodies[c.bodyB];
        for (int i = 0; i < c.count; ++i) {
            ContactPoint& p = c.points[i];
            if (p.relativeVelocity > -ctx.restitutionThreshold || p.maxNormalImpulse == 0.0f) { continue; }
            const float vn = dot(c.normal, velocityAt(b, p.anchorB) - velocityAt(a, p.anchorA));
            float impulse = -p.normalMass * (vn + c.restitution * p.relativeVelocity);
            const float total = std::fmax(p.normalImpulse + impulse, 0.0f);
            impulse = total - p.normalImpulse;
            p.normalImpulse = total;
            p.maxNormalImpulse = std::fmax(p.maxNormalImpulse, impulse);
            const Vec3 pv = c.normal * impulse;
            applyImpulse(a, -pv, p.anchorA);
            applyImpulse(b, pv, p.anchorB);
        }
    }
}

}
