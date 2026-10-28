#include <cmath>

#include "internal.h"

namespace twofold {

namespace {

constexpr float kPi = 3.14159265358979f;

Mat3 skew(const Vec3& v) {
    Mat3 m;
    m.r[0] = { 0.0f, -v.z, v.y };
    m.r[1] = { v.z, 0.0f, -v.x };
    m.r[2] = { -v.y, v.x, 0.0f };
    return m;
}

float trace(const Mat3& m) { return m.r[0].x + m.r[1].y + m.r[2].z; }

// Inverse within the directions m acts in. A body locked against a direction leaves m without rank there.
Mat3 pseudoInverse(const Mat3& m) {
    const float tr = trace(m);
    if (tr <= 0.0f) { return Mat3::zero(); }
    const float third = tr / 3.0f;
    if (dot(m.r[0], cross(m.r[1], m.r[2])) > 1e-3f * third * third * third) { return inverse(m); }

    float a[3][3], v[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            a[i][k] = 0.5f * (m.r[i][k] + m.r[k][i]);
            v[i][k] = i == k ? 1.0f : 0.0f;
        }
    }
    // Jacobi rotations until the off-diagonal is gone.
    for (int sweep = 0; sweep < 8; ++sweep) {
        if (a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2] <= 1e-14f * tr * tr) { break; }
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (a[p][q] == 0.0f) { continue; }
                const float theta = (a[q][q] - a[p][p]) / (2.0f * a[p][q]);
                const float t = (theta >= 0.0f ? 1.0f : -1.0f) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0f));
                const float c = 1.0f / std::sqrt(t * t + 1.0f), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const float kp = a[k][p], kq = a[k][q];
                    a[k][p] = c * kp - s * kq;
                    a[k][q] = s * kp + c * kq;
                }
                for (int k = 0; k < 3; ++k) {
                    const float pk = a[p][k], qk = a[q][k];
                    a[p][k] = c * pk - s * qk;
                    a[q][k] = s * pk + c * qk;
                }
                for (int k = 0; k < 3; ++k) {
                    const float kp = v[k][p], kq = v[k][q];
                    v[k][p] = c * kp - s * kq;
                    v[k][q] = s * kp + c * kq;
                }
            }
        }
    }
    const float top = maxf(a[0][0], maxf(a[1][1], a[2][2]));
    Mat3 out = Mat3::zero();
    for (int i = 0; i < 3; ++i) {
        if (a[i][i] <= 1e-6f * top) { continue; }
        const Vec3 e{ v[0][i], v[1][i], v[2][i] };
        out = out + Mat3::outer(e, e) * (1.0f / a[i][i]);
    }
    return out;
}

Body& bodyOf(SolverData& s, uint32_t i) { return i == kWorldBody ? s.world : s.bodies[i]; }

Vec3 arm(const Body& b, const Vec3& local) { return b.xf.q.rotate(local - b.localCenter); }

Axis makeTurn(const Body& a, const Body& b, const Vec3& d) {
    Axis x;
    x.armA = x.armB = d;
    x.angularA = a.invInertiaWorld * d;
    x.angularB = b.invInertiaWorld * d;
    const float k = dot(d, x.angularA) + dot(d, x.angularB);
    x.mass = k > 1e-6f * (trace(a.invInertiaWorld) + trace(b.invInertiaWorld)) ? 1.0f / k : 0.0f;
    return x;
}

Axis flipped(const Axis& x) {
    Axis f = x;
    f.armA = -x.armA;
    f.armB = -x.armB;
    f.linearA = -x.linearA;
    f.angularA = -x.angularA;
    f.linearB = -x.linearB;
    f.angularB = -x.angularB;
    return f;
}

// How a point offset r from the center moves per unit impulse there.
Mat3 pointResponse(const Body& b, const Vec3& r) {
    const Mat3 s = skew(r);
    return b.invMassWorld - s * b.invInertiaWorld * s;
}

// Turns about the axis through the anchors. pointImpulse is the point impulse per unit turn that holds them.
Axis anchoredTurn(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& axis, Vec3& pointImpulse) {
    const Vec3 drift = cross(b.invInertiaWorld * axis, rb) + cross(a.invInertiaWorld * axis, ra);
    pointImpulse = pseudoInverse(pointResponse(a, ra) + pointResponse(b, rb)) * -drift;
    Axis x;
    x.armA = axis + cross(ra, pointImpulse);
    x.armB = axis + cross(rb, pointImpulse);
    x.linearA = a.invMassWorld * pointImpulse;
    x.angularA = a.invInertiaWorld * x.armA;
    x.linearB = b.invMassWorld * pointImpulse;
    x.angularB = b.invInertiaWorld * x.armB;
    const float k = dot(pointImpulse, x.linearA + x.linearB) + dot(x.armA, x.angularA) + dot(x.armB, x.angularB);
    x.mass = k > 1e-6f * (trace(a.invInertiaWorld) + trace(b.invInertiaWorld)) ? 1.0f / k : 0.0f;
    return x;
}

// Slides along the axis with the turn impulse that keeps the relative rotation fixed.
Axis lockedSlide(const Body& a, const Body& b, const Vec3& ra, const Vec3& rb, const Vec3& axis) {
    Axis x = makeAxis(a, b, ra, rb, axis);
    const Vec3 twist = pseudoInverse(a.invInertiaWorld + b.invInertiaWorld) * -(x.angularA + x.angularB);
    x.armA += twist;
    x.armB += twist;
    x.angularA = a.invInertiaWorld * x.armA;
    x.angularB = b.invInertiaWorld * x.armB;
    const float k = dot(axis, x.linearA + x.linearB) + dot(x.armA, x.angularA) + dot(x.armB, x.angularB);
    x.mass = k > 1e-6f * (a.invMass + b.invMass) ? 1.0f / k : 0.0f;
    return x;
}

void pushPoint(Body& a, Body& b, const Vec3& ra, const Vec3& rb, const Vec3& impulse) {
    a.v -= a.invMassWorld * impulse;
    a.w -= a.invInertiaWorld * cross(ra, impulse);
    b.v += b.invMassWorld * impulse;
    b.w += b.invInertiaWorld * cross(rb, impulse);
}

void pushTurn(Body& a, Body& b, const Vec3& impulse) {
    a.w -= a.invInertiaWorld * impulse;
    b.w += b.invInertiaWorld * impulse;
}

Softness softFor(const SolverContext& ctx, bool useBias) { return useBias ? ctx.jointSoft : Softness{ 0.0f, 1.0f, 0.0f }; }

// Holds the points together along the directions proj keeps.
void solvePoint(Body& a, Body& b, const Vec3& ra, const Vec3& rb, const Mat3& proj, const Vec3& error, Vec3& total,
                const Softness& soft) {
    const Mat3 k = proj * (pointResponse(a, ra) + pointResponse(b, rb)) * proj;
    const Vec3 speed = proj * (b.v + cross(b.w, rb) - a.v - cross(a.w, ra));
    const Vec3 old = proj * total;
    const Vec3 impulse = pseudoInverse(k) * (speed + error * soft.biasRate) * -soft.massScale - old * soft.impulseScale;
    total = old + impulse;
    pushPoint(a, b, ra, rb, impulse);
}

// Stops the bodies turning relative to each other about the directions proj keeps.
void solveTurn(Body& a, Body& b, const Mat3& proj, const Vec3& error, Vec3& total, const Softness& soft) {
    const Mat3 k = proj * (a.invInertiaWorld + b.invInertiaWorld) * proj;
    const Vec3 speed = proj * (b.w - a.w);
    const Vec3 old = proj * total;
    const Vec3 impulse = pseudoInverse(k) * (speed + error * soft.biasRate) * -soft.massScale - old * soft.impulseScale;
    total = old + impulse;
    pushTurn(a, b, impulse);
}

void solveEqual(Body& a, Body& b, const Axis& x, const Vec3& d, float error, float& total, const Softness& soft) {
    const float impulse =
        -x.mass * soft.massScale * (speedAlong(a, b, d, x) + soft.biasRate * error) - soft.impulseScale * total;
    total += impulse;
    push(a, b, x, impulse);
}

// Keeps error at or above zero. While it is above zero, the bodies may close only that far this substep.
void solveLimit(Body& a, Body& b, const Axis& x, const Vec3& d, float error, float& total, const SolverContext& ctx,
                bool useBias) {
    float bias = 0.0f, massScale = 1.0f, impulseScale = 0.0f;
    if (error > 0.0f) {
        bias = error * ctx.invH;
    } else if (useBias) {
        bias = ctx.jointSoft.biasRate * error;
        massScale = ctx.jointSoft.massScale;
        impulseScale = ctx.jointSoft.impulseScale;
    }
    float impulse = -x.mass * massScale * (speedAlong(a, b, d, x) + bias) - impulseScale * total;
    const float next = maxf(total + impulse, 0.0f);
    impulse = next - total;
    total = next;
    push(a, b, x, impulse);
}

void solveMotor(Body& a, Body& b, const Axis& x, const Vec3& d, float speed, float maxImpulse, float& total) {
    const float old = total;
    total = minf(maxf(old - x.mass * (speedAlong(a, b, d, x) - speed), -maxImpulse), maxImpulse);
    push(a, b, x, total - old);
}

Vec3 turnError(const Joint& j, const Body& a, const Body& b) {
    Quat e = b.xf.q * (a.xf.q * j.restTurn).conjugate();
    if (e.w < 0.0f) { e = -e; }
    return e.vec() * 2.0f;
}

float hingeAngle(const Joint& j, const Body& a, const Body& b) {
    const Vec3 axis = a.xf.q.rotate(j.axisA);
    const Vec3 ra = a.xf.q.rotate(j.refA);
    const Vec3 rb = b.xf.q.rotate(j.refB);
    return std::atan2(dot(cross(ra, rb), axis), dot(ra, rb));
}

Mat3 across(const Vec3& axis) { return Mat3{} - Mat3::outer(axis, axis); }

// Holds the relative turn of the bodies with an angular impulse of up to maxImpulse.
void solveTurnFriction(Body& a, Body& b, float maxImpulse, Vec3& total) {
    const Vec3 old = total;
    Vec3 next = old - pseudoInverse(a.invInertiaWorld + b.invInertiaWorld) * (b.w - a.w);
    const float len = length(next);
    if (len > maxImpulse) { next = next * (maxImpulse / len); }
    total = next;
    pushTurn(a, b, next - old);
}

// Body B's turn in body A's joint frame, split into a swing that tilts the x axis and a twist about it.
struct SwingTwist {
    Quat frame;
    Vec3 axisA, axisB;
    Quat swing;
    float twist = 0.0f;
};

SwingTwist swingTwist(const Joint& j, const Body& a, const Body& b) {
    SwingTwist st;
    st.frame = a.xf.q * j.frameA;
    const Quat frameB = b.xf.q * j.frameB;
    st.axisA = st.frame.rotate({ 1, 0, 0 });
    st.axisB = frameB.rotate({ 1, 0, 0 });
    Quat r = st.frame.conjugate() * frameB;
    if (r.w < 0.0f) { r = -r; }
    const float s = std::sqrt(r.x * r.x + r.w * r.w);
    const Quat t = s > 0.0f ? Quat{ r.x / s, 0.0f, 0.0f, r.w / s } : Quat{};
    st.swing = r * t.conjugate();
    st.twist = 2.0f * std::atan2(t.x, t.w);
    return st;
}

// The twist changes at the relative spin along rate. False near a half turn of swing, where it is undefined.
bool twistRate(const SwingTwist& st, Vec3& rate) {
    const float c = 1.0f + dot(st.axisA, st.axisB);
    if (c < 1e-3f) { return false; }
    rate = (st.axisA + st.axisB) / c;
    return true;
}

bool swingFree(const Joint& j) { return j.swingY >= 1.0f && j.swingZ >= 1.0f; }

// How much further the swing may grow before it meets the cone, in radians, and the relative spin along rate
// that changes this. False while the swing is too small to have a direction.
bool swingRoom(const Joint& j, const SwingTwist& st, float& error, Vec3& rate) {
    const float len = std::sqrt(st.swing.y * st.swing.y + st.swing.z * st.swing.z);
    if (len < 1e-6f) { return false; }
    const float cy = st.swing.y / len, cz = st.swing.z / len;
    // A turn about y tilts the axis toward z, so swingZ bounds it.
    const float ey = cy / j.swingZ, ez = cz / j.swingY;
    const float reach = 1.0f / std::sqrt(ey * ey + ez * ez);
    const Vec3 k = st.frame.rotate({ 0.0f, cy, cz });
    rate = -k;
    float allowed = kPi;
    if (reach < 1.0f) {
        allowed = 2.0f * std::asin(reach);
        // Spinning about axis A turns the swing direction, and with it how far the cone reaches.
        const float slope = 2.0f * reach * reach * reach * (ey * cz / j.swingZ - ez * cy / j.swingY) /
                            std::sqrt(1.0f - reach * reach);
        const Vec3 side = cross(st.axisB, k);
        const float det = dot(st.axisA, side);
        if (std::fabs(det) > 1e-6f) { rate += side * (slope / det); }
    }
    error = allowed - 2.0f * std::atan2(len, st.swing.w);
    return true;
}

void warmStartJoint(Body& a, Body& b, const Joint& j) {
    const Vec3 ra = arm(a, j.localA), rb = arm(b, j.localB);
    const Vec3 d = (b.com + rb) - (a.com + ra);
    float axial = j.axialImpulse + j.motorImpulse + j.lowerImpulse - j.upperImpulse;
    if (j.maxFriction > 0.0f) { axial += j.frictionImpulse; }
    switch (j.type) {
        case JointType::Point:
            pushPoint(a, b, ra, rb, j.linearImpulse);
            break;
        case JointType::Fixed:
            pushPoint(a, b, ra, rb, j.linearImpulse);
            pushTurn(a, b, j.angularImpulse);
            break;
        case JointType::Hinge: {
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            Vec3 along;
            pushPoint(a, b, ra, rb, j.linearImpulse);
            pushTurn(a, b, across(axis) * j.angularImpulse);
            push(a, b, anchoredTurn(a, b, ra, rb, axis, along), axial);
            break;
        }
        case JointType::Slider: {
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            pushPoint(a, b, ra + d, rb, across(axis) * j.linearImpulse);
            pushTurn(a, b, j.angularImpulse);
            push(a, b, lockedSlide(a, b, ra + d, rb, axis), axial);
            break;
        }
        case JointType::Distance: {
            const float len = length(d);
            if (len > 0.0f) { push(a, b, makeAxis(a, b, ra, rb, d / len), axial); }
            break;
        }
        case JointType::SwingTwist: {
            pushPoint(a, b, ra, rb, j.linearImpulse);
            pushTurn(a, b, j.frictionTurn);
            const SwingTwist st = swingTwist(j, a, b);
            Vec3 rate, along;
            if (axial != 0.0f && twistRate(st, rate)) { push(a, b, anchoredTurn(a, b, ra, rb, rate, along), axial); }
            float error = 0.0f;
            if (j.swingImpulse != 0.0f && swingRoom(j, st, error, rate)) {
                push(a, b, anchoredTurn(a, b, ra, rb, rate, along), j.swingImpulse);
            }
            break;
        }
    }
}

void solveJoint(Body& a, Body& b, Joint& j, const SolverContext& ctx, bool useBias) {
    const Softness soft = softFor(ctx, useBias);
    const Vec3 ra = arm(a, j.localA), rb = arm(b, j.localB);
    const Vec3 d = (b.com + rb) - (a.com + ra);
    switch (j.type) {
        case JointType::Point:
            solvePoint(a, b, ra, rb, Mat3{}, d, j.linearImpulse, soft);
            break;
        case JointType::Fixed:
            solveTurn(a, b, Mat3{}, turnError(j, a, b), j.angularImpulse, soft);
            solvePoint(a, b, ra, rb, Mat3{}, d, j.linearImpulse, soft);
            break;
        case JointType::Hinge: {
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            Vec3 along;
            const Axis turn = anchoredTurn(a, b, ra, rb, axis, along);
            j.rowForce[0] = along;
            j.rowTorque[0] = axis;
            if (j.maxFriction > 0.0f) { solveMotor(a, b, turn, along, 0.0f, j.maxFriction * ctx.h, j.frictionImpulse); }
            if (j.maxMotor > 0.0f) { solveMotor(a, b, turn, along, j.motorSpeed, j.maxMotor * ctx.h, j.motorImpulse); }
            if (j.limited) {
                const float angle = hingeAngle(j, a, b);
                solveLimit(a, b, turn, along, angle - j.lower, j.lowerImpulse, ctx, useBias);
                solveLimit(a, b, flipped(turn), -along, j.upper - angle, j.upperImpulse, ctx, useBias);
            }
            solveTurn(a, b, across(axis), cross(axis, b.xf.q.rotate(j.axisB)), j.angularImpulse, soft);
            solvePoint(a, b, ra, rb, Mat3{}, d, j.linearImpulse, soft);
            break;
        }
        case JointType::Slider: {
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            const Axis slide = lockedSlide(a, b, ra + d, rb, axis);
            j.rowForce[0] = axis;
            j.rowTorque[0] = slide.armB - cross(rb, axis);
            if (j.maxFriction > 0.0f) { solveMotor(a, b, slide, axis, 0.0f, j.maxFriction * ctx.h, j.frictionImpulse); }
            if (j.maxMotor > 0.0f) { solveMotor(a, b, slide, axis, j.motorSpeed, j.maxMotor * ctx.h, j.motorImpulse); }
            if (j.limited) {
                const float t = dot(d, axis);
                solveLimit(a, b, slide, axis, t - j.lower, j.lowerImpulse, ctx, useBias);
                solveLimit(a, b, flipped(slide), -axis, j.upper - t, j.upperImpulse, ctx, useBias);
            }
            solveTurn(a, b, Mat3{}, turnError(j, a, b), j.angularImpulse, soft);
            const Mat3 side = across(axis);
            solvePoint(a, b, ra + d, rb, side, side * d, j.linearImpulse, soft);
            break;
        }
        case JointType::Distance: {
            const float len = length(d);
            if (len <= 0.0f) { break; }
            const Vec3 u = d / len;
            const Axis x = makeAxis(a, b, ra, rb, u);
            j.rowForce[0] = u;
            if (j.limited ? j.lower >= j.upper : j.hertz == 0.0f) {
                solveEqual(a, b, x, u, len - (j.limited ? j.lower : j.length), j.axialImpulse, soft);
                break;
            }
            if (j.hertz > 0.0f) {
                const Softness spring = makeSoft(j.hertz, j.dampingRatio, ctx.h);
                solveEqual(a, b, x, u, len - j.length, j.axialImpulse, spring);
            }
            if (j.limited) {
                solveLimit(a, b, x, u, len - j.lower, j.lowerImpulse, ctx, useBias);
                solveLimit(a, b, flipped(x), -u, j.upper - len, j.upperImpulse, ctx, useBias);
            }
            break;
        }
        case JointType::SwingTwist: {
            if (j.maxFriction > 0.0f) { solveTurnFriction(a, b, j.maxFriction * ctx.h, j.frictionTurn); }
            const SwingTwist st = swingTwist(j, a, b);
            Vec3 rate, along;
            if (j.limited && twistRate(st, rate)) {
                const Axis turn = anchoredTurn(a, b, ra, rb, rate, along);
                j.rowForce[0] = along;
                j.rowTorque[0] = rate;
                solveLimit(a, b, turn, along, st.twist - j.lower, j.lowerImpulse, ctx, useBias);
                solveLimit(a, b, flipped(turn), -along, j.upper - st.twist, j.upperImpulse, ctx, useBias);
            } else {
                j.lowerImpulse = j.upperImpulse = 0.0f;
            }
            float error = 0.0f;
            if (!swingFree(j) && swingRoom(j, st, error, rate)) {
                const Axis turn = anchoredTurn(a, b, ra, rb, rate, along);
                j.rowForce[1] = along;
                j.rowTorque[1] = rate;
                solveLimit(a, b, turn, along, error, j.swingImpulse, ctx, useBias);
            } else {
                j.swingImpulse = 0.0f;
            }
            solvePoint(a, b, ra, rb, Mat3{}, d, j.linearImpulse, soft);
            break;
        }
    }
}

void shift(Body& b, const Vec3& move, const Vec3& turn) {
    if (b.type != MotionType::Dynamic || !b.awake) { return; }
    b.com += move;
    b.xf.q = integrate(b.xf.q, turn, 1.0f);
    b.xf.p = b.com - b.xf.q.rotate(b.localCenter);
}

void movePoint(Body& a, Body& b, const Vec3& ra, const Vec3& rb, const Mat3& proj, const Vec3& error) {
    const Mat3 k = proj * (pointResponse(a, ra) + pointResponse(b, rb)) * proj;
    const Vec3 p = pseudoInverse(k) * (proj * error) * -1.0f;
    shift(a, -(a.invMassWorld * p), -(a.invInertiaWorld * cross(ra, p)));
    shift(b, b.invMassWorld * p, b.invInertiaWorld * cross(rb, p));
}

void moveTurn(Body& a, Body& b, const Mat3& proj, const Vec3& error) {
    const Mat3 k = proj * (a.invInertiaWorld + b.invInertiaWorld) * proj;
    const Vec3 l = pseudoInverse(k) * (proj * error) * -1.0f;
    shift(a, {}, -(a.invInertiaWorld * l));
    shift(b, {}, b.invInertiaWorld * l);
}

void moveAlong(Body& a, Body& b, const Axis& x, float error) {
    const float l = -x.mass * error;
    shift(a, -(x.linearA * l), -(x.angularA * l));
    shift(b, x.linearB * l, x.angularB * l);
}

void settleJoint(Body& a, Body& b, const Joint& j) {
    Vec3 ra = arm(a, j.localA), rb = arm(b, j.localB);
    Vec3 d = (b.com + rb) - (a.com + ra);
    auto measure = [&]() {
        ra = arm(a, j.localA);
        rb = arm(b, j.localB);
        d = (b.com + rb) - (a.com + ra);
    };
    switch (j.type) {
        case JointType::Point:
            movePoint(a, b, ra, rb, Mat3{}, d);
            break;
        case JointType::Fixed:
            moveTurn(a, b, Mat3{}, turnError(j, a, b));
            measure();
            movePoint(a, b, ra, rb, Mat3{}, d);
            break;
        case JointType::Hinge: {
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            if (j.limited) {
                const float angle = hingeAngle(j, a, b);
                const Axis turn = makeTurn(a, b, axis);
                if (angle < j.lower) { moveAlong(a, b, turn, angle - j.lower); }
                if (angle > j.upper) { moveAlong(a, b, turn, angle - j.upper); }
            }
            moveTurn(a, b, across(axis), cross(axis, b.xf.q.rotate(j.axisB)));
            measure();
            movePoint(a, b, ra, rb, Mat3{}, d);
            break;
        }
        case JointType::Slider: {
            moveTurn(a, b, Mat3{}, turnError(j, a, b));
            measure();
            const Vec3 axis = a.xf.q.rotate(j.axisA);
            if (j.limited) {
                const float t = dot(d, axis);
                const Axis slide = lockedSlide(a, b, ra + d, rb, axis);
                if (t < j.lower) { moveAlong(a, b, slide, t - j.lower); }
                if (t > j.upper) { moveAlong(a, b, slide, t - j.upper); }
                measure();
            }
            const Mat3 side = across(axis);
            movePoint(a, b, ra + d, rb, side, d);
            break;
        }
        case JointType::Distance: {
            const float len = length(d);
            if (len <= 0.0f) { break; }
            const Axis x = makeAxis(a, b, ra, rb, d / len);
            if (j.limited ? j.lower >= j.upper : j.hertz == 0.0f) {
                moveAlong(a, b, x, len - (j.limited ? j.lower : j.length));
            } else if (j.limited && len < j.lower) {
                moveAlong(a, b, x, len - j.lower);
            } else if (j.limited && len > j.upper) {
                moveAlong(a, b, x, len - j.upper);
            }
            break;
        }
        case JointType::SwingTwist: {
            Vec3 rate;
            if (j.limited) {
                const SwingTwist st = swingTwist(j, a, b);
                if (twistRate(st, rate)) {
                    const Axis turn = makeTurn(a, b, rate);
                    if (st.twist < j.lower) { moveAlong(a, b, turn, st.twist - j.lower); }
                    if (st.twist > j.upper) { moveAlong(a, b, turn, st.twist - j.upper); }
                }
            }
            float error = 0.0f;
            if (!swingFree(j) && swingRoom(j, swingTwist(j, a, b), error, rate) && error < 0.0f) {
                moveAlong(a, b, makeTurn(a, b, rate), error);
            }
            measure();
            movePoint(a, b, ra, rb, Mat3{}, d);
            break;
        }
    }
}

}

void settleJoints(SolverData& s, const std::vector<uint32_t>& active, int iterations) {
    for (int i = 0; i < iterations; ++i) {
        for (uint32_t ji : active) {
            const Joint& j = s.joints[ji];
            settleJoint(bodyOf(s, j.bodyA), bodyOf(s, j.bodyB), j);
        }
    }
}

void warmStartJoints(SolverData& s, const std::vector<uint32_t>& active) {
    for (uint32_t ji : active) {
        const Joint& j = s.joints[ji];
        warmStartJoint(bodyOf(s, j.bodyA), bodyOf(s, j.bodyB), j);
    }
}

void solveJoints(SolverData& s, const std::vector<uint32_t>& active, const SolverContext& ctx, bool useBias) {
    for (uint32_t ji : active) {
        Joint& j = s.joints[ji];
        solveJoint(bodyOf(s, j.bodyA), bodyOf(s, j.bodyB), j, ctx, useBias);
        if (!useBias) { addReaction(j); }
    }
}

float jointValue(const Joint& j, const Body& a, const Body& b) {
    const Vec3 d = (b.com + arm(b, j.localB)) - (a.com + arm(a, j.localA));
    switch (j.type) {
        case JointType::Hinge: return hingeAngle(j, a, b);
        case JointType::Slider: return dot(d, a.xf.q.rotate(j.axisA));
        case JointType::Distance: return length(d);
        case JointType::SwingTwist: return swingTwist(j, a, b).twist;
        default: return 0.0f;
    }
}

void addReaction(Joint& j) {
    Vec3 force = j.linearImpulse, torque;
    const float axial = j.axialImpulse + j.motorImpulse + j.lowerImpulse - j.upperImpulse + j.frictionImpulse;
    switch (j.type) {
        case JointType::Point:
            break;
        case JointType::Fixed:
            torque = j.angularImpulse;
            break;
        case JointType::Hinge:
        case JointType::Slider:
            force += j.rowForce[0] * axial;
            torque = j.angularImpulse + j.rowTorque[0] * axial;
            break;
        case JointType::Distance:
            force = j.rowForce[0] * axial;
            break;
        case JointType::SwingTwist:
            force += j.rowForce[0] * axial + j.rowForce[1] * j.swingImpulse;
            torque = j.rowTorque[0] * axial + j.rowTorque[1] * j.swingImpulse + j.frictionTurn;
            break;
    }
    j.force += force;
    j.torque += torque;
}

}
