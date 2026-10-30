#include "lanes.h"

namespace twofold {

namespace {

constexpr uint32_t kAll = 0xffffffffu;

bool moving(const Body& b) { return b.type == MotionType::Dynamic && b.awake; }

struct LaneMat {
    LaneVec r[3];
};

LaneVec operator*(const LaneMat& m, const LaneVec& v) { return { dot(m.r[0], v), dot(m.r[1], v), dot(m.r[2], v) }; }

struct LaneMass {
    LaneMat linear, angular;
    Lanes inverse;
};

LaneMass gatherMass(const std::vector<Body>& bodies, const uint32_t* index) {
    LaneMass m;
    for (int r = 0; r < 3; ++r) {
        m.linear.r[r] = gatherVec([&](int l) -> const Vec3& { return bodies[index[l]].invMassWorld.r[r]; });
        m.angular.r[r] = gatherVec([&](int l) -> const Vec3& { return bodies[index[l]].invInertiaWorld.r[r]; });
    }
    m.inverse = gather([&](int l) { return bodies[index[l]].invMass; });
    return m;
}

struct LaneDirection {
    LaneVec armA, armB, linearA, angularA, linearB, angularB;
    Lanes mass;
};

// makeAxis in every lane.
LaneDirection makeDirection(const LaneMass& a, const LaneMass& b, const LaneVec& ra, const LaneVec& rb, const LaneVec& d) {
    LaneDirection x;
    x.armA = cross(ra, d);
    x.armB = cross(rb, d);
    x.linearA = a.linear * d;
    x.angularA = a.angular * x.armA;
    x.linearB = b.linear * d;
    x.angularB = b.angular * x.armB;
    const Lanes k = dot(d, x.linearA) + dot(d, x.linearB) + dot(x.armA, x.angularA) + dot(x.armB, x.angularB);
    x.mass = select(k > lanes(1e-6f) * (a.inverse + b.inverse), lanes(1.0f) / k, lanes(0.0f));
    return x;
}

void storeAxis(LaneAxis& out, const LaneDirection& x) {
    storeVec(out.armA, x.armA);
    storeVec(out.armB, x.armB);
    storeVec(out.linearA, x.linearA);
    storeVec(out.angularA, x.angularA);
    storeVec(out.linearB, x.linearB);
    storeVec(out.angularB, x.angularB);
    store(out.mass, x.mass);
}

LaneDirection loadAxis(const LaneAxis& x) {
    return { loadVec(x.armA), loadVec(x.armB), loadVec(x.linearA), loadVec(x.angularA), loadVec(x.linearB), loadVec(x.angularB),
             load(x.mass) };
}

void prepareBatch(SolverData& s, const SolverContext& ctx, ContactLanes& L) {
    const Contact* c[kLanes];
    L.points = 0;
    for (int l = 0; l < kLanes; ++l) {
        const bool real = l < L.size;
        if (!real) { L.contact[l] = L.contact[0]; }
        c[l] = &s.contacts[L.contact[l]];
        const Body& a = s.bodies[c[l]->bodyA];
        const Body& b = s.bodies[c[l]->bodyB];
        L.bodyA[l] = c[l]->bodyA;
        L.bodyB[l] = c[l]->bodyB;
        L.moveA[l] = real && moving(a) ? kAll : 0u;
        L.moveB[l] = real && moving(b) ? kAll : 0u;
        L.count[l] = real ? (uint32_t)c[l]->count : 0u;
        L.friction[l] = real ? c[l]->friction : 0.0f;
        const Softness& soft = moving(a) && moving(b) ? ctx.contactSoft : ctx.staticSoft;
        L.biasRate[l] = soft.biasRate;
        L.massScale[l] = soft.massScale;
        L.impulseScale[l] = soft.impulseScale;
        L.points = (int)L.count[l] > L.points ? (int)L.count[l] : L.points;
    }
    L.even = true;
    for (int l = 0; l < L.size; ++l) { L.even = L.even && (int)L.count[l] == L.points; }
    if (L.points == 0) { return; }

    const LaneVec normal = gatherVec([&](int l) -> const Vec3& { return c[l]->normal; });
    const LaneVec tangent1 = gatherVec([&](int l) -> const Vec3& { return c[l]->tangent1; });
    const LaneVec tangent2 = gatherVec([&](int l) -> const Vec3& { return c[l]->tangent2; });
    storeVec(L.normal, normal);
    storeVec(L.tangent1, tangent1);
    storeVec(L.tangent2, tangent2);
    const LaneMass a = gatherMass(s.bodies, L.bodyA), b = gatherMass(s.bodies, L.bodyB);
    for (int i = 0; i < L.points; ++i) {
        for (int l = 0; l < kLanes; ++l) { L.used[i][l] = i < (int)L.count[l] ? kAll : 0u; }
        auto on = [&](int l) { return i < (int)L.count[l]; };
        const LaneVec ra = gatherVec([&](int l) { return on(l) ? c[l]->points[i].anchorA : Vec3{}; });
        const LaneVec rb = gatherVec([&](int l) { return on(l) ? c[l]->points[i].anchorB : Vec3{}; });
        storeVec(L.anchorA[i], ra);
        storeVec(L.anchorB[i], rb);
        store(L.adjustedSeparation[i], gather([&](int l) { return on(l) ? c[l]->points[i].adjustedSeparation : 0.0f; }));
        storeAxis(L.normalAxis[i], makeDirection(a, b, ra, rb, normal));
        storeAxis(L.tangentAxis[i][0], makeDirection(a, b, ra, rb, tangent1));
        storeAxis(L.tangentAxis[i][1], makeDirection(a, b, ra, rb, tangent2));
    }
}

struct LaneBody {
    LaneVec v, w, dp, q;
    Lanes qw;
    // Unturned this substep, so its anchors are taken as they are.
    LaneMask still;
};

LaneBody gatherBody(const std::vector<Body>& bodies, const uint32_t* index) {
    auto body = [&](int l) -> const Body& { return bodies[index[l]]; };
    LaneBody out;
    out.v = gatherVec([&](int l) -> const Vec3& { return body(l).v; });
    out.w = gatherVec([&](int l) -> const Vec3& { return body(l).w; });
    out.dp = gatherVec([&](int l) -> const Vec3& { return body(l).dp; });
    out.q = { gather([&](int l) { return body(l).dq.x; }), gather([&](int l) { return body(l).dq.y; }),
              gather([&](int l) { return body(l).dq.z; }) };
    out.qw = gather([&](int l) { return body(l).dq.w; });
    const Lanes zero = lanes(0.0f);
    out.still = (out.q.x == zero) & (out.q.y == zero) & (out.q.z == zero);
    return out;
}

void scatterBody(std::vector<Body>& bodies, const uint32_t* index, const uint32_t* move, const LaneBody& b) {
    float v[3][kLanes], w[3][kLanes];
    storeVec(v, b.v);
    storeVec(w, b.w);
    for (int l = 0; l < kLanes; ++l) {
        if (!move[l]) { continue; }
        Body& x = bodies[index[l]];
        x.v = { v[0][l], v[1][l], v[2][l] };
        x.w = { w[0][l], w[1][l], w[2][l] };
    }
}

void copyLane(LaneAxis& out, const LaneAxis& from, int l) {
    for (int k = 0; k < 3; ++k) {
        out.armA[k][l] = from.armA[k][l];
        out.armB[k][l] = from.armB[k][l];
        out.linearA[k][l] = from.linearA[k][l];
        out.angularA[k][l] = from.angularA[k][l];
        out.linearB[k][l] = from.linearB[k][l];
        out.angularB[k][l] = from.angularB[k][l];
    }
    out.mass[l] = from.mass[l];
}

Lanes speedAlong(const LaneBody& a, const LaneBody& b, const LaneVec& d, const LaneDirection& x) {
    return dot(d, b.v - a.v) + dot(x.armB, b.w) - dot(x.armA, a.w);
}

void push(LaneBody& a, LaneBody& b, const LaneDirection& x, Lanes impulse, LaneMask moveA, LaneMask moveB) {
    a.v = select(moveA, a.v - x.linearA * impulse, a.v);
    a.w = select(moveA, a.w - x.angularA * impulse, a.w);
    b.v = select(moveB, b.v + x.linearB * impulse, b.v);
    b.w = select(moveB, b.w + x.angularB * impulse, b.w);
}

LaneVec turned(const LaneBody& b, const LaneVec& r) { return select(b.still, r, rotate(b.q, b.qw, r)); }

// The point of each lane's contact, or null where the lane has no point i.
struct LanePoints {
    const SolverData& s;
    const ContactLanes& L;
    const ContactPoint* operator()(int l, int i) const {
        return i < (int)L.count[l] ? &s.contacts[L.contact[l]].points[i] : nullptr;
    }
};

void storeNormalImpulses(SolverData& s, const ContactLanes& L, const Lanes* normalImpulse, const Lanes* maxNormalImpulse) {
    float n[4][kLanes], m[4][kLanes];
    for (int i = 0; i < L.points; ++i) {
        store(n[i], normalImpulse[i]);
        store(m[i], maxNormalImpulse[i]);
    }
    for (int l = 0; l < L.size; ++l) {
        Contact& c = s.contacts[L.contact[l]];
        for (int i = 0; i < (int)L.count[l]; ++i) {
            c.points[i].normalImpulse = n[i][l];
            c.points[i].maxNormalImpulse = m[i][l];
        }
    }
}

// The lane form of warmStart.
void warmStartLanes(SolverData& s, const ContactLanes& L) {
    LaneBody a = gatherBody(s.bodies, L.bodyA);
    LaneBody b = gatherBody(s.bodies, L.bodyB);
    const LaneMask moveA = loadMask(L.moveA), moveB = loadMask(L.moveB);
    const LaneVec tangent1 = loadVec(L.tangent1), tangent2 = loadVec(L.tangent2);
    const LanePoints point{ s, L };
    for (int i = 0; i < L.points; ++i) {
        const Lanes normalImpulse = gather([&](int l) {
            const ContactPoint* p = point(l, i);
            return p ? p->normalImpulse : 0.0f;
        });
        const LaneVec held = gatherVec([&](int l) {
            const ContactPoint* p = point(l, i);
            return p ? p->tangentImpulse : Vec3{};
        });
        const LaneMask used = loadMask(L.used[i]);
        push(a, b, loadAxis(L.normalAxis[i]), normalImpulse, used & moveA, used & moveB);
        push(a, b, loadAxis(L.tangentAxis[i][0]), dot(held, tangent1), used & moveA, used & moveB);
        push(a, b, loadAxis(L.tangentAxis[i][1]), dot(held, tangent2), used & moveA, used & moveB);
    }
    scatterBody(s.bodies, L.bodyA, L.moveA, a);
    scatterBody(s.bodies, L.bodyB, L.moveB, b);
}

// The lane form of solveContacts, point for point.
void solveLanes(SolverData& s, const ContactLanes& L, const SolverContext& ctx, bool useBias) {
    LaneBody a = gatherBody(s.bodies, L.bodyA);
    LaneBody b = gatherBody(s.bodies, L.bodyB);
    const LaneMask moveA = loadMask(L.moveA), moveB = loadMask(L.moveB);
    const LaneVec normal = loadVec(L.normal), tangent1 = loadVec(L.tangent1), tangent2 = loadVec(L.tangent2);
    const Lanes zero = lanes(0.0f), one = lanes(1.0f);
    const Lanes invH = lanes(ctx.invH), pushout = lanes(-ctx.maxPushout);
    const Lanes biasRate = load(L.biasRate), massScaleSoft = load(L.massScale), impulseScaleSoft = load(L.impulseScale);
    const LanePoints point{ s, L };

    Lanes normalImpulse[4] = {}, maxNormalImpulse[4] = {};
    const LaneVec dp = b.dp - a.dp;
    for (int i = 0; i < L.points; ++i) {
        const LaneMask used = loadMask(L.used[i]);
        const Lanes held = gather([&](int l) {
            const ContactPoint* p = point(l, i);
            return p ? p->normalImpulse : 0.0f;
        });
        const Lanes peak = gather([&](int l) {
            const ContactPoint* p = point(l, i);
            return p ? p->maxNormalImpulse : 0.0f;
        });
        const LaneVec d = dp + turned(b, loadVec(L.anchorB[i])) - turned(a, loadVec(L.anchorA[i]));
        const Lanes sep = dot(d, normal) + load(L.adjustedSeparation[i]);
        const LaneMask ahead = sep > zero;
        Lanes bias = select(ahead, sep * invH, zero), massScale = one, impulseScale = zero;
        if (useBias) {
            bias = select(ahead, bias, laneMax(biasRate * sep, pushout));
            massScale = select(ahead, one, massScaleSoft);
            impulseScale = select(ahead, zero, impulseScaleSoft);
        }

        const LaneDirection x = loadAxis(L.normalAxis[i]);
        const Lanes vn = speedAlong(a, b, normal, x);
        Lanes impulse = -x.mass * massScale * (vn + bias) - impulseScale * held;
        const Lanes total = laneMax(held + impulse, zero);
        impulse = total - held;
        normalImpulse[i] = select(used, total, held);
        maxNormalImpulse[i] = select(used, laneMax(peak, impulse), peak);
        push(a, b, x, impulse, used & moveA, used & moveB);
    }

    // The relax pass takes friction points in reverse, as solveContacts does. Where the lanes' counts differ,
    // each lane finds its own point.
    float limits[4][kLanes];
    const bool even = useBias || L.even;
    if (!even) {
        for (int i = 0; i < L.points; ++i) { store(limits[i], normalImpulse[i]); }
    }
    auto order = [&](int l, int j) { return useBias ? j : (int)L.count[l] - 1 - j; };
    LaneVec tangentImpulse[4];
    const Lanes friction = load(L.friction);
    for (int j = 0; j < L.points; ++j) {
        const LaneMask used = loadMask(L.used[j]);
        LaneAxis mixed[2];
        const int at = useBias ? j : L.points - 1 - j;
        const LaneAxis* axes = L.tangentAxis[at];
        Lanes limitImpulse = normalImpulse[at];
        if (!even) {
            for (int l = 0; l < kLanes; ++l) {
                const int i = j < (int)L.count[l] ? order(l, j) : 0;
                for (int t = 0; t < 2; ++t) { copyLane(mixed[t], L.tangentAxis[i][t], l); }
            }
            axes = mixed;
            limitImpulse = gather([&](int l) { return j < (int)L.count[l] ? limits[order(l, j)][l] : 0.0f; });
        }
        const LaneDirection x1 = loadAxis(axes[0]);
        const LaneDirection x2 = loadAxis(axes[1]);
        const LaneVec held = gatherVec([&](int l) {
            const ContactPoint* p = j < (int)L.count[l] ? point(l, order(l, j)) : nullptr;
            return p ? p->tangentImpulse : Vec3{};
        });
        const Lanes old1 = dot(held, tangent1), old2 = dot(held, tangent2);
        Lanes t1 = old1 - x1.mass * speedAlong(a, b, tangent1, x1);
        Lanes t2 = old2 - x2.mass * speedAlong(a, b, tangent2, x2);
        const Lanes limit = friction * limitImpulse;
        const Lanes mag = laneSqrt(t1 * t1 + t2 * t2);
        const LaneMask over = mag > limit;
        const Lanes scale = select(mag > zero, limit / mag, zero);
        t1 = select(over, t1 * scale, t1);
        t2 = select(over, t2 * scale, t2);
        tangentImpulse[j] = select(used, tangent1 * t1 + tangent2 * t2, held);
        push(a, b, x1, t1 - old1, used & moveA, used & moveB);
        push(a, b, x2, t2 - old2, used & moveA, used & moveB);
    }

    storeNormalImpulses(s, L, normalImpulse, maxNormalImpulse);
    float t[4][3][kLanes];
    for (int j = 0; j < L.points; ++j) { storeVec(t[j], tangentImpulse[j]); }
    for (int l = 0; l < L.size; ++l) {
        Contact& c = s.contacts[L.contact[l]];
        for (int j = 0; j < (int)L.count[l]; ++j) { c.points[order(l, j)].tangentImpulse = { t[j][0][l], t[j][1][l], t[j][2][l] }; }
    }
    scatterBody(s.bodies, L.bodyA, L.moveA, a);
    scatterBody(s.bodies, L.bodyB, L.moveB, b);
}

// The lane form of applyRestitution.
void restituteLanes(SolverData& s, const ContactLanes& L, const SolverContext& ctx) {
    bool any = false;
    for (int l = 0; l < L.size; ++l) { any = any || s.contacts[L.contact[l]].restitution != 0.0f; }
    if (!any) { return; }
    LaneBody a = gatherBody(s.bodies, L.bodyA);
    LaneBody b = gatherBody(s.bodies, L.bodyB);
    const LaneMask moveA = loadMask(L.moveA), moveB = loadMask(L.moveB);
    const LaneVec normal = loadVec(L.normal);
    const Lanes zero = lanes(0.0f), threshold = lanes(-ctx.restitutionThreshold);
    const Lanes e = gather([&](int l) { return l < L.size ? s.contacts[L.contact[l]].restitution : 0.0f; });
    const LanePoints point{ s, L };
    Lanes normalImpulse[4] = {}, maxNormalImpulse[4] = {};
    for (int i = 0; i < L.points; ++i) {
        auto field = [&](float ContactPoint::*f) {
            return gather([&](int l) {
                const ContactPoint* p = point(l, i);
                return p ? p->*f : 0.0f;
            });
        };
        const Lanes rv = field(&ContactPoint::relativeVelocity);
        const Lanes held = field(&ContactPoint::normalImpulse), peak = field(&ContactPoint::maxNormalImpulse);
        const LaneMask on = loadMask(L.used[i]) & !(e == zero) & !(rv > threshold) & !(peak == zero);
        const LaneDirection x = loadAxis(L.normalAxis[i]);
        const Lanes vn = speedAlong(a, b, normal, x);
        Lanes impulse = -x.mass * (vn + e * rv);
        const Lanes total = laneMax(held + impulse, zero);
        impulse = total - held;
        normalImpulse[i] = select(on, total, held);
        maxNormalImpulse[i] = select(on, laneMax(peak, impulse), peak);
        push(a, b, x, impulse, on & moveA, on & moveB);
    }
    storeNormalImpulses(s, L, normalImpulse, maxNormalImpulse);
    scatterBody(s.bodies, L.bodyA, L.moveA, a);
    scatterBody(s.bodies, L.bodyB, L.moveB, b);
}

}

void colorContacts(SolverData& s, const std::vector<uint32_t>& list, ContactColors& colors) {
    colors.bodies.assign(s.bodies.size(), 0u);
    colors.of.resize(s.contacts.size());
    for (uint32_t ci : list) {
        Contact& c = s.contacts[ci];
        colors.of[ci] = (uint8_t)kColors;
        c.laned = false;
        if (c.planar || c.rollingResistance > 0.0f) { continue; }
        const bool ma = moving(s.bodies[c.bodyA]), mb = moving(s.bodies[c.bodyB]);
        const uint32_t taken = (ma ? colors.bodies[c.bodyA] : 0u) | (mb ? colors.bodies[c.bodyB] : 0u);
        if (taken == kAll) { continue; }
        int k = 0;
        while (taken & (1u << k)) { ++k; }
        if (ma) { colors.bodies[c.bodyA] |= 1u << k; }
        if (mb) { colors.bodies[c.bodyB] |= 1u << k; }
        colors.of[ci] = (uint8_t)k;
        c.laned = true;
    }
}

void batchContacts(SolverData& s, const std::vector<uint32_t>& list, const ContactColors& colors, const SolverContext& ctx,
                   ContactBatches& out) {
    out.lanes.clear();
    std::vector<uint32_t>& byColor = out.byColor;
    uint32_t sizes[kColors + 1] = {};
    for (uint32_t ci : list) { ++sizes[colors.of[ci]]; }
    uint32_t start[kColors + 1];
    uint32_t sum = 0;
    for (int k = 0; k <= kColors; ++k) {
        start[k] = sum;
        sum += sizes[k];
    }
    byColor.resize(list.size());
    for (int points = 4; points >= 0; --points) {
        for (uint32_t ci : list) {
            const uint8_t k = colors.of[ci];
            if (k == kColors ? points == 4 : s.contacts[ci].count == points) { byColor[start[k]++] = ci; }
        }
    }

    size_t n = 0;
    for (int k = 0; k < kColors; ++k) {
        for (uint32_t done = 0; done < sizes[k]; done += kLanes) {
            ContactLanes& L = out.lanes.emplace_back();
            const uint32_t take = sizes[k] - done < (uint32_t)kLanes ? sizes[k] - done : (uint32_t)kLanes;
            L.size = (int)take;
            for (uint32_t l = 0; l < take; ++l) { L.contact[l] = byColor[n++]; }
        }
    }
    out.rest.assign(byColor.begin() + (std::ptrdiff_t)n, byColor.end());
    prepareLanes(s, ctx, out);
}

void prepareLanes(SolverData& s, const SolverContext& ctx, ContactBatches& batches) {
    for (ContactLanes& L : batches.lanes) { prepareBatch(s, ctx, L); }
}

void warmStartBatches(SolverData& s, const ContactBatches& batches) {
    for (const ContactLanes& L : batches.lanes) {
        if (L.points > 0) { warmStartLanes(s, L); }
    }
}

void solveBatches(SolverData& s, const ContactBatches& batches, const SolverContext& ctx, bool useBias) {
    for (const ContactLanes& L : batches.lanes) {
        if (L.points > 0) { solveLanes(s, L, ctx, useBias); }
    }
    solveContacts(s, batches.rest, ctx, useBias);
}

void restituteBatches(SolverData& s, const ContactBatches& batches, const SolverContext& ctx) {
    for (const ContactLanes& L : batches.lanes) {
        if (L.points > 0) { restituteLanes(s, L, ctx); }
    }
}

}
