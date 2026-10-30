#pragma once

#include <cmath>
#include <cstdint>

#include "internal.h"

#if defined(TWOFOLD_SSE2)
#include <emmintrin.h>
#elif defined(TWOFOLD_NEON)
#include <arm_neon.h>
#endif

namespace twofold {

// Each lane takes exactly the operations a float would, rounded the same way, so every lane width gives the
// same bits. No fused multiply-add and no approximate reciprocals.
#if defined(TWOFOLD_SSE2)

struct Lanes {
    __m128 v;
};
struct LaneMask {
    __m128 m;
};

inline Lanes lanes(float s) { return { _mm_set1_ps(s) }; }
inline Lanes load(const float* p) { return { _mm_loadu_ps(p) }; }
inline void store(float* p, Lanes a) { _mm_storeu_ps(p, a.v); }
inline LaneMask loadMask(const uint32_t* p) { return { _mm_castsi128_ps(_mm_loadu_si128((const __m128i*)p)) }; }
inline Lanes operator+(Lanes a, Lanes b) { return { _mm_add_ps(a.v, b.v) }; }
inline Lanes operator-(Lanes a, Lanes b) { return { _mm_sub_ps(a.v, b.v) }; }
inline Lanes operator*(Lanes a, Lanes b) { return { _mm_mul_ps(a.v, b.v) }; }
inline Lanes operator/(Lanes a, Lanes b) { return { _mm_div_ps(a.v, b.v) }; }
inline Lanes operator-(Lanes a) { return { _mm_xor_ps(a.v, _mm_set1_ps(-0.0f)) }; }
inline Lanes laneSqrt(Lanes a) { return { _mm_sqrt_ps(a.v) }; }
inline LaneMask operator>(Lanes a, Lanes b) { return { _mm_cmpgt_ps(a.v, b.v) }; }
inline LaneMask operator==(Lanes a, Lanes b) { return { _mm_cmpeq_ps(a.v, b.v) }; }
inline LaneMask isNan(Lanes a) { return { _mm_cmpunord_ps(a.v, a.v) }; }
inline LaneMask operator&(LaneMask a, LaneMask b) { return { _mm_and_ps(a.m, b.m) }; }
inline LaneMask operator!(LaneMask a) { return { _mm_xor_ps(a.m, _mm_castsi128_ps(_mm_set1_epi32(-1))) }; }
inline Lanes select(LaneMask m, Lanes a, Lanes b) { return { _mm_or_ps(_mm_and_ps(m.m, a.v), _mm_andnot_ps(m.m, b.v)) }; }
template <class F> Lanes gather(F f) { return { _mm_setr_ps(f(0), f(1), f(2), f(3)) }; }
// Bit l set where lane l is.
inline int bits(LaneMask m) { return _mm_movemask_ps(m.m); }

#elif defined(TWOFOLD_NEON)

struct Lanes {
    float32x4_t v;
};
struct LaneMask {
    uint32x4_t m;
};

inline Lanes lanes(float s) { return { vdupq_n_f32(s) }; }
inline Lanes load(const float* p) { return { vld1q_f32(p) }; }
inline void store(float* p, Lanes a) { vst1q_f32(p, a.v); }
inline LaneMask loadMask(const uint32_t* p) { return { vld1q_u32(p) }; }
inline Lanes operator+(Lanes a, Lanes b) { return { vaddq_f32(a.v, b.v) }; }
inline Lanes operator-(Lanes a, Lanes b) { return { vsubq_f32(a.v, b.v) }; }
inline Lanes operator*(Lanes a, Lanes b) { return { vmulq_f32(a.v, b.v) }; }
inline Lanes operator/(Lanes a, Lanes b) { return { vdivq_f32(a.v, b.v) }; }
inline Lanes operator-(Lanes a) { return { vnegq_f32(a.v) }; }
inline Lanes laneSqrt(Lanes a) { return { vsqrtq_f32(a.v) }; }
inline LaneMask operator>(Lanes a, Lanes b) { return { vcgtq_f32(a.v, b.v) }; }
inline LaneMask operator==(Lanes a, Lanes b) { return { vceqq_f32(a.v, b.v) }; }
inline LaneMask isNan(Lanes a) { return { vmvnq_u32(vceqq_f32(a.v, a.v)) }; }
inline LaneMask operator&(LaneMask a, LaneMask b) { return { vandq_u32(a.m, b.m) }; }
inline LaneMask operator!(LaneMask a) { return { vmvnq_u32(a.m) }; }
inline Lanes select(LaneMask m, Lanes a, Lanes b) { return { vbslq_f32(m.m, a.v, b.v) }; }
template <class F> Lanes gather(F f) {
    float32x4_t r = vdupq_n_f32(f(0));
    r = vsetq_lane_f32(f(1), r, 1);
    r = vsetq_lane_f32(f(2), r, 2);
    return { vsetq_lane_f32(f(3), r, 3) };
}
inline int bits(LaneMask m) {
    const uint32x4_t one = vandq_u32(m.m, vdupq_n_u32(1));
    return (int)(vgetq_lane_u32(one, 0) | vgetq_lane_u32(one, 1) << 1 | vgetq_lane_u32(one, 2) << 2 |
                 vgetq_lane_u32(one, 3) << 3);
}

#else

struct Lanes {
    float v;
};
struct LaneMask {
    bool m;
};

inline Lanes lanes(float s) { return { s }; }
inline Lanes load(const float* p) { return { *p }; }
inline void store(float* p, Lanes a) { *p = a.v; }
inline LaneMask loadMask(const uint32_t* p) { return { *p != 0 }; }
inline Lanes operator+(Lanes a, Lanes b) { return { a.v + b.v }; }
inline Lanes operator-(Lanes a, Lanes b) { return { a.v - b.v }; }
inline Lanes operator*(Lanes a, Lanes b) { return { a.v * b.v }; }
inline Lanes operator/(Lanes a, Lanes b) { return { a.v / b.v }; }
inline Lanes operator-(Lanes a) { return { -a.v }; }
inline Lanes laneSqrt(Lanes a) { return { std::sqrt(a.v) }; }
inline LaneMask operator>(Lanes a, Lanes b) { return { a.v > b.v }; }
inline LaneMask operator==(Lanes a, Lanes b) { return { a.v == b.v }; }
inline LaneMask isNan(Lanes a) { return { a.v != a.v }; }
inline LaneMask operator&(LaneMask a, LaneMask b) { return { a.m && b.m }; }
inline LaneMask operator!(LaneMask a) { return { !a.m }; }
inline Lanes select(LaneMask m, Lanes a, Lanes b) { return m.m ? a : b; }
template <class F> Lanes gather(F f) { return { f(0) }; }
inline int bits(LaneMask m) { return m.m ? 1 : 0; }

#endif

// maxf in every lane.
inline Lanes laneMax(Lanes a, Lanes b) {
    Lanes r = select(a > b, a, b);
    r = select(a == b, select(a == lanes(0.0f), a + b, a), r);
    r = select(isNan(b), a, r);
    return select(isNan(a), b, r);
}

// minf in every lane.
inline Lanes laneMin(Lanes a, Lanes b) {
    Lanes r = select(b > a, a, b);
    r = select(a == b, select(a == lanes(0.0f), -(-a + -b), a), r);
    r = select(isNan(b), a, r);
    return select(isNan(a), b, r);
}

struct LaneVec {
    Lanes x, y, z;
};

inline LaneVec lanes(const Vec3& v) { return { lanes(v.x), lanes(v.y), lanes(v.z) }; }
inline LaneVec loadVec(const float (*p)[kLanes]) { return { load(p[0]), load(p[1]), load(p[2]) }; }
// Lane l takes f(l), called once per lane.
template <class F> LaneVec gatherVec(F f) {
    Vec3 v[kLanes];
    for (int l = 0; l < kLanes; ++l) { v[l] = f(l); }
    return { gather([&](int l) { return v[l].x; }), gather([&](int l) { return v[l].y; }), gather([&](int l) { return v[l].z; }) };
}
inline void storeVec(float (*p)[kLanes], const LaneVec& v) {
    store(p[0], v.x);
    store(p[1], v.y);
    store(p[2], v.z);
}
inline LaneVec operator+(const LaneVec& a, const LaneVec& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline LaneVec operator-(const LaneVec& a, const LaneVec& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline LaneVec operator*(const LaneVec& v, Lanes s) { return { v.x * s, v.y * s, v.z * s }; }
inline Lanes dot(const LaneVec& a, const LaneVec& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline LaneVec cross(const LaneVec& a, const LaneVec& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline LaneVec select(LaneMask m, const LaneVec& a, const LaneVec& b) {
    return { select(m, a.x, b.x), select(m, a.y, b.y), select(m, a.z, b.z) };
}

// Quat::rotate in every lane.
inline LaneVec rotate(const LaneVec& q, Lanes w, const LaneVec& v) {
    const LaneVec t = cross(q, v) * lanes(2.0f);
    return v + t * w + cross(q, t);
}

}
