#pragma once

#include <cmath>

namespace twofold {

struct Vec2 {
    float x = 0, y = 0;

    Vec2() = default;
    constexpr Vec2(float x_, float y_) : x(x_), y(y_) {}

    Vec2 operator+(Vec2 b) const { return { x + b.x, y + b.y }; }
    Vec2 operator-(Vec2 b) const { return { x - b.x, y - b.y }; }
    Vec2 operator-() const { return { -x, -y }; }
    Vec2 operator*(float s) const { return { x * s, y * s }; }
    Vec2& operator+=(Vec2 b) { x += b.x; y += b.y; return *this; }
    Vec2& operator-=(Vec2 b) { x -= b.x; y -= b.y; return *this; }
};

inline Vec2 operator*(float s, Vec2 v) { return v * s; }
inline float dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline float cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
inline float length(Vec2 v) { return std::sqrt(dot(v, v)); }
inline Vec2 normalize(Vec2 v) {
    const float l = length(v);
    return l > 0 ? v * (1.0f / l) : Vec2{};
}

struct Vec3 {
    float x = 0, y = 0, z = 0;

    Vec3() = default;
    constexpr Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    float operator[](int i) const { return (&x)[i]; }
    float& operator[](int i) { return (&x)[i]; }

    Vec3 operator+(const Vec3& b) const { return { x + b.x, y + b.y, z + b.z }; }
    Vec3 operator-(const Vec3& b) const { return { x - b.x, y - b.y, z - b.z }; }
    Vec3 operator-() const { return { -x, -y, -z }; }
    Vec3 operator*(float s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(float s) const { return { x / s, y / s, z / s }; }
    Vec3& operator+=(const Vec3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vec3& operator-=(const Vec3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
    Vec3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
    bool operator==(const Vec3& b) const { return x == b.x && y == b.y && z == b.z; }
    bool operator!=(const Vec3& b) const { return !(*this == b); }
};

inline Vec3 operator*(float s, const Vec3& v) { return v * s; }
inline float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}
inline float lengthSquared(const Vec3& v) { return dot(v, v); }
inline float length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) {
    const float l = length(v);
    return l > 0 ? v * (1.0f / l) : Vec3{};
}
inline Vec3 minVec(const Vec3& a, const Vec3& b) {
    return { a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y, a.z < b.z ? a.z : b.z };
}
inline Vec3 maxVec(const Vec3& a, const Vec3& b) {
    return { a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y, a.z > b.z ? a.z : b.z };
}

inline Vec3 perpendicular(const Vec3& n) {
    if (std::fabs(n.x) > 0.57735f) { return normalize(Vec3{ n.y, -n.x, 0.0f }); }
    return normalize(Vec3{ 0.0f, n.z, -n.y });
}

// Row-major 3x3 matrix.
struct Mat3 {
    Vec3 r[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

    static Mat3 zero() {
        Mat3 m;
        m.r[0] = m.r[1] = m.r[2] = Vec3{};
        return m;
    }
    static Mat3 diagonal(const Vec3& d) {
        Mat3 m = zero();
        m.r[0].x = d.x; m.r[1].y = d.y; m.r[2].z = d.z;
        return m;
    }
    static Mat3 outer(const Vec3& a, const Vec3& b) {
        Mat3 m;
        m.r[0] = b * a.x; m.r[1] = b * a.y; m.r[2] = b * a.z;
        return m;
    }

    Vec3 column(int i) const { return { r[0][i], r[1][i], r[2][i] }; }

    Vec3 operator*(const Vec3& v) const { return { dot(r[0], v), dot(r[1], v), dot(r[2], v) }; }
    Mat3 operator*(const Mat3& b) const {
        Mat3 m;
        for (int i = 0; i < 3; ++i) { m.r[i] = { dot(r[i], b.column(0)), dot(r[i], b.column(1)), dot(r[i], b.column(2)) }; }
        return m;
    }
    Mat3 operator+(const Mat3& b) const {
        Mat3 m;
        for (int i = 0; i < 3; ++i) { m.r[i] = r[i] + b.r[i]; }
        return m;
    }
    Mat3 operator-(const Mat3& b) const {
        Mat3 m;
        for (int i = 0; i < 3; ++i) { m.r[i] = r[i] - b.r[i]; }
        return m;
    }
    Mat3 operator*(float s) const {
        Mat3 m;
        for (int i = 0; i < 3; ++i) { m.r[i] = r[i] * s; }
        return m;
    }
};

inline Mat3 transpose(const Mat3& m) {
    Mat3 t;
    for (int i = 0; i < 3; ++i) { t.r[i] = m.column(i); }
    return t;
}

inline Mat3 inverse(const Mat3& m) {
    const Vec3 c0 = cross(m.r[1], m.r[2]);
    const Vec3 c1 = cross(m.r[2], m.r[0]);
    const Vec3 c2 = cross(m.r[0], m.r[1]);
    float det = dot(m.r[0], c0);
    if (det == 0.0f) { return Mat3::zero(); }
    det = 1.0f / det;
    Mat3 inv;
    inv.r[0] = Vec3{ c0.x, c1.x, c2.x } * det;
    inv.r[1] = Vec3{ c0.y, c1.y, c2.y } * det;
    inv.r[2] = Vec3{ c0.z, c1.z, c2.z } * det;
    return inv;
}

struct Quat {
    float x = 0, y = 0, z = 0, w = 1;

    Quat() = default;
    constexpr Quat(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}

    static Quat axisAngle(const Vec3& axis, float angle) {
        const Vec3 a = normalize(axis) * std::sin(0.5f * angle);
        return { a.x, a.y, a.z, std::cos(0.5f * angle) };
    }

    // Rotation that turns unit vector a onto unit vector b.
    static Quat fromTo(const Vec3& a, const Vec3& b) {
        const float d = dot(a, b);
        if (d < -0.999999f) { return axisAngle(perpendicular(a), 3.14159265358979f); }
        const Vec3 c = cross(a, b);
        const float s = std::sqrt((1.0f + d) * 2.0f);
        return { c.x / s, c.y / s, c.z / s, 0.5f * s };
    }

    Vec3 vec() const { return { x, y, z }; }

    Quat operator*(const Quat& b) const {
        return { w * b.x + x * b.w + y * b.z - z * b.y,
                 w * b.y - x * b.z + y * b.w + z * b.x,
                 w * b.z + x * b.y - y * b.x + z * b.w,
                 w * b.w - x * b.x - y * b.y - z * b.z };
    }
    Quat operator-() const { return { -x, -y, -z, -w }; }
    Quat conjugate() const { return { -x, -y, -z, w }; }

    Vec3 rotate(const Vec3& v) const {
        const Vec3 q = vec();
        const Vec3 t = cross(q, v) * 2.0f;
        return v + t * w + cross(q, t);
    }
    Vec3 unrotate(const Vec3& v) const { return conjugate().rotate(v); }

    Mat3 toMat3() const {
        const float xx = x * x, yy = y * y, zz = z * z;
        const float xy = x * y, xz = x * z, yz = y * z;
        const float wx = w * x, wy = w * y, wz = w * z;
        Mat3 m;
        m.r[0] = { 1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy) };
        m.r[1] = { 2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx) };
        m.r[2] = { 2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy) };
        return m;
    }
};

inline float dot(const Quat& a, const Quat& b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
inline Quat normalize(const Quat& q) {
    const float l = std::sqrt(dot(q, q));
    return l > 0 ? Quat{ q.x / l, q.y / l, q.z / l, q.w / l } : Quat{};
}

// Rotation matrix to quaternion. The matrix must be orthonormal.
inline Quat fromMat3(const Mat3& m) {
    const float tr = m.r[0].x + m.r[1].y + m.r[2].z;
    Quat q;
    if (tr > 0.0f) {
        const float s = std::sqrt(tr + 1.0f) * 2.0f;
        q = { (m.r[2].y - m.r[1].z) / s, (m.r[0].z - m.r[2].x) / s, (m.r[1].x - m.r[0].y) / s, 0.25f * s };
    } else if (m.r[0].x > m.r[1].y && m.r[0].x > m.r[2].z) {
        const float s = std::sqrt(1.0f + m.r[0].x - m.r[1].y - m.r[2].z) * 2.0f;
        q = { 0.25f * s, (m.r[0].y + m.r[1].x) / s, (m.r[0].z + m.r[2].x) / s, (m.r[2].y - m.r[1].z) / s };
    } else if (m.r[1].y > m.r[2].z) {
        const float s = std::sqrt(1.0f + m.r[1].y - m.r[0].x - m.r[2].z) * 2.0f;
        q = { (m.r[0].y + m.r[1].x) / s, 0.25f * s, (m.r[1].z + m.r[2].y) / s, (m.r[0].z - m.r[2].x) / s };
    } else {
        const float s = std::sqrt(1.0f + m.r[2].z - m.r[0].x - m.r[1].y) * 2.0f;
        q = { (m.r[0].z + m.r[2].x) / s, (m.r[1].z + m.r[2].y) / s, 0.25f * s, (m.r[1].x - m.r[0].y) / s };
    }
    return normalize(q);
}

// Angular velocity that turns a into b over dt, the short way round.
inline Vec3 angularVelocity(const Quat& a, const Quat& b, float dt) {
    Quat d = normalize(b * a.conjugate());
    if (d.w < 0) { d = -d; }
    const float s = length(d.vec());
    if (s < 1e-7f || dt <= 0) { return {}; }
    const float angle = 2.0f * std::atan2(s, d.w);
    return d.vec() * (angle / (s * dt));
}

// Advances q by angular velocity w over h.
inline Quat integrate(const Quat& q, const Vec3& w, float h) {
    const Quat dq = Quat{ w.x, w.y, w.z, 0.0f } * q;
    return normalize(Quat{ q.x + 0.5f * h * dq.x, q.y + 0.5f * h * dq.y, q.z + 0.5f * h * dq.z, q.w + 0.5f * h * dq.w });
}

struct Transform {
    Vec3 p;
    Quat q;

    Vec3 apply(const Vec3& v) const { return p + q.rotate(v); }
    Vec3 applyInverse(const Vec3& v) const { return q.unrotate(v - p); }
    Transform operator*(const Transform& b) const { return { apply(b.p), q * b.q }; }
};

}
