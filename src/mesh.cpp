#include <algorithm>
#include <cmath>

#include "internal.h"

namespace twofold {

namespace {

constexpr int kLeafSize = 4;

// A neighbor folding away by less than this sine leaves the shared edge without a normal of its own.
constexpr float kFlatEdge = 0.087f;

// A triangle whose widest corner has a smaller sine is a sliver, its normal left to rounding.
constexpr float kSliver = 1e-3f;

bool isFinite(const Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

bool lessXyz(const Vec3& a, const Vec3& b) {
    if (a.x != b.x) { return a.x < b.x; }
    if (a.y != b.y) { return a.y < b.y; }
    return a.z < b.z;
}

struct Builder {
    TriangleMesh::Data& d;
    std::vector<int> order;
    std::vector<Vec3> centers;

    void fit(TriangleMesh::Data::Node& n, int first, int count) {
        Vec3 c[3];
        d.corners(order[first], c);
        n.lo = n.hi = c[0];
        for (int i = first; i < first + count; ++i) {
            d.corners(order[i], c);
            for (const Vec3& v : c) {
                n.lo = minVec(n.lo, v);
                n.hi = maxVec(n.hi, v);
            }
        }
    }

    // Splits at the median along the widest spread of triangle centers.
    void build(int node, int first, int count) {
        fit(d.nodes[node], first, count);
        if (count <= kLeafSize) {
            d.nodes[node].first = first;
            d.nodes[node].count = count;
            return;
        }
        Vec3 lo = centers[order[first]], hi = lo;
        for (int i = first; i < first + count; ++i) {
            lo = minVec(lo, centers[order[i]]);
            hi = maxVec(hi, centers[order[i]]);
        }
        const Vec3 spread = hi - lo;
        const int axis = spread.x >= spread.y && spread.x >= spread.z ? 0 : spread.y >= spread.z ? 1 : 2;
        std::sort(order.begin() + first, order.begin() + first + count, [&](int a, int b) {
            const float ca = centers[a][axis], cb = centers[b][axis];
            return ca != cb ? ca < cb : a < b;
        });
        const int child = (int)d.nodes.size();
        d.nodes.emplace_back();
        d.nodes.emplace_back();
        d.nodes[node].first = child;
        d.nodes[node].count = 0;
        const int half = count / 2;
        build(child, first, half);
        build(child + 1, first + half, count - half);
    }
};

uint8_t edgeFlags(const Vec3& a, const Vec3& b, const Vec3& normal, const Vec3& far) {
    const Vec3 e = normalize(b - a);
    const Vec3 r = far - a;
    const float reach = length(r - e * dot(r, e));
    const float rise = dot(normal, r);
    uint8_t flags = 0;
    if (rise < -kFlatEdge * reach) { flags |= 1; }
    if (rise > kFlatEdge * reach) { flags |= 8; }
    return flags;
}

bool rayTriangle(const Vec3 c[3], const Vec3& o, const Vec3& d, float maxFraction, float& fraction) {
    const Vec3 e1 = c[1] - c[0], e2 = c[2] - c[0];
    const Vec3 p = cross(d, e2);
    const float det = dot(e1, p);
    if (det == 0.0f) { return false; }
    const float inv = 1.0f / det;
    const Vec3 s = o - c[0];
    const float u = dot(s, p) * inv;
    if (u < 0.0f || u > 1.0f) { return false; }
    const Vec3 q = cross(s, e1);
    const float v = dot(d, q) * inv;
    if (v < 0.0f || u + v > 1.0f) { return false; }
    const float t = dot(e2, q) * inv;
    if (t < 0.0f || t > maxFraction) { return false; }
    fraction = t;
    return true;
}

bool rayBox(const Vec3& lo, const Vec3& hi, const Vec3& o, const Vec3& d, float maxFraction) {
    float lower = 0.0f, upper = maxFraction;
    for (int i = 0; i < 3; ++i) {
        if (d[i] == 0.0f) {
            if (o[i] < lo[i] || o[i] > hi[i]) { return false; }
            continue;
        }
        const float inv = 1.0f / d[i];
        float t0 = (lo[i] - o[i]) * inv, t1 = (hi[i] - o[i]) * inv;
        if (t0 > t1) { std::swap(t0, t1); }
        lower = maxf(lower, t0);
        upper = minf(upper, t1);
        if (lower > upper) { return false; }
    }
    return true;
}

}

TriangleMesh::TriangleMesh(const std::vector<Vec3>& vertices, const std::vector<uint32_t>& indices)
    : data_(std::make_unique<Data>()) {
    Data& d = *data_;

    std::vector<uint32_t> sorted;
    sorted.reserve(vertices.size());
    for (uint32_t i = 0; i < (uint32_t)vertices.size(); ++i) {
        if (isFinite(vertices[i])) { sorted.push_back(i); }
    }
    std::sort(sorted.begin(), sorted.end(), [&](uint32_t a, uint32_t b) {
        if (lessXyz(vertices[a], vertices[b])) { return true; }
        if (lessXyz(vertices[b], vertices[a])) { return false; }
        return a < b;
    });
    std::vector<uint32_t> first(vertices.size(), UINT32_MAX);
    for (size_t i = 0; i < sorted.size(); ++i) {
        const bool same = i > 0 && !lessXyz(vertices[sorted[i - 1]], vertices[sorted[i]]);
        first[sorted[i]] = same ? first[sorted[i - 1]] : sorted[i];
    }

    std::vector<uint32_t> remap(vertices.size(), UINT32_MAX);
    std::vector<uint32_t> tris;
    for (size_t t = 0; t + 2 < indices.size(); t += 3) {
        uint32_t v[3];
        bool valid = true;
        for (int k = 0; k < 3; ++k) {
            if (indices[t + k] >= vertices.size() || first[indices[t + k]] == UINT32_MAX) {
                valid = false;
                break;
            }
            v[k] = first[indices[t + k]];
        }
        if (!valid || v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) { continue; }
        int wide = 0;
        float longest = 0.0f;
        for (int k = 0; k < 3; ++k) {
            const float opposite = lengthSquared(vertices[v[(k + 2) % 3]] - vertices[v[(k + 1) % 3]]);
            if (opposite > longest) {
                longest = opposite;
                wide = k;
            }
        }
        const Vec3 e1 = vertices[v[(wide + 1) % 3]] - vertices[v[wide]];
        const Vec3 e2 = vertices[v[(wide + 2) % 3]] - vertices[v[wide]];
        if (!(lengthSquared(cross(e1, e2)) > kSliver * kSliver * lengthSquared(e1) * lengthSquared(e2))) { continue; }
        for (int k = 0; k < 3; ++k) {
            if (remap[v[k]] == UINT32_MAX) {
                remap[v[k]] = (uint32_t)d.vertices.size();
                d.vertices.push_back(vertices[v[k]]);
            }
            tris.push_back(remap[v[k]]);
        }
    }
    const int count = (int)tris.size() / 3;
    if (count == 0) { return; }

    std::vector<Vec3> normals(count);
    std::vector<uint8_t> edges(count, 0);
    for (int t = 0; t < count; ++t) {
        const Vec3 a = d.vertices[tris[3 * t]], b = d.vertices[tris[3 * t + 1]], c = d.vertices[tris[3 * t + 2]];
        normals[t] = normalize(cross(b - a, c - a));
    }

    // Edges keyed by their two vertices. An edge with exactly one neighbor gets its own normal only on the
    // side where the neighbor folds away. Open and shared by more than two triangles, it gets one on both.
    struct EdgeRef {
        uint64_t key;
        int tri;
        int corner;
    };
    std::vector<EdgeRef> refs;
    refs.reserve(tris.size());
    for (int t = 0; t < count; ++t) {
        for (int k = 0; k < 3; ++k) {
            const uint32_t a = tris[3 * t + k], b = tris[3 * t + (k + 1) % 3];
            refs.push_back({ ((uint64_t)std::min(a, b) << 32) | std::max(a, b), t, k });
        }
    }
    std::sort(refs.begin(), refs.end(), [](const EdgeRef& x, const EdgeRef& y) {
        if (x.key != y.key) { return x.key < y.key; }
        return x.tri != y.tri ? x.tri < y.tri : x.corner < y.corner;
    });
    for (size_t i = 0; i < refs.size();) {
        size_t j = i;
        while (j < refs.size() && refs[j].key == refs[i].key) { ++j; }
        for (size_t k = i; k < j; ++k) {
            const EdgeRef& r = refs[k];
            if (j - i != 2) {
                edges[r.tri] |= (uint8_t)(9u << r.corner);
                continue;
            }
            const EdgeRef& o = refs[k == i ? i + 1 : i];
            const Vec3 a = d.vertices[tris[3 * r.tri + r.corner]];
            const Vec3 b = d.vertices[tris[3 * r.tri + (r.corner + 1) % 3]];
            const Vec3 far = d.vertices[tris[3 * o.tri + (o.corner + 2) % 3]];
            edges[r.tri] |= (uint8_t)(edgeFlags(a, b, normals[r.tri], far) << r.corner);
        }
        i = j;
    }

    d.indices = tris;
    Builder builder{ d, std::vector<int>(count), std::vector<Vec3>(count) };
    for (int t = 0; t < count; ++t) {
        builder.order[t] = t;
        Vec3 c[3];
        d.corners(t, c);
        builder.centers[t] = (c[0] + c[1] + c[2]) * (1.0f / 3.0f);
    }
    d.nodes.reserve(2 * (size_t)count / kLeafSize + 1);
    d.nodes.emplace_back();
    builder.build(0, 0, count);

    d.indices.resize(tris.size());
    d.normals.resize(count);
    d.edges.resize(count);
    for (int i = 0; i < count; ++i) {
        const int t = builder.order[i];
        for (int k = 0; k < 3; ++k) { d.indices[3 * i + k] = tris[3 * t + k]; }
        d.normals[i] = normals[t];
        d.edges[i] = edges[t];
    }
}

TriangleMesh::~TriangleMesh() = default;

const std::vector<Vec3>& TriangleMesh::vertices() const { return data_->vertices; }

const std::vector<uint32_t>& TriangleMesh::indices() const { return data_->indices; }

int TriangleMesh::triangleCount() const { return (int)data_->normals.size(); }

const TriangleMesh::Data* meshOf(const Shape& s) {
    const Mesh* m = std::get_if<Mesh>(&s.geometry);
    return m && m->mesh && !m->mesh->data().nodes.empty() ? &m->mesh->data() : nullptr;
}

void boundsInMesh(const Transform& xm, const Vec3& lo, const Vec3& hi, Vec3& meshLo, Vec3& meshHi) {
    const Vec3 c = xm.applyInverse((lo + hi) * 0.5f);
    const Vec3 e = (hi - lo) * 0.5f;
    const Mat3 r = xm.q.conjugate().toMat3();
    Vec3 reach;
    for (int i = 0; i < 3; ++i) {
        reach[i] = std::fabs(r.r[i].x) * e.x + std::fabs(r.r[i].y) * e.y + std::fabs(r.r[i].z) * e.z;
    }
    meshLo = c - reach;
    meshHi = c + reach;
}

bool rayMesh(const TriangleMesh::Data& d, const Vec3& o, const Vec3& dir, float maxFraction, float& fraction, Vec3& normal) {
    int stack[64];
    int top = 0;
    stack[top++] = 0;
    float best = maxFraction;
    int hit = -1;
    while (top > 0) {
        const TriangleMesh::Data::Node& n = d.nodes[stack[--top]];
        if (!rayBox(n.lo, n.hi, o, dir, best)) { continue; }
        if (n.count == 0) {
            stack[top++] = n.first + 1;
            stack[top++] = n.first;
            continue;
        }
        for (int t = n.first; t < n.first + n.count; ++t) {
            Vec3 c[3];
            d.corners(t, c);
            float f;
            if (rayTriangle(c, o, dir, best, f) && (hit < 0 || f < best)) {
                best = f;
                hit = t;
            }
        }
    }
    if (hit < 0) { return false; }
    fraction = best;
    normal = dot(d.normals[hit], dir) > 0.0f ? -d.normals[hit] : d.normals[hit];
    return true;
}

}
