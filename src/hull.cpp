#include <algorithm>
#include <cfloat>
#include <cmath>
#include <unordered_map>

#include "twofold/hull.h"

namespace twofold {

namespace {

// Neighboring faces closer than this cosine merge into one polygon.
constexpr float kCoplanar = 0.9998f;

struct Face {
    int v[3];
    Vec3 normal;
    float offset = 0.0f;
    std::vector<int> outside;
    bool alive = true;
};

uint64_t edgeKey(int a, int b) { return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }

class Builder {
public:
    Builder(const std::vector<Vec3>& points, float tolerance) : p_(points), tol_(tolerance) {}

    void addFace(int a, int b, int c) {
        Face f;
        f.v[0] = a;
        f.v[1] = b;
        f.v[2] = c;
        f.normal = normalize(cross(p_[b] - p_[a], p_[c] - p_[a]));
        f.offset = dot(f.normal, p_[a]);
        const int fi = (int)faces_.size();
        for (int i = 0; i < 3; ++i) { edges_[edgeKey(f.v[i], f.v[(i + 1) % 3])] = fi; }
        faces_.push_back(std::move(f));
    }

    void assign(const std::vector<int>& points, const std::vector<int>& candidates) {
        for (int pi : points) {
            int best = -1;
            float far = tol_;
            for (int fi : candidates) {
                const float d = distance(faces_[fi], pi);
                if (d > far) { far = d; best = fi; }
            }
            if (best >= 0) { faces_[best].outside.push_back(pi); }
        }
    }

    void build(const int simplex[4], int maxVertices) {
        Vec3 center;
        for (int i = 0; i < 4; ++i) { center += p_[simplex[i]] * 0.25f; }
        const int tris[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 } };
        for (const auto& t : tris) {
            int a = simplex[t[0]], b = simplex[t[1]], c = simplex[t[2]];
            if (dot(cross(p_[b] - p_[a], p_[c] - p_[a]), p_[a] - center) < 0.0f) { std::swap(b, c); }
            addFace(a, b, c);
        }
        std::vector<int> rest;
        for (int i = 0; i < (int)p_.size(); ++i) {
            if (i != simplex[0] && i != simplex[1] && i != simplex[2] && i != simplex[3]) { rest.push_back(i); }
        }
        assign(rest, { 0, 1, 2, 3 });

        while (vertexCount() < maxVertices) {
            int face = -1, eye = -1;
            float far = 0.0f;
            for (int fi = 0; fi < (int)faces_.size(); ++fi) {
                if (!faces_[fi].alive) { continue; }
                for (int pi : faces_[fi].outside) {
                    const float d = distance(faces_[fi], pi);
                    if (d > far) { far = d; face = fi; eye = pi; }
                }
            }
            if (face < 0) { break; }
            addPoint(face, eye);
        }
    }

    HullMesh mesh() const {
        std::vector<int> region(faces_.size(), -1);
        std::vector<std::vector<int>> loops;
        for (int seed = 0; seed < (int)faces_.size(); ++seed) {
            if (!faces_[seed].alive || region[seed] >= 0) { continue; }
            const int id = (int)loops.size();
            std::vector<int> members{ seed };
            region[seed] = id;
            for (size_t k = 0; k < members.size(); ++k) {
                const Face& f = faces_[members[k]];
                for (int i = 0; i < 3; ++i) {
                    const int twin = edges_.at(edgeKey(f.v[(i + 1) % 3], f.v[i]));
                    if (region[twin] >= 0 || !coplanar(faces_[seed], faces_[twin])) { continue; }
                    region[twin] = id;
                    members.push_back(twin);
                }
            }
            std::vector<int> loop;
            if (!boundary(members, region, id, loop)) {
                for (size_t k = 1; k < members.size(); ++k) { region[members[k]] = -2; }
                const Face& f = faces_[seed];
                loop.assign(f.v, f.v + 3);
            }
            loops.push_back(loop);
        }
        for (int fi = 0; fi < (int)faces_.size(); ++fi) {
            if (region[fi] == -2) { loops.push_back({ faces_[fi].v[0], faces_[fi].v[1], faces_[fi].v[2] }); }
        }

        HullMesh out;
        std::vector<int> remap(p_.size(), -1);
        for (std::vector<int>& loop : loops) {
            for (int& v : loop) {
                if (remap[v] < 0) {
                    remap[v] = (int)out.vertices.size();
                    out.vertices.push_back(p_[v]);
                }
                v = remap[v];
            }
            out.faces.push_back(std::move(loop));
        }
        return out;
    }

private:
    float distance(const Face& f, int pi) const { return dot(f.normal, p_[pi]) - f.offset; }

    bool coplanar(const Face& seed, const Face& f) const {
        if (dot(seed.normal, f.normal) < kCoplanar) { return false; }
        for (int v : f.v) {
            if (std::fabs(distance(seed, v)) > tol_) { return false; }
        }
        return true;
    }

    int vertexCount() const {
        std::vector<bool> used(p_.size(), false);
        int n = 0;
        for (const Face& f : faces_) {
            if (!f.alive) { continue; }
            for (int v : f.v) {
                if (!used[v]) { used[v] = true; ++n; }
            }
        }
        return n;
    }

    void addPoint(int start, int eye) {
        std::vector<int> visible{ start };
        std::vector<bool> seen(faces_.size(), false);
        seen[start] = true;
        std::vector<std::pair<int, int>> horizon;
        for (size_t k = 0; k < visible.size(); ++k) {
            const Face& f = faces_[visible[k]];
            for (int i = 0; i < 3; ++i) {
                const int a = f.v[i], b = f.v[(i + 1) % 3];
                const int twin = edges_.at(edgeKey(b, a));
                if (seen[twin]) { continue; }
                if (distance(faces_[twin], eye) > tol_) {
                    seen[twin] = true;
                    visible.push_back(twin);
                } else {
                    horizon.push_back({ a, b });
                }
            }
        }

        std::vector<int> orphans;
        for (int fi : visible) {
            Face& f = faces_[fi];
            f.alive = false;
            for (int i = 0; i < 3; ++i) { edges_.erase(edgeKey(f.v[i], f.v[(i + 1) % 3])); }
            for (int pi : f.outside) {
                if (pi != eye) { orphans.push_back(pi); }
            }
            f.outside.clear();
        }
        std::vector<int> created;
        for (const auto& e : horizon) {
            created.push_back((int)faces_.size());
            addFace(e.first, e.second, eye);
        }
        assign(orphans, created);
    }

    // Walks the outline of a merged region. Fails unless it is a single simple loop.
    bool boundary(const std::vector<int>& members, const std::vector<int>& region, int id, std::vector<int>& loop) const {
        std::unordered_map<int, int> next;
        for (int fi : members) {
            const Face& f = faces_[fi];
            for (int i = 0; i < 3; ++i) {
                const int a = f.v[i], b = f.v[(i + 1) % 3];
                if (region[edges_.at(edgeKey(b, a))] == id) { continue; }
                if (!next.emplace(a, b).second) { return false; }
            }
        }
        if (next.empty()) { return false; }
        int start = next.begin()->first;
        for (const auto& e : next) { start = std::min(start, e.first); }
        int v = start;
        do {
            loop.push_back(v);
            auto it = next.find(v);
            if (it == next.end() || loop.size() > next.size()) { return false; }
            v = it->second;
        } while (v != start);
        return loop.size() == next.size();
    }

    const std::vector<Vec3>& p_;
    float tol_;
    std::vector<Face> faces_;
    std::unordered_map<uint64_t, int> edges_;
};

std::vector<Vec3> thicken(const std::vector<Vec3>& points, const Vec3* axes, int count, float t) {
    std::vector<Vec3> out;
    out.reserve(points.size() << count);
    for (const Vec3& p : points) {
        for (int corner = 0; corner < (1 << count); ++corner) {
            Vec3 q = p;
            for (int i = 0; i < count; ++i) { q += axes[i] * ((corner >> i) & 1 ? t : -t); }
            out.push_back(q);
        }
    }
    return out;
}

}

HullMesh computeHull(const std::vector<Vec3>& points, int maxVertices) {
    if (points.empty()) { return {}; }
    maxVertices = std::max(maxVertices, 4);

    Vec3 lo = points[0], hi = points[0];
    int extreme[6] = {};
    for (int i = 0; i < (int)points.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            if (points[i][k] < points[extreme[2 * k]][k]) { extreme[2 * k] = i; }
            if (points[i][k] > points[extreme[2 * k + 1]][k]) { extreme[2 * k + 1] = i; }
        }
        lo = minVec(lo, points[i]);
        hi = maxVec(hi, points[i]);
    }
    const float size = length(hi - lo);
    float reach = 0.0f;
    for (int k = 0; k < 3; ++k) { reach += std::max(std::fabs(lo[k]), std::fabs(hi[k])); }
    const float tol = std::max(1e-4f * size, 3.0f * FLT_EPSILON * reach);
    const float slab = std::max(0.005f * size, 1e-4f);

    int s[4] = { extreme[0], extreme[0], 0, 0 };
    float best = 0.0f;
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            const float d = lengthSquared(points[extreme[i]] - points[extreme[j]]);
            if (d > best) { best = d; s[0] = extreme[i]; s[1] = extreme[j]; }
        }
    }
    if (best <= tol * tol) {
        const Vec3 axes[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        return computeHull(thicken(points, axes, 3, slab), maxVertices);
    }

    const Vec3 dir = normalize(points[s[1]] - points[s[0]]);
    best = 0.0f;
    for (int i = 0; i < (int)points.size(); ++i) {
        const float d = length(cross(points[i] - points[s[0]], dir));
        if (d > best) { best = d; s[2] = i; }
    }
    if (best <= tol) {
        const Vec3 u = perpendicular(dir);
        const Vec3 axes[2] = { u, cross(dir, u) };
        return computeHull(thicken(points, axes, 2, slab), maxVertices);
    }

    const Vec3 n = normalize(cross(points[s[1]] - points[s[0]], points[s[2]] - points[s[0]]));
    best = 0.0f;
    for (int i = 0; i < (int)points.size(); ++i) {
        const float d = std::fabs(dot(points[i] - points[s[0]], n));
        if (d > best) { best = d; s[3] = i; }
    }
    if (best <= tol) { return computeHull(thicken(points, &n, 1, slab), maxVertices); }

    Builder builder(points, tol);
    builder.build(s, maxVertices);
    return builder.mesh();
}

}
