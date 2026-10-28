#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iterator>
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
    Builder(const std::vector<Vec3>& points, float tolerance, float size)
        : p_(points), tol_(tolerance), fold_(tolerance / size), reclaimed_(points.size(), false) {}

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
            float far = tol_, top = -FLT_MAX;
            for (int fi : candidates) {
                const float d = distance(faces_[fi], pi);
                top = std::max(top, d);
                if (d > far) { far = d; best = fi; }
            }
            if (best >= 0) {
                faces_[best].outside.push_back(pi);
            } else if (top > 0.0f && !reclaimed_[pi]) {
                near_.push_back(pi);
            }
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
            if (face < 0) {
                if (reclaim()) { continue; }
                break;
            }
            if (!addPoint(face, eye)) {
                std::vector<int>& out = faces_[face].outside;
                out.erase(std::find(out.begin(), out.end(), eye));
            }
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
            if (!boundary(members, region, id, loop) || !convex(loop, faces_[seed].normal)) {
                for (size_t k = 1; k < members.size(); ++k) { region[members[k]] = -1; }
                loop = grow(seed, region, id);
            }
            loops.push_back(loop);
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

    // A point dropped within tolerance of the surface can still end up outside the finished hull. Each point is
    // tested again once.
    bool reclaim() {
        std::vector<int> points, alive;
        points.swap(near_);
        for (int pi : points) { reclaimed_[pi] = true; }
        for (int fi = 0; fi < (int)faces_.size(); ++fi) {
            if (faces_[fi].alive) { alive.push_back(fi); }
        }
        assign(points, alive);
        return std::any_of(alive.begin(), alive.end(), [&](int fi) { return !faces_[fi].outside.empty(); });
    }

    // Whether the face across edge a-b of a visible face goes too. It does if the eye is clearly above it, if
    // the new face over a-b would fold back over it, or if that new face would be a sliver or face inward.
    bool sees(const Face& from, const Face& f, int a, int b, int eye) const {
        const float d = distance(f, eye);
        if (d > tol_) { return true; }
        const Vec3 edge = normalize(p_[b] - p_[a]);
        const Vec3 r = p_[eye] - p_[a];
        const float h = length(r - edge * dot(r, edge));
        if (h <= tol_) { return true; }
        if (distance(from, eye) <= tol_ && dot(r, cross(from.normal, edge)) <= 0.0f) { return true; }
        if (d <= 0.0f) { return false; }
        return dot(r, cross(f.normal, edge)) <= 0.0f || d > fold_ * h;
    }

    // Fails, changing nothing, unless the horizon is a single loop.
    bool addPoint(int start, int eye) {
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
                if (sees(f, faces_[twin], a, b, eye)) {
                    seen[twin] = true;
                    visible.push_back(twin);
                } else {
                    horizon.push_back({ a, b });
                }
            }
        }
        horizon.erase(std::remove_if(horizon.begin(), horizon.end(),
                                     [&](const std::pair<int, int>& e) { return seen[edges_.at(edgeKey(e.second, e.first))]; }),
                      horizon.end());
        if (!singleLoop(horizon)) { return false; }

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
        return true;
    }

    static bool singleLoop(const std::vector<std::pair<int, int>>& horizon) {
        for (size_t i = 0; i < horizon.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (horizon[i].first == horizon[j].first) { return false; }
            }
        }
        int v = horizon.empty() ? -1 : horizon[0].second;
        for (size_t steps = 1; steps < horizon.size(); ++steps) {
            size_t k = 0;
            while (k < horizon.size() && horizon[k].first != v) { ++k; }
            if (k == horizon.size() || k == 0) { return false; }
            v = horizon[k].second;
        }
        return !horizon.empty() && v == horizon[0].first;
    }

    // Corners may bend back by no more than the tolerance.
    bool convex(const std::vector<int>& loop, const Vec3& normal) const {
        const size_t n = loop.size();
        for (size_t i = 0; i < n; ++i) {
            const Vec3& a = p_[loop[i]];
            const Vec3& b = p_[loop[(i + 1) % n]];
            const Vec3& c = p_[loop[(i + 2) % n]];
            if (dot(cross(b - a, c - b), normal) < -tol_ * length(c - a)) { return false; }
        }
        return true;
    }

    // Grows the region one face at a time, keeping its outline a single convex loop.
    std::vector<int> grow(int seed, std::vector<int>& region, int id) const {
        std::vector<int> members{ seed }, loop;
        boundary(members, region, id, loop);
        for (size_t k = 0; k < members.size(); ++k) {
            for (int i = 0; i < 3; ++i) {
                const Face& f = faces_[members[k]];
                const int twin = edges_.at(edgeKey(f.v[(i + 1) % 3], f.v[i]));
                if (region[twin] >= 0 || !coplanar(faces_[seed], faces_[twin])) { continue; }
                region[twin] = id;
                members.push_back(twin);
                std::vector<int> grown;
                if (boundary(members, region, id, grown) && convex(grown, faces_[seed].normal)) {
                    loop.swap(grown);
                } else {
                    region[twin] = -1;
                    members.pop_back();
                }
            }
        }
        return loop;
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
    // Sine of the steepest fold back that leaves no corner more than the tolerance outside.
    float fold_;
    std::vector<Face> faces_;
    std::unordered_map<uint64_t, int> edges_;
    std::vector<int> near_;
    std::vector<bool> reclaimed_;
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

// A degenerate cloud is thickened once. The slab is well over the tolerance, so the thickened cloud never is.
HullMesh hullOf(const std::vector<Vec3>& points, int maxVertices, bool thickenFlat) {
    if (points.empty()) { return {}; }

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
    if (!std::isfinite(size) || !std::isfinite(reach)) { return {}; }
    const float tol = std::max({ 1e-4f * size, 3.0f * FLT_EPSILON * reach, 1e-8f });
    const float slab = std::max({ 0.005f * size, 1e-4f, 10.0f * tol });

    int s[4] = { extreme[0], extreme[0], 0, 0 };
    float best = 0.0f;
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            const float d = lengthSquared(points[extreme[i]] - points[extreme[j]]);
            if (d > best) { best = d; s[0] = extreme[i]; s[1] = extreme[j]; }
        }
    }
    if (best <= tol * tol) {
        if (!thickenFlat) { return {}; }
        const Vec3 axes[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        return hullOf(thicken(points, axes, 3, slab), maxVertices, false);
    }

    const Vec3 dir = normalize(points[s[1]] - points[s[0]]);
    best = 0.0f;
    for (int i = 0; i < (int)points.size(); ++i) {
        const float d = length(cross(points[i] - points[s[0]], dir));
        if (d > best) { best = d; s[2] = i; }
    }
    if (best <= tol) {
        if (!thickenFlat) { return {}; }
        const Vec3 u = perpendicular(dir);
        const Vec3 axes[2] = { u, cross(dir, u) };
        return hullOf(thicken(points, axes, 2, slab), maxVertices, false);
    }

    const Vec3 n = normalize(cross(points[s[1]] - points[s[0]], points[s[2]] - points[s[0]]));
    best = 0.0f;
    for (int i = 0; i < (int)points.size(); ++i) {
        const float d = std::fabs(dot(points[i] - points[s[0]], n));
        if (d > best) { best = d; s[3] = i; }
    }
    if (best <= tol) {
        if (!thickenFlat) { return {}; }
        return hullOf(thicken(points, &n, 1, slab), maxVertices, false);
    }

    Builder builder(points, tol, size);
    builder.build(s, maxVertices);
    return builder.mesh();
}

bool finite(const Vec3& p) { return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z); }

}

HullMesh computeHull(const std::vector<Vec3>& points, int maxVertices) {
    maxVertices = std::max(maxVertices, 4);
    if (std::all_of(points.begin(), points.end(), finite)) { return hullOf(points, maxVertices, true); }
    std::vector<Vec3> kept;
    std::copy_if(points.begin(), points.end(), std::back_inserter(kept), finite);
    return hullOf(kept, maxVertices, true);
}

}
