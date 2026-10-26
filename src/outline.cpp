#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>

#include "twofold/outline.h"

namespace twofold {

namespace {

// Corners whose turn has a smaller sine than this count as straight.
constexpr float kFlat = 1e-5f;

float lengthSquared(Vec2 v) { return dot(v, v); }

bool flat(Vec2 a, Vec2 b, Vec2 c) {
    const Vec2 e = b - a, f = c - b;
    return std::fabs(cross(e, f)) <= kFlat * std::sqrt(lengthSquared(e) * lengthSquared(f));
}

bool turnsLeft(Vec2 a, Vec2 b, Vec2 c) { return cross(b - a, c - b) > 0.0f && !flat(a, b, c); }

double signedArea(const std::vector<Vec2>& poly) {
    double a = 0.0;
    for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
        a += (double)poly[j].x * poly[i].y - (double)poly[i].x * poly[j].y;
    }
    return 0.5 * a;
}

std::vector<Vec2> removeFlat(std::vector<Vec2> poly) {
    bool changed = true;
    while (changed && poly.size() >= 3) {
        changed = false;
        std::vector<Vec2> out;
        for (size_t i = 0; i < poly.size(); ++i) {
            const Vec2 prev = out.empty() ? poly.back() : out.back();
            if (flat(prev, poly[i], poly[(i + 1) % poly.size()])) {
                changed = true;
                continue;
            }
            out.push_back(poly[i]);
        }
        poly.swap(out);
    }
    if (poly.size() < 3) { poly.clear(); }
    return poly;
}

float segmentDistanceSquared(Vec2 p, Vec2 a, Vec2 b) {
    const Vec2 d = b - a, r = p - a;
    const float len = lengthSquared(d);
    const float t = len > 0.0f ? std::min(std::max(dot(r, d) / len, 0.0f), 1.0f) : 0.0f;
    return lengthSquared(r - d * t);
}

// Line through the points strictly between first and last, which may run past the end of the loop.
void fitLine(const std::vector<Vec2>& loop, size_t first, size_t last, Vec2& origin, Vec2& dir) {
    const size_t n = loop.size();
    const Vec2 a = loop[first % n], chord = loop[last % n] - a;
    if (last - first < 3) {
        origin = a;
        dir = normalize(chord);
        return;
    }
    Vec2 c;
    for (size_t i = first + 1; i < last; ++i) { c += loop[i % n]; }
    c = c * (1.0f / (float)(last - first - 1));
    float xx = 0.0f, xy = 0.0f, yy = 0.0f;
    for (size_t i = first + 1; i < last; ++i) {
        const Vec2 r = loop[i % n] - c;
        xx += r.x * r.x;
        xy += r.x * r.y;
        yy += r.y * r.y;
    }
    const float angle = 0.5f * std::atan2(2.0f * xy, xx - yy);
    origin = c;
    dir = { std::cos(angle), std::sin(angle) };
    if (dot(dir, chord) < 0.0f) { dir = -dir; }
}

uint64_t pointKey(Vec2 p) {
    uint32_t x, y;
    std::memcpy(&x, &p.x, 4);
    std::memcpy(&y, &p.y, 4);
    return ((uint64_t)x << 32) | y;
}

uint64_t edgeKey(int a, int b) { return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; }

struct Contour {
    std::vector<Vec2> points;
    Vec2 lo, hi;
    int depth = 0;
    int parent = -1;
};

bool contains(const Contour& c, Vec2 p) {
    if (p.x < c.lo.x || p.x > c.hi.x || p.y < c.lo.y || p.y > c.hi.y) { return false; }
    bool inside = false;
    const std::vector<Vec2>& q = c.points;
    for (size_t i = 0, j = q.size() - 1; i < q.size(); j = i++) {
        if ((q[i].y > p.y) != (q[j].y > p.y) && p.x < q[i].x + (q[j].x - q[i].x) * (p.y - q[i].y) / (q[j].y - q[i].y)) {
            inside = !inside;
        }
    }
    return inside;
}

bool insideTriangle(Vec2 a, Vec2 b, Vec2 c, Vec2 p) {
    return cross(b - a, p - a) >= 0.0f && cross(c - b, p - b) >= 0.0f && cross(a - c, p - c) >= 0.0f;
}

// True when direction d from vertex i of a counter-clockwise polygon points into the interior.
bool opensToward(const std::vector<Vec2>& poly, size_t i, Vec2 d) {
    const size_t n = poly.size();
    const Vec2 in = poly[(i + n - 1) % n] - poly[i], out = poly[(i + 1) % n] - poly[i];
    if (cross(out, in) >= 0.0f) { return cross(out, d) > 0.0f && cross(d, in) > 0.0f; }
    return !(cross(in, d) >= 0.0f && cross(d, out) >= 0.0f);
}

// Joins a clockwise hole to a counter-clockwise outline through a zero-width cut to a visible vertex.
void bridge(std::vector<Vec2>& outer, const std::vector<Vec2>& hole) {
    size_t m = 0;
    for (size_t i = 1; i < hole.size(); ++i) {
        if (hole[i].x > hole[m].x || (hole[i].x == hole[m].x && hole[i].y < hole[m].y)) { m = i; }
    }
    const Vec2 mp = hole[m];
    const size_t n = outer.size();

    float hitX = INFINITY;
    size_t edge = n;
    for (size_t i = 0; i < n; ++i) {
        const Vec2 a = outer[i], b = outer[(i + 1) % n];
        if (a.y == b.y || (a.y > mp.y) == (b.y > mp.y)) { continue; }
        const float x = a.x + (b.x - a.x) * (mp.y - a.y) / (b.y - a.y);
        if (x >= mp.x && x < hitX) {
            hitX = x;
            edge = i;
        }
    }
    size_t p = 0;
    if (edge == n) {
        float best = INFINITY;
        for (size_t i = 0; i < n; ++i) {
            const float d = lengthSquared(outer[i] - mp);
            if (d < best) { best = d; p = i; }
        }
    } else {
        const Vec2 hit{ hitX, mp.y };
        const size_t a = edge, b = (edge + 1) % n;
        p = outer[a].x > outer[b].x ? a : b;
        const Vec2 pp = outer[p];
        float bestAngle = INFINITY, bestDist = INFINITY;
        for (size_t i = 0; i < n && pp.y != mp.y; ++i) {
            const Vec2 v = outer[i];
            if (i == p || (v.x == pp.x && v.y == pp.y)) { continue; }
            if (turnsLeft(outer[(i + n - 1) % n], v, outer[(i + 1) % n])) { continue; }
            const bool in = pp.y > mp.y ? insideTriangle(mp, hit, pp, v) : insideTriangle(mp, pp, hit, v);
            if (!in) { continue; }
            const Vec2 d = v - mp;
            const float angle = std::fabs(d.y) / std::sqrt(lengthSquared(d));
            const float dist = lengthSquared(d);
            if (angle < bestAngle || (angle == bestAngle && dist < bestDist)) {
                bestAngle = angle;
                bestDist = dist;
                p = i;
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        if (outer[i].x == outer[p].x && outer[i].y == outer[p].y && opensToward(outer, i, mp - outer[i])) {
            p = i;
            break;
        }
    }

    std::vector<Vec2> out;
    out.reserve(n + hole.size() + 2);
    out.insert(out.end(), outer.begin(), outer.begin() + (std::ptrdiff_t)p + 1);
    for (size_t k = 0; k <= hole.size(); ++k) { out.push_back(hole[(m + k) % hole.size()]); }
    out.insert(out.end(), outer.begin() + (std::ptrdiff_t)p, outer.end());
    outer.swap(out);
}

struct Triangle {
    int v[3];
};

void earClip(const std::vector<Vec2>& poly, std::vector<Triangle>& tris) {
    const int n = (int)poly.size();
    std::vector<int> prev(n), next(n);
    for (int i = 0; i < n; ++i) {
        prev[i] = (i + n - 1) % n;
        next[i] = (i + 1) % n;
    }
    auto same = [&](int i, int j) { return poly[i].x == poly[j].x && poly[i].y == poly[j].y; };
    auto remove = [&](int i) {
        next[prev[i]] = next[i];
        prev[next[i]] = prev[i];
    };
    auto isEar = [&](int i) {
        const int a = prev[i], c = next[i];
        for (int j = next[c]; j != a; j = next[j]) {
            if (same(j, a) || same(j, i) || same(j, c)) { continue; }
            if (insideTriangle(poly[a], poly[i], poly[c], poly[j])) { return false; }
        }
        return true;
    };

    int count = n, i = 0, stall = 0;
    while (count > 2) {
        const int a = prev[i], c = next[i];
        if (flat(poly[a], poly[i], poly[c])) {
            remove(i);
            --count;
            i = c;
            stall = 0;
            continue;
        }
        if (turnsLeft(poly[a], poly[i], poly[c]) && isEar(i)) {
            tris.push_back({ { a, i, c } });
            remove(i);
            --count;
            i = c;
            stall = 0;
            continue;
        }
        i = c;
        if (++stall <= count) { continue; }
        int best = i;
        float sharpest = -INFINITY;
        int j = i;
        do {
            const float t = cross(poly[j] - poly[prev[j]], poly[next[j]] - poly[j]);
            if (t > sharpest) { sharpest = t; best = j; }
            j = next[j];
        } while (j != i);
        if (sharpest > 0.0f) { tris.push_back({ { prev[best], best, next[best] } }); }
        i = next[best];
        remove(best);
        --count;
        stall = 0;
    }
}

// Hertel-Mehlhorn: drops diagonals, longest first, wherever both ends stay convex.
void mergeConvex(const std::vector<Vec2>& poly, const std::vector<Triangle>& tris, int maxVertices,
                 std::vector<std::vector<Vec2>>& pieces) {
    std::unordered_map<uint64_t, int> ids;
    std::vector<Vec2> points;
    std::vector<int> canon(poly.size());
    for (size_t i = 0; i < poly.size(); ++i) {
        auto it = ids.emplace(pointKey(poly[i]), (int)points.size());
        if (it.second) { points.push_back(poly[i]); }
        canon[i] = it.first->second;
    }

    std::vector<std::vector<int>> loops;
    std::unordered_map<uint64_t, int> owner;
    for (const Triangle& t : tris) {
        const int id = (int)loops.size();
        loops.push_back({ canon[t.v[0]], canon[t.v[1]], canon[t.v[2]] });
        for (int k = 0; k < 3; ++k) { owner[edgeKey(loops[id][k], loops[id][(k + 1) % 3])] = id; }
    }

    std::vector<std::pair<int, int>> diagonals;
    for (const std::vector<int>& l : loops) {
        for (int k = 0; k < 3; ++k) {
            const int a = l[k], b = l[(k + 1) % 3];
            if (a < b && owner.count(edgeKey(b, a))) { diagonals.push_back({ a, b }); }
        }
    }
    std::sort(diagonals.begin(), diagonals.end(), [&](const std::pair<int, int>& x, const std::pair<int, int>& y) {
        const float lx = lengthSquared(points[x.first] - points[x.second]);
        const float ly = lengthSquared(points[y.first] - points[y.second]);
        return lx > ly || (lx == ly && x < y);
    });

    std::vector<bool> alive(loops.size(), true);
    for (const std::pair<int, int>& d : diagonals) {
        const int a = d.first, b = d.second;
        auto ox = owner.find(edgeKey(a, b)), oy = owner.find(edgeKey(b, a));
        if (ox == owner.end() || oy == owner.end() || ox->second == oy->second) { continue; }
        const int x = ox->second, y = oy->second;
        const std::vector<int>& lx = loops[x];
        const std::vector<int>& ly = loops[y];
        const int nx = (int)lx.size(), ny = (int)ly.size();
        const int ia = (int)(std::find(lx.begin(), lx.end(), a) - lx.begin());
        const int ib = (int)(std::find(ly.begin(), ly.end(), b) - ly.begin());

        const Vec2 pa = points[a], pb = points[b];
        const Vec2 beforeA = points[lx[(ia + nx - 1) % nx]], afterA = points[ly[(ib + 2) % ny]];
        const Vec2 beforeB = points[ly[(ib + ny - 1) % ny]], afterB = points[lx[(ia + 2) % nx]];
        const bool flatA = flat(beforeA, pa, afterA), flatB = flat(beforeB, pb, afterB);
        if (!flatA && !turnsLeft(beforeA, pa, afterA)) { continue; }
        if (!flatB && !turnsLeft(beforeB, pb, afterB)) { continue; }
        if (nx + ny - 2 - (int)flatA - (int)flatB > maxVertices) { continue; }

        std::vector<int> merged;
        merged.reserve(nx + ny - 2);
        for (int k = 1; k <= nx; ++k) { merged.push_back(lx[(ia + k) % nx]); }
        for (int k = 2; k < ny; ++k) { merged.push_back(ly[(ib + k) % ny]); }
        owner.erase(edgeKey(a, b));
        owner.erase(edgeKey(b, a));
        for (size_t k = 0; k < merged.size(); ++k) {
            auto it = owner.find(edgeKey(merged[k], merged[(k + 1) % merged.size()]));
            if (it != owner.end()) { it->second = x; }
        }
        loops[x] = std::move(merged);
        loops[y].clear();
        alive[y] = false;
    }

    for (size_t i = 0; i < loops.size(); ++i) {
        if (!alive[i]) { continue; }
        std::vector<Vec2> piece;
        for (int v : loops[i]) { piece.push_back(points[v]); }
        piece = removeFlat(std::move(piece));
        if (!piece.empty() && signedArea(piece) > 0.0) { pieces.push_back(std::move(piece)); }
    }
}

}

std::vector<std::vector<Vec2>> traceAlpha(const AlphaImage& image, uint8_t threshold) {
    std::vector<std::vector<Vec2>> contours;
    if (!image.pixels || image.width <= 0 || image.height <= 0) { return contours; }
    const int w = image.width, h = image.height;
    const size_t px = (size_t)image.pixelStride;
    const size_t row = image.rowStride > 0 ? (size_t)image.rowStride : (size_t)w * px;
    const int inside = std::max<int>(threshold, 1);
    const float level = (float)inside - 0.5f;
    auto value = [&](int x, int y) -> float {
        if (x < 0 || y < 0 || x >= w || y >= h) { return 0.0f; }
        return (float)image.pixels[(size_t)y * row + (size_t)x * px];
    };
    const size_t cols = (size_t)w + 2;
    // Edges of the sample grid: the horizontal one to the right of a sample, then the vertical one below it.
    auto edgeId = [&](int x, int y, int k) {
        if (k == 1) { ++x; }
        if (k == 2) { ++y; }
        return ((size_t)(y + 1) * cols + (size_t)(x + 1)) * 2 + (k & 1);
    };
    std::vector<uint8_t> visited(cols * ((size_t)h + 2) * 2, 0);

    // Cell corners run clockwise from the top left and edge k joins corner k to corner k + 1. A contour
    // enters a cell on an edge that leaves the filled side and exits on the next edge that returns to it.
    auto follow = [&](int x, int y, int k) {
        std::vector<Vec2> loop;
        for (size_t id = edgeId(x, y, k); !visited[id]; id = edgeId(x, y, k)) {
            visited[id] = 1;
            const float c[4] = { value(x, y), value(x + 1, y), value(x + 1, y + 1), value(x, y + 1) };
            int mask = 0;
            for (int i = 0; i < 4; ++i) { mask |= (c[i] > level ? 1 : 0) << i; }
            const Vec2 corner[4] = { { x + 0.5f, y + 0.5f }, { x + 1.5f, y + 0.5f }, { x + 1.5f, y + 1.5f }, { x + 0.5f, y + 1.5f } };
            const int k1 = (k + 1) & 3;
            loop.push_back(corner[k] + (corner[k1] - corner[k]) * ((level - c[k]) / (c[k1] - c[k])));
            int e = k1;
            if (mask == 5 || mask == 10) {
                if (c[0] + c[1] + c[2] + c[3] <= 4.0f * level) { e = (k + 3) & 3; }
            } else {
                while ((mask >> e & 1) || !(mask >> ((e + 1) & 3) & 1)) { e = (e + 1) & 3; }
            }
            static const int dx[4] = { 0, 1, 0, -1 }, dy[4] = { -1, 0, 1, 0 };
            x += dx[e];
            y += dy[e];
            k = (e + 2) & 3;
        }
        if (loop.size() >= 3) { contours.push_back(std::move(loop)); }
    };

    std::vector<uint8_t> above(cols, 0), below(cols, 0);
    int aboveCount = 0, belowCount = 0;
    for (int y = -1; y < h; ++y) {
        above.swap(below);
        aboveCount = belowCount;
        belowCount = 0;
        if (y + 1 < h) {
            const uint8_t* p = image.pixels + (size_t)(y + 1) * row;
            for (int x = 0; x < w; ++x, p += px) {
                below[x + 1] = *p >= inside ? 1 : 0;
                belowCount += below[x + 1];
            }
        } else {
            std::fill(below.begin(), below.end(), (uint8_t)0);
        }
        if (aboveCount == 0 && belowCount == 0) { continue; }
        for (int x = -1; x < w; ++x) {
            const int mask = above[x + 1] | above[x + 2] << 1 | below[x + 2] << 2 | below[x + 1] << 3;
            if (mask == 0 || mask == 15) { continue; }
            for (int k = 0; k < 4; ++k) {
                if ((mask >> k & 1) && !(mask >> ((k + 1) & 3) & 1) && !visited[edgeId(x, y, k)]) { follow(x, y, k); }
            }
        }
    }
    return contours;
}

namespace {

// Douglas-Peucker on a closed contour, with corners moved to where fitted straight runs meet.
// Edges that cross other edges can be refined one at a time until the result is simple again.
class Simplifier {
public:
    Simplifier(const std::vector<Vec2>& contour, float tolerance)
        : c_(contour), n_(contour.size()), tolSq_(tolerance * tolerance), keep_(n_, false), pinned_(n_, false) {
        if (n_ < 4) {
            out_ = c_;
            return;
        }
        size_t far = 0;
        float best = -1.0f;
        for (size_t i = 1; i < n_; ++i) {
            const float d = lengthSquared(c_[i] - c_[0]);
            if (d > best) { best = d; far = i; }
        }
        keep_[0] = keep_[far] = true;
        std::vector<std::pair<size_t, size_t>> stack{ { 0, far }, { far, n_ } };
        while (!stack.empty()) {
            const size_t a = stack.back().first, b = stack.back().second;
            stack.pop_back();
            const size_t worst = farthest(a, b, tolSq_);
            if (worst == a) { continue; }
            keep_[worst] = true;
            stack.push_back({ a, worst });
            stack.push_back({ worst, b });
        }
        build();
    }

    const std::vector<Vec2>& points() const { return out_; }

    // Brings the edge from output vertex v to the next one closer to the contour. False once it
    // already follows the contour exactly.
    bool refine(size_t v) {
        if (corners_.empty()) { return false; }
        const size_t g = corners_.size();
        const Corner& a = corners_[v];
        const Corner& b = corners_[(v + 1) % g];
        bool pinned = false;
        for (const Corner* c : { &a, &b }) {
            if (c->first == c->last && !moved(*c)) { continue; }
            for (size_t slot = c->first; slot <= c->last; ++slot) { pinned_[kept_[slot % m_]] = true; }
            pinned = true;
        }
        if (!pinned) {
            const size_t last = v + 1 < g ? b.first : b.first + m_;
            const size_t worst = farthest(index(a.last), index(last), -1.0f);
            if (worst == index(a.last)) { return false; }
            keep_[worst % n_] = true;
        }
        build();
        return true;
    }

private:
    struct Corner {
        size_t first, last;
        Vec2 at;
    };

    // Contour point strictly between a and b farthest from their chord, or a when none is beyond floor.
    size_t farthest(size_t a, size_t b, float floor) const {
        const Vec2 pa = c_[a % n_], pb = c_[b % n_];
        size_t worst = a;
        for (size_t i = a + 1; i < b; ++i) {
            const float e = segmentDistanceSquared(c_[i % n_], pa, pb);
            if (e > floor) { floor = e; worst = i; }
        }
        return worst;
    }

    // Slots count kept vertices around the loop and may run past m, so contour indices keep increasing.
    size_t index(size_t slot) const { return kept_[slot % m_] + slot / m_ * n_; }
    Vec2 home(size_t slot) const { return c_[kept_[slot % m_]]; }
    bool moved(const Corner& c) const { return c.at.x != home(c.first).x || c.at.y != home(c.first).y; }

    void build() {
        kept_.clear();
        for (size_t i = 0; i < n_; ++i) {
            if (keep_[i]) { kept_.push_back(i); }
        }
        m_ = kept_.size();
        corners_.clear();
        out_.clear();
        if (m_ < 3) {
            for (size_t k : kept_) { out_.push_back(c_[k]); }
            return;
        }

        std::vector<Vec2> origin(m_), dir(m_);
        std::vector<bool> fitted(m_);
        for (size_t j = 0; j < m_; ++j) {
            fitLine(c_, index(j), index(j + 1), origin[j], dir[j]);
            fitted[j] = index(j + 1) - index(j) >= 3;
        }
        size_t start = 0;
        while (start < m_ && !fitted[(start + m_ - 1) % m_]) { ++start; }
        if (start == m_) { start = 0; }

        // Up to two short spans between fitted runs, such as a corner cut by tracing, collapse into the
        // point where the runs meet.
        for (size_t j = start; j < start + m_;) {
            const size_t in = (j + m_ - 1) % m_;
            size_t k = 0;
            while (k < 3 && !fitted[(j + k) % m_]) { ++k; }
            Corner c{ j, j, home(j) };
            const size_t out = (j + k) % m_;
            const float s = cross(dir[in], dir[out]);
            bool free = fitted[in] && k < 3 && j + k < start + m_ && std::fabs(s) > 1e-6f;
            for (size_t slot = j; free && slot <= j + k; ++slot) { free = !pinned_[kept_[slot % m_]]; }
            if (free) {
                const Vec2 q = origin[in] + dir[in] * (cross(origin[out] - origin[in], dir[out]) / s);
                bool near = true;
                for (size_t slot = j; slot <= j + k; ++slot) { near = near && lengthSquared(q - home(slot)) <= 4.0f * tolSq_; }
                if (near) { c = { j, j + k, q }; }
            }
            corners_.push_back(c);
            j = c.last + 1;
        }

        const float limit = tolSq_ * 1.0001f + 1e-12f;
        for (bool split = true; split;) {
            split = false;
            const size_t g = corners_.size();
            for (size_t c = 0; c < g && !split; ++c) {
                const Corner& a = corners_[c];
                const Corner& prev = corners_[(c + g - 1) % g];
                const Corner& next = corners_[(c + 1) % g];
                const size_t nextFirst = c + 1 < g ? next.first : next.first + m_;
                bool bad = false;
                for (size_t i = index(a.first) + 1; i < index(a.last) && !bad; ++i) {
                    bad = std::min(segmentDistanceSquared(c_[i % n_], prev.at, a.at),
                                   segmentDistanceSquared(c_[i % n_], a.at, next.at)) > limit;
                }
                for (size_t i = index(a.last); i <= index(nextFirst) && !bad; ++i) {
                    bad = segmentDistanceSquared(c_[i % n_], a.at, next.at) > limit;
                }
                if (!bad) { continue; }
                for (size_t r : { std::max(c, (c + 1) % g), std::min(c, (c + 1) % g) }) {
                    const Corner old = corners_[r];
                    if (old.first == old.last && !moved(old)) { continue; }
                    std::vector<Corner> singles;
                    for (size_t slot = old.first; slot <= old.last; ++slot) { singles.push_back({ slot, slot, home(slot) }); }
                    corners_.erase(corners_.begin() + (std::ptrdiff_t)r);
                    corners_.insert(corners_.begin() + (std::ptrdiff_t)r, singles.begin(), singles.end());
                    split = true;
                }
            }
        }
        for (const Corner& c : corners_) { out_.push_back(c.at); }
    }

    const std::vector<Vec2>& c_;
    size_t n_;
    float tolSq_;
    std::vector<bool> keep_, pinned_;
    std::vector<size_t> kept_;
    size_t m_ = 0;
    std::vector<Corner> corners_;
    std::vector<Vec2> out_;
};

bool segmentsMeet(Vec2 p, Vec2 q, Vec2 r, Vec2 s) {
    const float d1 = cross(q - p, r - p), d2 = cross(q - p, s - p);
    const float d3 = cross(s - r, p - r), d4 = cross(s - r, q - r);
    if (((d1 > 0.0f && d2 < 0.0f) || (d1 < 0.0f && d2 > 0.0f)) && ((d3 > 0.0f && d4 < 0.0f) || (d3 < 0.0f && d4 > 0.0f))) {
        return true;
    }
    auto on = [](Vec2 a, Vec2 b, Vec2 c, float d) {
        return d == 0.0f && std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= c.y &&
               c.y <= std::max(a.y, b.y);
    };
    return on(p, q, r, d1) || on(p, q, s, d2) || on(r, s, p, d3) || on(r, s, q, d4);
}

// Simplifies every contour, then refines edges that touch or cross until no two edges meet.
std::vector<std::vector<Vec2>> simplifyTogether(const std::vector<std::vector<Vec2>>& contours, float tolerance) {
    std::vector<Simplifier> parts;
    parts.reserve(contours.size());
    for (const std::vector<Vec2>& c : contours) { parts.emplace_back(c, tolerance); }

    struct Edge {
        uint32_t contour, index;
        Vec2 a, b;
        float lo, hi;
    };
    for (bool progress = true; progress;) {
        std::vector<Edge> edges;
        for (uint32_t ci = 0; ci < parts.size(); ++ci) {
            const std::vector<Vec2>& p = parts[ci].points();
            for (uint32_t i = 0; i < p.size(); ++i) {
                const Vec2 a = p[i], b = p[(i + 1) % p.size()];
                edges.push_back({ ci, i, a, b, std::min(a.x, b.x), std::max(a.x, b.x) });
            }
        }
        std::sort(edges.begin(), edges.end(), [](const Edge& x, const Edge& y) {
            return x.lo < y.lo || (x.lo == y.lo && (x.contour < y.contour || (x.contour == y.contour && x.index < y.index)));
        });
        std::vector<std::pair<uint32_t, uint32_t>> crossing;
        for (size_t i = 0; i < edges.size(); ++i) {
            const Edge& e = edges[i];
            for (size_t j = i + 1; j < edges.size() && edges[j].lo <= e.hi; ++j) {
                const Edge& f = edges[j];
                if (std::max(e.a.y, e.b.y) < std::min(f.a.y, f.b.y) || std::max(f.a.y, f.b.y) < std::min(e.a.y, e.b.y)) {
                    continue;
                }
                if (e.contour == f.contour) {
                    const uint32_t size = (uint32_t)parts[e.contour].points().size();
                    if ((e.index + 1) % size == f.index || (f.index + 1) % size == e.index) { continue; }
                }
                if (!segmentsMeet(e.a, e.b, f.a, f.b)) { continue; }
                crossing.push_back({ e.contour, e.index });
                crossing.push_back({ f.contour, f.index });
            }
        }
        std::sort(crossing.begin(), crossing.end());
        crossing.erase(std::unique(crossing.begin(), crossing.end()), crossing.end());
        // Refining changes a contour's edge numbering, so each contour refines one edge per pass.
        progress = false;
        for (size_t k = crossing.size(); k-- > 0;) {
            const uint32_t ci = crossing[k].first;
            if (!parts[ci].refine(crossing[k].second)) { continue; }
            progress = true;
            while (k > 0 && crossing[k - 1].first == ci) { --k; }
        }
    }
    std::vector<std::vector<Vec2>> out;
    for (const Simplifier& s : parts) { out.push_back(s.points()); }
    return out;
}

}

std::vector<Vec2> simplifyContour(const std::vector<Vec2>& contour, float tolerance) {
    return simplifyTogether({ contour }, tolerance)[0];
}

std::vector<std::vector<Vec2>> convexDecompose(const std::vector<std::vector<Vec2>>& contours, int maxVertices) {
    maxVertices = std::max(maxVertices, 3);
    std::vector<Contour> cs;
    for (const std::vector<Vec2>& in : contours) {
        Contour c;
        c.points = removeFlat(in);
        if (c.points.empty() || signedArea(c.points) == 0.0) { continue; }
        c.lo = c.hi = c.points[0];
        for (const Vec2& p : c.points) {
            c.lo = { std::min(c.lo.x, p.x), std::min(c.lo.y, p.y) };
            c.hi = { std::max(c.hi.x, p.x), std::max(c.hi.y, p.y) };
        }
        cs.push_back(std::move(c));
    }

    const int count = (int)cs.size();
    std::vector<std::vector<int>> containers(count);
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < count; ++j) {
            if (i != j && contains(cs[j], cs[i].points[0])) { containers[i].push_back(j); }
        }
        cs[i].depth = (int)containers[i].size();
    }
    for (int i = 0; i < count; ++i) {
        if (cs[i].depth % 2 == 0) { continue; }
        for (int j : containers[i]) {
            if (cs[j].depth == cs[i].depth - 1) { cs[i].parent = j; }
        }
    }

    std::vector<std::vector<Vec2>> pieces;
    for (int i = 0; i < count; ++i) {
        if (cs[i].depth % 2 != 0) { continue; }
        std::vector<Vec2> poly = cs[i].points;
        if (signedArea(poly) < 0.0) { std::reverse(poly.begin(), poly.end()); }
        std::vector<std::vector<Vec2>> holes;
        for (int j = 0; j < count; ++j) {
            if (cs[j].parent != i) { continue; }
            holes.push_back(cs[j].points);
            if (signedArea(holes.back()) > 0.0) { std::reverse(holes.back().begin(), holes.back().end()); }
        }
        std::stable_sort(holes.begin(), holes.end(), [](const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
            float xa = a[0].x, xb = b[0].x;
            for (const Vec2& p : a) { xa = std::max(xa, p.x); }
            for (const Vec2& p : b) { xb = std::max(xb, p.x); }
            return xa > xb;
        });
        for (const std::vector<Vec2>& hole : holes) { bridge(poly, hole); }

        std::vector<Triangle> tris;
        earClip(poly, tris);
        mergeConvex(poly, tris, maxVertices, pieces);
    }
    return pieces;
}

std::vector<std::vector<Vec2>> alphaPieces(const AlphaImage& image, uint8_t threshold, float tolerance, int maxVertices) {
    return convexDecompose(simplifyTogether(traceAlpha(image, threshold), tolerance), maxVertices);
}

}
