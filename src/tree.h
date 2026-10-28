#pragma once

#include <cstdint>
#include <vector>

#include "twofold/math.h"

namespace twofold {

// Bounding volume tree kept balanced by rotations. Each leaf holds an enlarged box, so it only moves
// when its item leaves the box.
class BoundsTree {
public:
    int insert(const Vec3& lo, const Vec3& hi, uint32_t item);
    void remove(int leaf);
    // Refits the leaf to lo, hi grown by margin if they left its box, and returns whether it did.
    bool update(int leaf, const Vec3& lo, const Vec3& hi, float margin);
    bool overlap(int a, int b) const { return overlaps(nodes_[a], nodes_[b].lo, nodes_[b].hi); }

    // Calls visit(item) for every other leaf whose box overlaps this leaf's.
    template <class Visit>
    void query(int leaf, Visit&& visit) const {
        const Node& n = nodes_[leaf];
        query(n.lo, n.hi, [&](uint32_t item) {
            if (item != n.item) { visit(item); }
        });
    }

    // Calls visit(item) for every leaf whose box overlaps lo, hi.
    template <class Visit>
    void query(const Vec3& lo, const Vec3& hi, Visit&& visit) const {
        traverse([&](const Node& n) { return overlaps(n, lo, hi); }, visit);
    }

    // Calls visit(item) for every leaf whose box the segment origin + t * d, t in [0, 1], passes through.
    template <class Visit>
    void raycast(const Vec3& origin, const Vec3& d, Visit&& visit) const {
        traverse([&](const Node& n) { return segmentHits(n, origin, d); }, visit);
    }

    // Calls visit(item) for every leaf whose box, grown by extent on every side, the segment passes through.
    template <class Visit>
    void boxcast(const Vec3& origin, const Vec3& d, const Vec3& extent, Visit&& visit) const {
        traverse([&](const Node& n) { return segmentHits({ n.lo - extent, n.hi + extent }, origin, d); }, visit);
    }

private:
    struct Node {
        Vec3 lo, hi;
        int parent = -1;
        int left = -1;
        int right = -1;
        int height = 0;
        uint32_t item = 0;
    };

    std::vector<Node> nodes_;
    int root_ = -1;
    int free_ = -1;

    static bool overlaps(const Node& n, const Vec3& lo, const Vec3& hi) {
        return n.lo.x <= hi.x && lo.x <= n.hi.x && n.lo.y <= hi.y && lo.y <= n.hi.y && n.lo.z <= hi.z &&
               lo.z <= n.hi.z;
    }
    static bool segmentHits(const Node& n, const Vec3& origin, const Vec3& d);

    template <class Test, class Visit>
    void traverse(Test&& test, Visit&& visit) const {
        if (root_ < 0) { return; }
        // A balanced tree of any size that fits in memory is far shallower than this.
        int stack[128];
        int count = 0;
        stack[count++] = root_;
        while (count > 0) {
            const Node& n = nodes_[stack[--count]];
            if (!test(n)) { continue; }
            if (n.left < 0) {
                visit(n.item);
            } else {
                stack[count++] = n.left;
                stack[count++] = n.right;
            }
        }
    }

    int allocate();
    void release(int i);
    void insertLeaf(int leaf);
    void removeLeaf(int leaf);
    void refit(int i);
    int balance(int a);
};

}
