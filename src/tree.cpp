#include "tree.h"

#include <algorithm>
#include <cmath>

namespace twofold {

namespace {

float area(const Vec3& lo, const Vec3& hi) {
    const Vec3 e = hi - lo;
    return e.x * e.y + e.y * e.z + e.z * e.x;
}

}

bool BoundsTree::segmentHits(const Node& n, const Vec3& origin, const Vec3& d) {
    float t0 = 0.0f, t1 = 1.0f;
    for (int k = 0; k < 3; ++k) {
        if (d[k] == 0.0f) {
            if (origin[k] < n.lo[k] || origin[k] > n.hi[k]) { return false; }
            continue;
        }
        const float inv = 1.0f / d[k];
        float ta = (n.lo[k] - origin[k]) * inv, tb = (n.hi[k] - origin[k]) * inv;
        if (ta > tb) { std::swap(ta, tb); }
        t0 = std::fmax(t0, ta);
        t1 = std::fmin(t1, tb);
        if (t0 > t1) { return false; }
    }
    return true;
}

int BoundsTree::allocate() {
    if (free_ < 0) {
        nodes_.emplace_back();
        return (int)nodes_.size() - 1;
    }
    const int i = free_;
    free_ = nodes_[i].parent;
    nodes_[i] = Node{};
    return i;
}

void BoundsTree::release(int i) {
    nodes_[i].parent = free_;
    nodes_[i].height = -1;
    free_ = i;
}

int BoundsTree::insert(const Vec3& lo, const Vec3& hi, uint32_t item) {
    const int leaf = allocate();
    nodes_[leaf].lo = lo;
    nodes_[leaf].hi = hi;
    nodes_[leaf].item = item;
    insertLeaf(leaf);
    return leaf;
}

void BoundsTree::remove(int leaf) {
    removeLeaf(leaf);
    release(leaf);
}

bool BoundsTree::update(int leaf, const Vec3& lo, const Vec3& hi, float margin) {
    Node& n = nodes_[leaf];
    if (n.lo.x <= lo.x && n.lo.y <= lo.y && n.lo.z <= lo.z && hi.x <= n.hi.x && hi.y <= n.hi.y && hi.z <= n.hi.z) {
        return false;
    }
    removeLeaf(leaf);
    const Vec3 g{ margin, margin, margin };
    nodes_[leaf].lo = lo - g;
    nodes_[leaf].hi = hi + g;
    insertLeaf(leaf);
    return true;
}

// Walks down to the sibling that adds the least area to the tree, as in Box2D.
void BoundsTree::insertLeaf(int leaf) {
    if (root_ < 0) {
        root_ = leaf;
        nodes_[leaf].parent = -1;
        return;
    }
    const Vec3 lo = nodes_[leaf].lo, hi = nodes_[leaf].hi;
    int index = root_;
    while (nodes_[index].left >= 0) {
        const Node& n = nodes_[index];
        const float here = area(n.lo, n.hi);
        const float combined = area(minVec(n.lo, lo), maxVec(n.hi, hi));
        const float inherited = 2.0f * (combined - here);
        float cost[2];
        const int child[2] = { n.left, n.right };
        for (int k = 0; k < 2; ++k) {
            const Node& c = nodes_[child[k]];
            const float grown = area(minVec(c.lo, lo), maxVec(c.hi, hi));
            cost[k] = (c.left < 0 ? grown : grown - area(c.lo, c.hi)) + inherited;
        }
        if (2.0f * combined < cost[0] && 2.0f * combined < cost[1]) { break; }
        index = cost[0] < cost[1] ? child[0] : child[1];
    }

    const int sibling = index;
    const int oldParent = nodes_[sibling].parent;
    const int parent = allocate();
    nodes_[parent].parent = oldParent;
    nodes_[parent].left = sibling;
    nodes_[parent].right = leaf;
    nodes_[sibling].parent = parent;
    nodes_[leaf].parent = parent;
    if (oldParent < 0) {
        root_ = parent;
    } else if (nodes_[oldParent].left == sibling) {
        nodes_[oldParent].left = parent;
    } else {
        nodes_[oldParent].right = parent;
    }

    for (int i = parent; i >= 0; i = nodes_[i].parent) {
        i = balance(i);
        refit(i);
    }
}

void BoundsTree::removeLeaf(int leaf) {
    if (leaf == root_) {
        root_ = -1;
        return;
    }
    const int parent = nodes_[leaf].parent;
    const int grand = nodes_[parent].parent;
    const int sibling = nodes_[parent].left == leaf ? nodes_[parent].right : nodes_[parent].left;
    release(parent);
    if (grand < 0) {
        root_ = sibling;
        nodes_[sibling].parent = -1;
        return;
    }
    if (nodes_[grand].left == parent) {
        nodes_[grand].left = sibling;
    } else {
        nodes_[grand].right = sibling;
    }
    nodes_[sibling].parent = grand;
    for (int i = grand; i >= 0; i = nodes_[i].parent) {
        i = balance(i);
        refit(i);
    }
}

void BoundsTree::refit(int i) {
    Node& n = nodes_[i];
    const Node& l = nodes_[n.left];
    const Node& r = nodes_[n.right];
    n.lo = minVec(l.lo, r.lo);
    n.hi = maxVec(l.hi, r.hi);
    n.height = 1 + std::max(l.height, r.height);
}

// Rotates the taller child up when the heights of a's children differ by more than one. Returns the
// node now in a's place.
int BoundsTree::balance(int a) {
    if (nodes_[a].left < 0 || nodes_[a].height < 2) { return a; }
    const int b = nodes_[a].left, c = nodes_[a].right;
    const int diff = nodes_[c].height - nodes_[b].height;
    if (diff >= -1 && diff <= 1) { return a; }

    const bool rightUp = diff > 1;
    const int up = rightUp ? c : b;
    const int f = nodes_[up].left, g = nodes_[up].right;

    nodes_[up].left = a;
    nodes_[up].parent = nodes_[a].parent;
    nodes_[a].parent = up;
    const int p = nodes_[up].parent;
    if (p < 0) {
        root_ = up;
    } else if (nodes_[p].left == a) {
        nodes_[p].left = up;
    } else {
        nodes_[p].right = up;
    }

    const bool keepF = nodes_[f].height > nodes_[g].height;
    const int kept = keepF ? f : g, moved = keepF ? g : f;
    nodes_[up].right = kept;
    if (rightUp) {
        nodes_[a].right = moved;
    } else {
        nodes_[a].left = moved;
    }
    nodes_[moved].parent = a;
    refit(a);
    refit(up);
    return up;
}

}
