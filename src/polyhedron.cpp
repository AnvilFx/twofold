#include <map>
#include <utility>

#include "internal.h"

namespace twofold {

Polyhedron makePolyhedron(std::vector<Vec3> vertices, const std::vector<std::vector<int>>& faces) {
    Polyhedron p;
    p.vertices = std::move(vertices);
    for (const Vec3& v : p.vertices) { p.centroid += v; }
    if (!p.vertices.empty()) {
        p.centroid = p.centroid / (float)p.vertices.size();
        p.lo = p.hi = p.vertices[0];
    }
    for (const Vec3& v : p.vertices) {
        p.lo = minVec(p.lo, v);
        p.hi = maxVec(p.hi, v);
    }

    std::map<std::pair<int, int>, int> directed;
    for (const std::vector<int>& loop : faces) {
        Polyhedron::Face f;
        f.first = (int)p.indices.size();
        f.count = (int)loop.size();
        // Newell's method.
        Vec3 n;
        for (int i = 0; i < f.count; ++i) {
            const Vec3& a = p.vertices[loop[i]];
            const Vec3& b = p.vertices[loop[(i + 1) % f.count]];
            n.x += (a.y - b.y) * (a.z + b.z);
            n.y += (a.z - b.z) * (a.x + b.x);
            n.z += (a.x - b.x) * (a.y + b.y);
        }
        p.normals.push_back(normalize(n));
        p.offsets.push_back(dot(p.normals.back(), p.vertices[loop[0]]));
        const int fi = (int)p.faces.size();
        for (int i = 0; i < f.count; ++i) {
            p.indices.push_back(loop[i]);
            directed[{ loop[i], loop[(i + 1) % f.count] }] = fi;
        }
        p.faces.push_back(f);
    }

    for (const auto& e : directed) {
        const int a = e.first.first, b = e.first.second;
        if (a > b) { continue; }
        auto twin = directed.find({ b, a });
        if (twin == directed.end()) { continue; }
        Polyhedron::Edge edge;
        edge.v0 = a;
        edge.v1 = b;
        edge.f0 = e.second;
        edge.f1 = twin->second;
        edge.length = length(p.vertices[b] - p.vertices[a]);
        p.edges.push_back(edge);
        p.arcs.push_back(cross(p.normals[edge.f1], p.normals[edge.f0]));
    }
    return p;
}

Polyhedron makePrism(const std::vector<Vec2>& polygon, float zMin, float zMax) {
    const int n = (int)polygon.size();
    std::vector<Vec3> verts;
    verts.reserve(2 * n);
    for (const Vec2& v : polygon) { verts.push_back({ v.x, v.y, zMin }); }
    for (const Vec2& v : polygon) { verts.push_back({ v.x, v.y, zMax }); }

    std::vector<std::vector<int>> faces;
    std::vector<int> top, bottom;
    for (int i = 0; i < n; ++i) {
        top.push_back(n + i);
        bottom.push_back(n - 1 - i);
        const int j = (i + 1) % n;
        faces.push_back({ i, j, n + j, n + i });
    }
    faces.push_back(top);
    faces.push_back(bottom);
    return makePolyhedron(std::move(verts), faces);
}

Polyhedron makeBox(const Vec3& h) {
    return makePrism({ { -h.x, -h.y }, { h.x, -h.y }, { h.x, h.y }, { -h.x, h.y } }, -h.z, h.z);
}

MassProperties polyhedronMass(const Polyhedron& poly, float density) {
    MassProperties mp;
    const Vec3 ref = poly.centroid;
    float volume = 0.0f;
    Vec3 center;
    Mat3 cov = Mat3::zero();
    for (const Polyhedron::Face& f : poly.faces) {
        const Vec3 a = poly.vertices[poly.indices[f.first]] - ref;
        for (int i = 1; i + 1 < f.count; ++i) {
            const Vec3 b = poly.vertices[poly.indices[f.first + i]] - ref;
            const Vec3 c = poly.vertices[poly.indices[f.first + i + 1]] - ref;
            const float det = dot(a, cross(b, c));
            volume += det / 6.0f;
            center += (a + b + c) * (det / 24.0f);
            const Vec3 s = a + b + c;
            cov = cov + (Mat3::outer(a, a) + Mat3::outer(b, b) + Mat3::outer(c, c) + Mat3::outer(s, s)) * (det / 120.0f);
        }
    }
    if (volume <= 0.0f) { return mp; }
    center = center / volume;
    mp.mass = density * volume;
    const Mat3 cc = cov * density - Mat3::outer(center, center) * mp.mass;
    const float tr = cc.r[0].x + cc.r[1].y + cc.r[2].z;
    mp.inertia = Mat3::diagonal({ tr, tr, tr }) - cc;
    mp.center = center + ref;
    return mp;
}

}
