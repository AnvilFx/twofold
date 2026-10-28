#pragma once

#include <vector>

#include "math.h"

namespace twofold {

// Faces are vertex loops, counter-clockwise seen from outside.
struct HullMesh {
    std::vector<Vec3> vertices;
    std::vector<std::vector<int>> faces;
};

// Convex hull with at most maxVertices corners, keeping the points that reach farthest out.
// Flat, collinear and single-point clouds are thickened into a thin slab. Non-finite points are skipped, and a
// cloud with none left gives an empty mesh.
HullMesh computeHull(const std::vector<Vec3>& points, int maxVertices = 64);

}
