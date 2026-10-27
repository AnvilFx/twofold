#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "math.h"

namespace twofold {

// Triangles indexed three to a face, with the tree and edge data collision needs. Build once and share it
// between shapes. Vertices at the same position are merged and degenerate triangles are dropped.
class TriangleMesh {
public:
    TriangleMesh(const std::vector<Vec3>& vertices, const std::vector<uint32_t>& indices);
    ~TriangleMesh();
    TriangleMesh(const TriangleMesh&) = delete;
    TriangleMesh& operator=(const TriangleMesh&) = delete;

    const std::vector<Vec3>& vertices() const;
    // Three per triangle, in the mesh's own order.
    const std::vector<uint32_t>& indices() const;
    int triangleCount() const;

    struct Data;
    const Data& data() const { return *data_; }

private:
    std::unique_ptr<Data> data_;
};

}
