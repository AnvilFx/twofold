#pragma once

#include <vector>

#include "raylib.h"
#include "twofold/twofold.h"

struct Visual {
    twofold::BodyId body;
    twofold::ShapeId shape;
    twofold::Geometry geometry;
    twofold::Transform local;
    Color color;
    twofold::HullMesh hull;
};

struct Scene {
    twofold::World world;
    std::vector<twofold::BodyId> bodies;
    std::vector<Visual> visuals;
    std::vector<twofold::BodyId> driven;
    float dt = 1.0f / 60.0f;
    float time = 0.0f;
    void (*update)(Scene& scene, float dt) = nullptr;

    twofold::BodyId add(const twofold::BodyDesc& body, const twofold::ShapeDesc& shape, Color color);
    void attach(twofold::BodyId body, const twofold::ShapeDesc& shape, Color color);
};

struct SceneInfo {
    const char* name;
    const char* description;
    void (*build)(Scene& scene);
    twofold::Vec3 target;
    float distance;
    float yaw;
    float pitch;
};

const std::vector<SceneInfo>& sceneList();
