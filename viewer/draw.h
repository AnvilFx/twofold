#pragma once

#include "scene.h"

enum class ColorMode { Material, Motion, Island, Speed };

struct DrawOptions {
    ColorMode colorMode = ColorMode::Material;
    bool wireframe = false;
    bool selectedOnly = false;
    bool contacts = false;
    bool speculative = false;
    bool forces = false;
    bool friction = false;
    float forceScale = 0.05f;
    bool bounds = false;
    bool centers = false;
    bool velocities = false;
    bool planes = false;
    bool joints = true;
};

void drawScene(const Scene& scene, const DrawOptions& options, const twofold::Vec3& focus, twofold::BodyId selected);
void drawLine(const twofold::Vec3& a, const twofold::Vec3& b, Color color);
