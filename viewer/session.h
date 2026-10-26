#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <utility>
#include <vector>

#include "scene.h"

using Command = std::function<void(Scene&)>;

struct BodySample {
    twofold::Vec3 position;
    float speed = 0.0f;
    float angularSpeed = 0.0f;
};

struct FrameRecord {
    uint64_t hash = 0;
    twofold::StepProfile profile;
    float kineticEnergy = 0.0f;
    float maxPenetration = 0.0f;
    // Indexed like Scene::bodies.
    std::vector<BodySample> bodies;
};

// Runs a scene and records every change made to it, so any recorded frame can be restored
// exactly and replays can be checked against the original run.
class Session {
public:
    void load(int sceneIndex, const twofold::WorldSettings& settings, float dt);
    void reset();

    // Applies a change at the current frame and drops the recorded frames after it.
    void issue(const Command& command);
    void step();
    void seek(int frame);
    // Replays from a freshly built scene. Returns the first frame that differs from the record, or -1.
    int verify() const;

    Scene& scene() { return scene_; }
    const Scene& scene() const { return scene_; }
    int sceneIndex() const { return index_; }
    int frame() const { return frame_; }
    int firstFrame() const { return base_; }
    int lastFrame() const { return base_ + (int)records_.size() - 1; }
    const FrameRecord& record(int frame) const { return records_[frame - base_]; }
    // First frame where stepping through recorded frames gave a different state, or -1.
    int divergence() const { return divergence_; }

private:
    struct Logged {
        int frame;
        Command apply;
    };

    Scene build() const;
    void advance(Scene& scene, size_t& cursor, int& frame) const;
    void truncate();
    void trim();

    int index_ = 0;
    twofold::WorldSettings settings_;
    float dt_ = 1.0f / 60.0f;
    Scene scene_;
    int frame_ = 0;
    int base_ = 0;
    int divergence_ = -1;
    std::vector<Logged> commands_;
    size_t cursor_ = 0;
    std::deque<FrameRecord> records_;
    std::vector<std::pair<int, Scene>> snapshots_;
};
