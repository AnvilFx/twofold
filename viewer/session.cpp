#include "session.h"

#include <algorithm>
#include <cmath>

using namespace twofold;

namespace {

constexpr int kSnapshotInterval = 30;
constexpr size_t kMaxRecords = 3600;

void mix(uint64_t& h, const void* data, size_t size) {
    const unsigned char* p = (const unsigned char*)data;
    for (size_t i = 0; i < size; ++i) {
        h ^= p[i];
        h *= 1099511628211ull;
    }
}

uint64_t stateHash(const World& w) {
    uint64_t h = 14695981039346656037ull;
    for (BodyId b : w.bodies()) {
        const twofold::Transform xf = w.transform(b);
        const Vec3 v = w.linearVelocity(b), av = w.angularVelocity(b);
        const bool awake = w.isAwake(b);
        mix(h, &b.index, sizeof(b.index));
        mix(h, &xf, sizeof(xf));
        mix(h, &v, sizeof(v));
        mix(h, &av, sizeof(av));
        mix(h, &awake, sizeof(awake));
    }
    return h;
}

FrameRecord capture(const Scene& s) {
    const World& w = s.world;
    FrameRecord r;
    r.hash = stateHash(w);
    r.profile = w.profile();
    r.bodies.resize(s.bodies.size());
    for (size_t i = 0; i < s.bodies.size(); ++i) {
        const BodyId b = s.bodies[i];
        if (!w.isValid(b)) { continue; }
        const Vec3 v = w.linearVelocity(b), av = w.angularVelocity(b);
        r.bodies[i] = { w.centerOfMass(b), length(v), length(av) };
        if (w.motionType(b) != MotionType::Dynamic) { continue; }
        const MassProperties mp = w.massProperties(b);
        const Mat3 rot = w.rotation(b).toMat3();
        const Mat3 inertia = rot * mp.inertia * transpose(rot);
        r.kineticEnergy += 0.5f * mp.mass * dot(v, v) + 0.5f * dot(av, inertia * av);
    }
    for (const ContactInfo& c : w.contacts()) {
        for (int i = 0; i < c.pointCount; ++i) { r.maxPenetration = std::fmax(r.maxPenetration, -c.separations[i]); }
    }
    return r;
}

}

void Session::load(int sceneIndex, const WorldSettings& settings, float dt) {
    const int n = (int)sceneList().size();
    index_ = (sceneIndex % n + n) % n;
    settings_ = settings;
    dt_ = dt;
    reset();
}

void Session::reset() {
    scene_ = build();
    frame_ = 0;
    base_ = 0;
    divergence_ = -1;
    commands_.clear();
    cursor_ = 0;
    records_.assign(1, capture(scene_));
    snapshots_.clear();
    snapshots_.emplace_back(0, scene_);
}

Scene Session::build() const {
    Scene s;
    s.dt = dt_;
    s.world.setSettings(settings_);
    sceneList()[index_].build(s);
    return s;
}

void Session::advance(Scene& s, size_t& cursor, int& frame) const {
    while (cursor < commands_.size() && commands_[cursor].frame <= frame) { commands_[cursor++].apply(s); }
    if (s.update) { s.update(s, s.dt); }
    s.world.step(s.dt);
    s.time += s.dt;
    ++frame;
}

void Session::truncate() {
    while (!commands_.empty() && commands_.back().frame > frame_) { commands_.pop_back(); }
    cursor_ = std::min(cursor_, commands_.size());
    records_.resize((size_t)(frame_ - base_ + 1));
    while (!snapshots_.empty() && snapshots_.back().first > frame_) { snapshots_.pop_back(); }
    if (divergence_ > frame_) { divergence_ = -1; }
}

void Session::trim() {
    while (records_.size() > kMaxRecords) {
        records_.pop_front();
        ++base_;
    }
    while (snapshots_.size() > 1 && snapshots_[1].first <= base_) { snapshots_.erase(snapshots_.begin()); }
}

void Session::issue(const Command& command) {
    truncate();
    while (cursor_ < commands_.size()) { commands_[cursor_++].apply(scene_); }
    command(scene_);
    commands_.push_back({ frame_, command });
    cursor_ = commands_.size();
}

void Session::step() {
    advance(scene_, cursor_, frame_);
    FrameRecord r = capture(scene_);
    const size_t i = (size_t)(frame_ - base_);
    if (i < records_.size()) {
        if (records_[i].hash != r.hash && divergence_ < 0) { divergence_ = frame_; }
        records_[i] = std::move(r);
    } else {
        records_.push_back(std::move(r));
    }
    if (frame_ % kSnapshotInterval == 0) {
        auto it = std::find_if(snapshots_.begin(), snapshots_.end(), [this](const auto& s) { return s.first >= frame_; });
        if (it != snapshots_.end() && it->first == frame_) {
            it->second = scene_;
        } else {
            snapshots_.emplace(it, frame_, scene_);
        }
    }
    trim();
}

void Session::seek(int target) {
    target = std::clamp(target, firstFrame(), lastFrame());
    if (target == frame_) { return; }
    auto it = std::find_if(snapshots_.rbegin(), snapshots_.rend(), [target](const auto& s) { return s.first <= target; });
    if (it == snapshots_.rend()) { return; }
    if (target < frame_ || it->first > frame_) {
        scene_ = it->second;
        frame_ = it->first;
        cursor_ = (size_t)(std::lower_bound(commands_.begin(), commands_.end(), frame_,
                                            [](const Logged& c, int f) { return c.frame < f; }) - commands_.begin());
    }
    while (frame_ < target) { step(); }
}

int Session::verify() const {
    Scene s = build();
    size_t cursor = 0;
    int frame = 0;
    if (base_ == 0 && stateHash(s.world) != records_.front().hash) { return 0; }
    while (frame < frame_) {
        advance(s, cursor, frame);
        if (frame >= base_ && stateHash(s.world) != record(frame).hash) { return frame; }
    }
    return -1;
}
