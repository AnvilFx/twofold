#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "draw.h"
#include "imgui.h"
#include "raylib.h"
#include "rlImGui.h"
#include "rlgl.h"
#include "session.h"

using namespace twofold;

namespace {

constexpr float kPi = 3.14159265358979f;
constexpr float kRayLength = 1000.0f;
constexpr int kPlotFrames = 300;
constexpr int kMaxStepsPerFrame = 8;
const int kRates[] = { 30, 60, 120, 240 };
const ImVec4 kWarning{ 1.0f, 0.45f, 0.35f, 1.0f };

Vector3 rl(const Vec3& v) { return { v.x, v.y, v.z }; }
Vec3 tf(const Vector3& v) { return { v.x, v.y, v.z }; }

// Layout lengths are written for a 100% display scale.
float px(float v) { return v * ImGui::GetStyle().FontScaleDpi; }

// TakeScreenshot scales the read by the display scale, which this window does not use.
void saveScreenshot(const char* file) {
    Image image{ rlReadScreenPixels(GetRenderWidth(), GetRenderHeight()), GetRenderWidth(), GetRenderHeight(), 1,
                 PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
    ExportImage(image, file);
    UnloadImage(image);
}

const char* geometryName(const Geometry& g) {
    static const char* names[] = { "Sphere", "Box", "Cylinder", "Prism", "Hull", "Half-space", "Mesh" };
    return names[g.index()];
}

const char* typeName(MotionType t) {
    static const char* names[] = { "Static", "Kinematic", "Dynamic" };
    return names[(int)t];
}

std::string describe(const Geometry& g) {
    char text[128];
    if (const Sphere* s = std::get_if<Sphere>(&g)) {
        std::snprintf(text, sizeof(text), "Sphere  r %.3f", s->radius);
    } else if (const Box* b = std::get_if<Box>(&g)) {
        std::snprintf(text, sizeof(text), "Box  half %.3f %.3f %.3f", b->halfExtents.x, b->halfExtents.y, b->halfExtents.z);
    } else if (const Cylinder* c = std::get_if<Cylinder>(&g)) {
        std::snprintf(text, sizeof(text), "Cylinder  r %.3f  half height %.3f", c->radius, c->halfHeight);
    } else if (const Prism* p = std::get_if<Prism>(&g)) {
        std::snprintf(text, sizeof(text), "Prism  %d sides  z %.3f to %.3f", (int)p->polygon.size(), p->zMin, p->zMax);
    } else if (const Hull* h = std::get_if<Hull>(&g)) {
        std::snprintf(text, sizeof(text), "Hull  %d points  at most %d corners", (int)h->points.size(), h->maxVertices);
    } else if (const twofold::Mesh* m = std::get_if<twofold::Mesh>(&g)) {
        std::snprintf(text, sizeof(text), "Mesh  %d triangles", m->mesh ? m->mesh->triangleCount() : 0);
    } else {
        std::snprintf(text, sizeof(text), "Half-space");
    }
    return text;
}

struct Orbit {
    Vec3 target;
    float distance = 20.0f;
    float yaw = 0.5f;
    float pitch = 0.3f;

    Camera3D camera() const {
        Camera3D cam{};
        const Vec3 dir{ std::cos(pitch) * std::sin(yaw), std::sin(pitch), std::cos(pitch) * std::cos(yaw) };
        cam.position = rl(target + dir * distance);
        cam.target = rl(target);
        cam.up = { 0, 1, 0 };
        cam.fovy = 45.0f;
        cam.projection = CAMERA_PERSPECTIVE;
        return cam;
    }

    void update() {
        if (ImGui::GetIO().WantCaptureMouse) { return; }
        const Vector2 d = GetMouseDelta();
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            yaw -= d.x * 0.006f;
            pitch = std::clamp(pitch + d.y * 0.006f, -1.5f, 1.5f);
        }
        if (IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
            const Vec3 forward = normalize(target - tf(camera().position));
            const Vec3 right = normalize(cross(forward, { 0, 1, 0 }));
            const Vec3 up = cross(right, forward);
            target += (right * -d.x + up * d.y) * (distance * 0.0015f);
        }
        distance = std::clamp(distance * std::pow(0.9f, GetMouseWheelMove()), 2.0f, 200.0f);
    }
};

struct Grab {
    BodyId body;
    Vec3 local;
    Vec3 target;
    float distance = 0.0f;
    bool active = false;
};

struct BodyRow {
    BodyId id;
    MotionType type;
    const char* shape;
    bool awake;
    uint32_t island;
    float speed;
    float height;
};

struct Panel {
    ImVec2 pos;
    ImVec2 size;
    float scale = 0.0f;
};

enum BodyFilter { FilterAll, FilterDynamic, FilterAwake, FilterPlanar };

void plot(const char* label, const std::vector<float>& values, const char* format, float current) {
    float top = 0.0f;
    for (float v : values) { top = std::fmax(top, v); }
    char overlay[64];
    std::snprintf(overlay, sizeof(overlay), format, current);
    ImGui::TextUnformatted(label);
    ImGui::PlotLines((std::string("##") + label).c_str(), values.data(), (int)values.size(), 0, overlay, 0.0f,
                     top > 0.0f ? top * 1.15f : 1.0f, { -1.0f, px(48) });
}

struct App {
    Session session;
    Orbit orbit;
    DrawOptions draw;
    WorldSettings settings;
    int rate = 1;
    bool playing = true;
    float timeScale = 1.0f;
    float accumulator = 0.0f;
    bool labels = false;
    BodyId selected;
    Grab grab;
    float shotRadius = 0.3f;
    float shotSpeed = 30.0f;
    float shotDensity = 6.0f;
    int bodyFilter = FilterAll;
    std::string verifyResult;
    ImGuiStyle baseStyle;
    std::map<std::string, Panel> panels;

    Scene& scene() { return session.scene(); }
    World& world() { return session.scene().world; }

    // Follows the display scale of the monitor the window is on. The window keeps its top edge and
    // horizontal center, within that monitor.
    void updateScale() {
        ImGuiStyle& style = ImGui::GetStyle();
        const float scale = GetWindowScaleDPI().x;
        if (scale == style.FontScaleDpi) { return; }
        const float ratio = scale / style.FontScaleDpi;
        style = baseStyle;
        style.ScaleAllSizes(scale);
        style.FontScaleDpi = scale;
        if (IsWindowMaximized() || IsWindowFullscreen()) { return; }
        const Vector2 pos = GetWindowPosition();
        const int w = GetScreenWidth(), h = GetScreenHeight();
        const int nw = (int)std::lround((float)w * ratio), nh = (int)std::lround((float)h * ratio);
        const int monitor = GetCurrentMonitor();
        const int left = (int)GetMonitorPosition(monitor).x;
        const int right = left + GetMonitorWidth(monitor) - nw;
        SetWindowSize(nw, nh);
        SetWindowPosition(std::max(left, std::min((int)pos.x + (w - nw) / 2, right)), (int)pos.y);
    }

    // A panel laid out at another display scale is scaled with it, so it keeps its place and proportions.
    // A zero height fits the panel to its content.
    bool beginPanel(const char* name, const ImVec2& pos, const ImVec2& size) {
        Panel& p = panels[name];
        const float scale = px(1.0f);
        if (p.scale > 0.0f && p.scale != scale) {
            const float r = scale / p.scale;
            ImGui::SetNextWindowPos({ p.pos.x * r, p.pos.y * r });
            ImGui::SetNextWindowSize({ p.size.x * r, size.y > 0.0f ? p.size.y * r : 0.0f });
        } else {
            ImGui::SetNextWindowPos(pos, ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(size, ImGuiCond_FirstUseEver);
        }
        const bool open = ImGui::Begin(name);
        p.pos = ImGui::GetWindowPos();
        if (open || p.scale != scale) { p.size = ImGui::GetWindowSize(); }
        p.scale = scale;
        return open;
    }

    void load(int index, bool resetCamera) {
        session.load(index, settings, 1.0f / (float)kRates[rate]);
        grab.active = false;
        selected = {};
        verifyResult.clear();
        if (resetCamera) {
            const SceneInfo& info = sceneList()[session.sceneIndex()];
            orbit.target = info.target;
            orbit.distance = info.distance;
            orbit.yaw = info.yaw;
            orbit.pitch = info.pitch;
        }
    }

    void reset() {
        session.reset();
        grab.active = false;
        verifyResult.clear();
    }

    void edit(const Command& command) { session.issue(command); }

    void seek(int frame) {
        playing = false;
        grab.active = false;
        session.seek(frame);
    }

    void applyGrab() {
        World& w = world();
        const BodyId b = grab.body;
        if (!w.isValid(b) || w.motionType(b) != MotionType::Dynamic) {
            grab.active = false;
            return;
        }
        const MassProperties mp = w.massProperties(b);
        const Vec3 point = w.transform(b).apply(grab.local);
        const Vec3 av = w.angularVelocity(b);
        const Vec3 v = w.linearVelocity(b) + cross(av, point - w.centerOfMass(b));
        const float omega = 2.0f * kPi * 4.0f;
        const Vec3 g = w.settings().gravity * w.gravityScale(b);
        Vec3 force = ((grab.target - point) * (omega * omega) - v * (1.6f * omega) - g) * mp.mass;
        const float limit = 400.0f * mp.mass * w.settings().lengthUnit;
        if (length(force) > limit) { force = normalize(force) * limit; }
        const Mat3 rot = w.rotation(b).toMat3();
        const Vec3 torque = (rot * mp.inertia * transpose(rot)) * av * -4.0f;
        edit([b, force, point, torque](Scene& s) {
            s.world.addForceAtPoint(b, force, point);
            s.world.addTorque(b, torque);
        });
    }

    void stepForward() {
        if (grab.active) { applyGrab(); }
        session.step();
    }

    void advance() {
        if (!playing) {
            accumulator = 0.0f;
            return;
        }
        const float dt = scene().dt;
        accumulator += std::fmin(GetFrameTime(), 0.1f) * timeScale;
        int steps = 0;
        while (accumulator >= dt && steps < kMaxStepsPerFrame) {
            stepForward();
            accumulator -= dt;
            ++steps;
        }
        if (steps == kMaxStepsPerFrame) { accumulator = 0.0f; }
    }

    void pick(const Camera3D& cam) {
        const Ray ray = GetScreenToWorldRay(GetMousePosition(), cam);
        const Vec3 dir = normalize(tf(ray.direction));
        RayHit hit;
        const World& w = world();
        if (!w.castRay(tf(ray.position), dir * kRayLength, hit) || std::holds_alternative<HalfSpace>(w.geometry(hit.shape))) {
            selected = {};
            return;
        }
        selected = hit.body;
        if (w.motionType(hit.body) == MotionType::Dynamic) {
            grab.body = hit.body;
            grab.local = w.transform(hit.body).applyInverse(hit.point);
            grab.distance = hit.fraction * kRayLength;
            grab.target = hit.point;
            grab.active = true;
        }
    }

    void shoot(const Camera3D& cam) {
        const Ray ray = GetScreenToWorldRay(GetMousePosition(), cam);
        const Vec3 dir = normalize(tf(ray.direction));
        BodyDesc bd;
        bd.position = tf(ray.position) + dir * 1.0f;
        bd.linearVelocity = dir * shotSpeed;
        ShapeDesc sd;
        sd.geometry = Sphere{ shotRadius };
        sd.density = shotDensity;
        edit([bd, sd](Scene& s) { s.add(bd, sd, { 235, 235, 240, 255 }); });
    }

    void focusSelected() {
        if (world().isValid(selected)) { orbit.target = world().centerOfMass(selected); }
    }

    void deleteSelected() {
        if (!world().isValid(selected)) { return; }
        const BodyId b = selected;
        edit([b](Scene& s) { s.world.destroyBody(b); });
        selected = {};
        grab.active = false;
    }

    void input(const Camera3D& cam) {
        const ImGuiIO& io = ImGui::GetIO();
        if (!io.WantCaptureKeyboard) {
            for (int k = 0; k < (int)sceneList().size() && k < 10; ++k) {
                if (IsKeyPressed(k < 9 ? KEY_ONE + k : KEY_ZERO)) { load(k, true); }
            }
            if (IsKeyPressed(KEY_R)) { reset(); }
            if (IsKeyPressed(KEY_SPACE)) { playing = !playing; }
            if (IsKeyPressed(KEY_RIGHT) || IsKeyPressedRepeat(KEY_RIGHT)) {
                playing = false;
                stepForward();
            }
            if (IsKeyPressed(KEY_LEFT) || IsKeyPressedRepeat(KEY_LEFT)) { seek(session.frame() - 1); }
            if (IsKeyPressed(KEY_HOME)) { seek(session.firstFrame()); }
            if (IsKeyPressed(KEY_END)) { seek(session.lastFrame()); }
            if (IsKeyPressed(KEY_B)) { shoot(cam); }
            if (IsKeyPressed(KEY_F)) { focusSelected(); }
            if (IsKeyPressed(KEY_C)) { draw.contacts = !draw.contacts; }
            if (IsKeyPressed(KEY_DELETE)) { deleteSelected(); }
            if (IsKeyPressed(KEY_ESCAPE)) { selected = {}; }
        }
        if (!io.WantCaptureMouse && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) { pick(cam); }
        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) { grab.active = false; }
        if (grab.active) {
            const Ray ray = GetScreenToWorldRay(GetMousePosition(), cam);
            grab.target = tf(ray.position) + normalize(tf(ray.direction)) * grab.distance;
        }
    }

    std::vector<float> history(float (*value)(const FrameRecord&)) const {
        std::vector<float> out;
        const int last = session.frame();
        for (int f = std::max(session.firstFrame(), last - kPlotFrames + 1); f <= last; ++f) {
            out.push_back(value(session.record(f)));
        }
        return out;
    }

    void timeline() {
        const SceneInfo& info = sceneList()[session.sceneIndex()];
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::BeginCombo("##scene", info.name)) {
            for (int i = 0; i < (int)sceneList().size(); ++i) {
                if (ImGui::Selectable(sceneList()[i].name, i == session.sceneIndex())) { load(i, true); }
            }
            ImGui::EndCombo();
        }
        ImGui::TextWrapped("%s", info.description);
        ImGui::Spacing();

        if (ImGui::Button("|<")) { seek(session.firstFrame()); }
        ImGui::SameLine();
        if (ImGui::Button("<")) { seek(session.frame() - 1); }
        ImGui::SameLine();
        if (ImGui::Button(playing ? "Pause" : "Play", { px(60), 0 })) { playing = !playing; }
        ImGui::SameLine();
        if (ImGui::Button(">")) {
            playing = false;
            stepForward();
        }
        ImGui::SameLine();
        if (ImGui::Button(">|")) { seek(session.lastFrame()); }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) { reset(); }

        int frame = session.frame();
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::SliderInt("##frame", &frame, session.firstFrame(), std::max(session.lastFrame(), session.firstFrame() + 1),
                             "frame %d")) {
            seek(frame);
        }
        ImGui::Text("t %.3f s   recorded %d to %d", scene().time, session.firstFrame(), session.lastFrame());
        ImGui::SetNextItemWidth(px(150));
        ImGui::SliderFloat("Speed", &timeScale, 0.02f, 4.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
        ImGui::SetNextItemWidth(px(150));
        if (ImGui::Combo("Step rate", &rate, "30 Hz\0" "60 Hz\0" "120 Hz\0" "240 Hz\0")) {
            const float dt = 1.0f / (float)kRates[rate];
            edit([dt](Scene& s) { s.dt = dt; });
        }

        if (session.divergence() >= 0) {
            ImGui::TextColored(kWarning, "Replay differs from the record at frame %d", session.divergence());
        }
        if (ImGui::Button("Verify determinism")) {
            const int bad = session.verify();
            char text[96];
            if (bad < 0) {
                std::snprintf(text, sizeof(text), "A fresh replay matches all %d frames", session.frame());
            } else {
                std::snprintf(text, sizeof(text), "A fresh replay differs at frame %d", bad);
            }
            verifyResult = text;
        }
        if (!verifyResult.empty()) {
            ImGui::SameLine();
            ImGui::TextUnformatted(verifyResult.c_str());
        }
    }

    void drawOptions() {
        int mode = (int)draw.colorMode;
        ImGui::SetNextItemWidth(px(150));
        if (ImGui::Combo("Color", &mode, "Material\0" "Motion state\0" "Island\0" "Speed\0")) { draw.colorMode = (ColorMode)mode; }
        ImGui::Checkbox("Wireframe", &draw.wireframe);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Labels", &labels);
        ImGui::Checkbox("Contacts", &draw.contacts);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Speculative", &draw.speculative);
        ImGui::Checkbox("Normal forces", &draw.forces);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Friction forces", &draw.friction);
        if (draw.forces || draw.friction) {
            ImGui::SetNextItemWidth(px(150));
            ImGui::SliderFloat("Length per newton", &draw.forceScale, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic);
        }
        ImGui::Checkbox("Bounds", &draw.bounds);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Centers and axes", &draw.centers);
        ImGui::Checkbox("Velocities", &draw.velocities);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Lock planes", &draw.planes);
        ImGui::Checkbox("Joints", &draw.joints);
        ImGui::SameLine(px(160));
        ImGui::Checkbox("Only the selected body", &draw.selectedOnly);
    }

    void worldSettings() {
        WorldSettings s = world().settings();
        bool changed = false;
        changed |= ImGui::SliderFloat("Gravity", &s.gravity.y, -30.0f, 0.0f, "%.2f");
        changed |= ImGui::SliderInt("Substeps", &s.substeps, 1, 16);
        changed |= ImGui::SliderFloat("Contact hertz", &s.contactHertz, 5.0f, 120.0f, "%.0f");
        changed |= ImGui::SliderFloat("Damping ratio", &s.contactDampingRatio, 0.0f, 20.0f, "%.1f");
        changed |= ImGui::SliderFloat("Max pushout", &s.maxPushoutSpeed, 0.1f, 20.0f, "%.1f");
        changed |= ImGui::SliderFloat("Bounce threshold", &s.restitutionThreshold, 0.0f, 5.0f, "%.2f");
        changed |= ImGui::SliderFloat("Sleep threshold", &s.sleepThreshold, 0.0f, 0.5f, "%.3f");
        changed |= ImGui::SliderFloat("Time to sleep", &s.timeToSleep, 0.0f, 3.0f, "%.2f s");
        changed |= ImGui::Checkbox("Continuous collision", &s.continuous);
        if (ImGui::Button("Defaults")) {
            s = WorldSettings{};
            changed = true;
        }
        if (changed) {
            settings = s;
            edit([s](Scene& sc) { sc.world.setSettings(s); });
        }
    }

    void controlPanel() {
        if (beginPanel("Twofold", { px(10), px(10) }, { px(360), 0 })) {
            timeline();
            if (ImGui::CollapsingHeader("Draw", ImGuiTreeNodeFlags_DefaultOpen)) { drawOptions(); }
            if (ImGui::CollapsingHeader("World")) { worldSettings(); }
            if (ImGui::CollapsingHeader("Shooting")) {
                ImGui::SliderFloat("Radius", &shotRadius, 0.05f, 2.0f, "%.2f");
                ImGui::SliderFloat("Speed##shot", &shotSpeed, 0.0f, 100.0f, "%.0f");
                ImGui::SliderFloat("Density", &shotDensity, 0.1f, 50.0f, "%.1f", ImGuiSliderFlags_Logarithmic);
            }
            if (ImGui::CollapsingHeader("Controls")) {
                ImGui::TextUnformatted(
                    "Left click        select, drag to pull\n"
                    "B                 shoot a ball\n"
                    "Space             play or pause\n"
                    "Right, Left       step forward, back\n"
                    "Home, End         first, last frame\n"
                    "R                 reset the scene\n"
                    "1 to 9, 0         change scene\n"
                    "F                 focus the selection\n"
                    "C                 contacts\n"
                    "Delete            delete the selection\n"
                    "Escape            clear the selection\n"
                    "Right drag        orbit\n"
                    "Middle drag       pan\n"
                    "Wheel             zoom");
            }
        }
        ImGui::End();
    }

    void statsPanel() {
        if (!beginPanel("Step", { (float)GetScreenWidth() - px(330), px(10) }, { px(320), 0 })) {
            ImGui::End();
            return;
        }
        const FrameRecord& r = session.record(session.frame());
        const StepProfile& p = r.profile;
        if (ImGui::BeginTable("phases", 2, ImGuiTableFlags_SizingStretchProp)) {
            const char* names[] = { "Broadphase", "Narrowphase", "Solver", "Continuous", "Sleep", "Total" };
            const float values[] = { p.broadphase, p.narrowphase, p.solver, p.continuous, p.sleep, p.total };
            for (int i = 0; i < 6; ++i) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(names[i]);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f ms", values[i]);
            }
            ImGui::EndTable();
        }
        ImGui::Text("Bodies %d   awake %d   islands %d", (int)world().bodies().size(), p.awakeBodies, p.awakeIslands);
        ImGui::Text("Pairs %d   manifolds %d   touching %d", p.candidatePairs, p.manifolds, p.touchingContacts);
        ImGui::Text("Swept pairs %d   impacts %d", p.sweptPairs, p.impacts);
        ImGui::Separator();
        const float mm = 1000.0f / world().settings().lengthUnit;
        plot("Step time", history([](const FrameRecord& f) { return f.profile.total; }), "%.3f ms", p.total);
        plot("Kinetic energy", history([](const FrameRecord& f) { return f.kineticEnergy; }), "%.2f J", r.kineticEnergy);
        std::vector<float> depth = history([](const FrameRecord& f) { return f.maxPenetration; });
        for (float& d : depth) { d *= mm; }
        plot("Deepest penetration", depth, "%.2f mm", r.maxPenetration * mm);
        plot("Awake bodies", history([](const FrameRecord& f) { return (float)f.profile.awakeBodies; }), "%.0f",
             (float)p.awakeBodies);
        ImGui::End();
    }

    std::vector<BodyRow> bodyRows() {
        const World& w = world();
        std::vector<BodyRow> rows;
        for (BodyId b : scene().bodies) {
            if (!w.isValid(b)) { continue; }
            const MotionType type = w.motionType(b);
            const bool awake = type == MotionType::Dynamic ? w.isAwake(b) : false;
            if (bodyFilter == FilterDynamic && type != MotionType::Dynamic) { continue; }
            if (bodyFilter == FilterAwake && !awake) { continue; }
            if (bodyFilter == FilterPlanar && (w.locks(b) & LockPlanar) != LockPlanar) { continue; }
            const std::vector<ShapeId> shapes = w.shapes(b);
            const char* shape = shapes.empty() ? "None" : shapes.size() > 1 ? "Compound" : geometryName(w.geometry(shapes[0]));
            rows.push_back({ b, type, shape, awake, w.island(b), length(w.linearVelocity(b)), w.centerOfMass(b).y });
        }
        return rows;
    }

    void bodiesPanel() {
        if (!beginPanel("Bodies", { px(10), (float)GetScreenHeight() - px(330) }, { px(460), px(320) })) {
            ImGui::End();
            return;
        }
        ImGui::SetNextItemWidth(px(150));
        ImGui::Combo("Show", &bodyFilter, "All\0" "Dynamic\0" "Awake\0" "Planar\0");
        std::vector<BodyRow> rows = bodyRows();
        const ImGuiTableFlags flags = ImGuiTableFlags_Sortable | ImGuiTableFlags_ScrollY | ImGuiTableFlags_RowBg |
                                      ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_Resizable;
        if (ImGui::BeginTable("bodies", 7, flags)) {
            ImGui::TableSetupScrollFreeze(0, 1);
            ImGui::TableSetupColumn("Body", ImGuiTableColumnFlags_DefaultSort);
            ImGui::TableSetupColumn("Type");
            ImGui::TableSetupColumn("Shape");
            ImGui::TableSetupColumn("State");
            ImGui::TableSetupColumn("Island");
            ImGui::TableSetupColumn("Speed");
            ImGui::TableSetupColumn("Height");
            ImGui::TableHeadersRow();

            if (const ImGuiTableSortSpecs* specs = ImGui::TableGetSortSpecs(); specs && specs->SpecsCount > 0) {
                const int column = specs->Specs[0].ColumnIndex;
                const bool ascending = specs->Specs[0].SortDirection == ImGuiSortDirection_Ascending;
                auto key = [column](const BodyRow& r) -> double {
                    switch (column) {
                        case 1: return (double)r.type;
                        case 2: return (double)(unsigned char)r.shape[0];
                        case 3: return r.awake ? 1.0 : 0.0;
                        case 4: return (double)r.island;
                        case 5: return (double)r.speed;
                        case 6: return (double)r.height;
                        default: return (double)r.id.index;
                    }
                };
                std::stable_sort(rows.begin(), rows.end(), [&](const BodyRow& a, const BodyRow& b) {
                    return ascending ? key(a) < key(b) : key(a) > key(b);
                });
            }

            ImGuiListClipper clipper;
            clipper.Begin((int)rows.size());
            while (clipper.Step()) {
                for (int i = clipper.DisplayStart; i < clipper.DisplayEnd; ++i) {
                    const BodyRow& r = rows[i];
                    ImGui::PushID((int)r.id.index);
                    ImGui::TableNextRow();
                    ImGui::TableNextColumn();
                    char label[32];
                    std::snprintf(label, sizeof(label), "%u", r.id.index);
                    if (ImGui::Selectable(label, r.id == selected, ImGuiSelectableFlags_SpanAllColumns)) { selected = r.id; }
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(typeName(r.type));
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(r.shape);
                    ImGui::TableNextColumn();
                    ImGui::TextUnformatted(r.type != MotionType::Dynamic ? "" : r.awake ? "awake" : "asleep");
                    ImGui::TableNextColumn();
                    if (r.island != UINT32_MAX) { ImGui::Text("%u", r.island); }
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", r.speed);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", r.height);
                    ImGui::PopID();
                }
            }
            ImGui::EndTable();
        }
        ImGui::End();
    }

    void bodyHistory(BodyId b) {
        const std::vector<BodyId>& bodies = scene().bodies;
        const auto it = std::find(bodies.begin(), bodies.end(), b);
        if (it == bodies.end()) { return; }
        const size_t slot = (size_t)(it - bodies.begin());
        std::vector<float> height, speed, spin;
        const int last = session.frame();
        for (int f = std::max(session.firstFrame(), last - kPlotFrames + 1); f <= last; ++f) {
            const FrameRecord& r = session.record(f);
            const BodySample s = slot < r.bodies.size() ? r.bodies[slot] : BodySample{};
            height.push_back(s.position.y);
            speed.push_back(s.speed);
            spin.push_back(s.angularSpeed);
        }
        float low = 0.0f;
        for (float h : height) { low = std::fmin(low, h); }
        for (float& h : height) { h -= low; }
        const World& w = world();
        plot("Height", height, "%.3f", w.centerOfMass(b).y);
        plot("Speed", speed, "%.3f", length(w.linearVelocity(b)));
        plot("Angular speed", spin, "%.3f rad/s", length(w.angularVelocity(b)));
    }

    void contactTable(BodyId b) {
        const World& w = world();
        const float toForce = 1.0f / scene().dt;
        const ImGuiTableFlags flags = ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_SizingStretchProp;
        if (!ImGui::BeginTable("contacts", 5, flags)) { return; }
        ImGui::TableSetupColumn("Other");
        ImGui::TableSetupColumn("Points");
        ImGui::TableSetupColumn("Min sep (mm)");
        ImGui::TableSetupColumn("Force (N)");
        ImGui::TableSetupColumn("Friction");
        ImGui::TableHeadersRow();
        BodyId next = b;
        for (const ContactInfo& c : w.contacts(true)) {
            if (c.bodyA != b && c.bodyB != b) { continue; }
            const BodyId other = c.bodyA == b ? c.bodyB : c.bodyA;
            float minSep = 1e9f, force = 0.0f;
            for (int i = 0; i < c.pointCount; ++i) {
                minSep = std::fmin(minSep, c.separations[i]);
                force += c.normalImpulses[i] * toForce;
            }
            ImGui::PushID((int)other.index);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            char label[32];
            std::snprintf(label, sizeof(label), "%u", other.index);
            if (ImGui::SmallButton(label)) { next = other; }
            ImGui::TableNextColumn();
            ImGui::Text("%d", c.pointCount);
            ImGui::TableNextColumn();
            if (minSep > 0.0f) {
                ImGui::TextDisabled("%.2f", minSep * 1000.0f);
            } else {
                ImGui::Text("%.2f", minSep * 1000.0f);
            }
            ImGui::TableNextColumn();
            ImGui::Text("%.2f", force);
            ImGui::TableNextColumn();
            ImGui::Text("%.2f", c.friction);
            ImGui::PopID();
        }
        ImGui::EndTable();
        selected = next;
    }

    void inspector() {
        if (selected == BodyId{}) { return; }
        if (!beginPanel("Inspector", { (float)GetScreenWidth() - px(380), px(420) }, { px(370), px(480) })) {
            ImGui::End();
            return;
        }
        World& w = world();
        const BodyId b = selected;
        if (!w.isValid(b)) {
            ImGui::Text("Body %u does not exist at this frame.", b.index);
            ImGui::End();
            return;
        }

        ImGui::Text("Body %u", b.index);
        ImGui::SameLine();
        ImGui::TextDisabled("generation %u", b.generation);
        int type = (int)w.motionType(b);
        ImGui::SetNextItemWidth(px(120));
        if (ImGui::Combo("Type", &type, "Static\0" "Kinematic\0" "Dynamic\0")) {
            edit([b, type](Scene& s) { s.world.setMotionType(b, (MotionType)type); });
        }
        if (w.motionType(b) == MotionType::Dynamic) {
            ImGui::Text("%s   island %u", w.isAwake(b) ? "Awake" : "Asleep", w.island(b));
            if (ImGui::Button("Wake")) { edit([b](Scene& s) { s.world.wake(b); }); }
            ImGui::SameLine();
            if (ImGui::Button("Sleep")) { edit([b](Scene& s) { s.world.sleep(b); }); }
            ImGui::SameLine();
            bool canSleep = w.canSleep(b);
            if (ImGui::Checkbox("Can sleep", &canSleep)) { edit([b, canSleep](Scene& s) { s.world.setCanSleep(b, canSleep); }); }
        }
        if (ImGui::Button("Focus")) { focusSelected(); }
        ImGui::SameLine();
        if (ImGui::Button("Delete")) {
            deleteSelected();
            ImGui::End();
            return;
        }

        ImGui::SeparatorText("Motion");
        Vec3 p = w.position(b);
        if (ImGui::DragFloat3("Position", &p.x, 0.01f)) {
            const Quat q = w.rotation(b);
            edit([b, p, q](Scene& s) { s.world.setTransform(b, p, q); });
        }
        const Quat q = w.rotation(b);
        const float angle = 2.0f * std::acos(std::clamp(q.w, -1.0f, 1.0f));
        const float sinHalf = std::sqrt(std::fmax(0.0f, 1.0f - q.w * q.w));
        const Vec3 axis = sinHalf > 1e-6f ? Vec3{ q.x, q.y, q.z } / sinHalf : Vec3{ 1, 0, 0 };
        ImGui::Text("Rotation %.1f deg about %.3f %.3f %.3f", angle * 180.0f / kPi, axis.x, axis.y, axis.z);
        Vec3 v = w.linearVelocity(b);
        if (ImGui::DragFloat3("Velocity", &v.x, 0.01f)) { edit([b, v](Scene& s) { s.world.setLinearVelocity(b, v); }); }
        Vec3 av = w.angularVelocity(b);
        if (ImGui::DragFloat3("Angular", &av.x, 0.01f)) { edit([b, av](Scene& s) { s.world.setAngularVelocity(b, av); }); }

        ImGui::SeparatorText("Mass");
        const MassProperties mp = w.massProperties(b);
        ImGui::Text("Mass %.4f   center %.3f %.3f %.3f", mp.mass, mp.center.x, mp.center.y, mp.center.z);
        ImGui::Text("Inertia %.4f %.4f %.4f", mp.inertia.r[0].x, mp.inertia.r[1].y, mp.inertia.r[2].z);
        float gravityScale = w.gravityScale(b);
        if (ImGui::DragFloat("Gravity scale", &gravityScale, 0.01f)) {
            edit([b, gravityScale](Scene& s) { s.world.setGravityScale(b, gravityScale); });
        }
        float damping[2] = { w.linearDamping(b), w.angularDamping(b) };
        if (ImGui::DragFloat2("Damping", damping, 0.01f, 0.0f, 100.0f)) {
            const float linear = damping[0], angular = damping[1];
            edit([b, linear, angular](Scene& s) { s.world.setDamping(b, linear, angular); });
        }

        ImGui::SeparatorText("Locks");
        const uint8_t locks = w.locks(b);
        static const char* names[] = { "Move X", "Move Y", "Move Z", "Turn X", "Turn Y", "Turn Z" };
        for (int i = 0; i < 6; ++i) {
            bool on = (locks & (1 << i)) != 0;
            if (i % 3) { ImGui::SameLine(); }
            if (ImGui::Checkbox(names[i], &on)) {
                const uint8_t next = (uint8_t)(on ? (locks | (1 << i)) : (locks & ~(1 << i)));
                const Quat frame = w.lockFrame(b);
                edit([b, next, frame](Scene& s) { s.world.setLocks(b, next, frame); });
            }
        }

        ImGui::SeparatorText("Shapes");
        for (ShapeId s : w.shapes(b)) {
            ImGui::PushID((int)s.index);
            ImGui::TextUnformatted(describe(w.geometry(s)).c_str());
            float friction = w.friction(s), restitution = w.restitution(s), rolling = w.rollingResistance(s);
            if (ImGui::SliderFloat("Friction", &friction, 0.0f, 1.5f, "%.2f")) {
                edit([s, friction](Scene& sc) { sc.world.setFriction(s, friction); });
            }
            if (ImGui::SliderFloat("Restitution", &restitution, 0.0f, 1.0f, "%.2f")) {
                edit([s, restitution](Scene& sc) { sc.world.setRestitution(s, restitution); });
            }
            if (ImGui::SliderFloat("Rolling resistance", &rolling, 0.0f, 0.5f, "%.3f")) {
                edit([s, rolling](Scene& sc) { sc.world.setRollingResistance(s, rolling); });
            }
            ImGui::PopID();
        }

        ImGui::SeparatorText("Contacts");
        contactTable(b);
        ImGui::SeparatorText("History");
        bodyHistory(b);
        ImGui::End();
    }

    void bodyLabels(const Camera3D& cam) {
        ImDrawList* list = ImGui::GetBackgroundDrawList();
        const World& w = world();
        for (BodyId b : scene().bodies) {
            if (!w.isValid(b) || w.motionType(b) == MotionType::Static) { continue; }
            const Vector2 s = GetWorldToScreen(rl(w.centerOfMass(b)), cam);
            char text[16];
            std::snprintf(text, sizeof(text), "%u", b.index);
            const ImU32 color = b == selected ? IM_COL32(255, 220, 60, 255) : IM_COL32(255, 255, 255, 220);
            list->AddText({ s.x - px(6), s.y - px(7) }, color, text);
        }
    }

    void render(const Camera3D& cam, const char* screenshot) {
        BeginDrawing();
        ClearBackground({ 28, 30, 36, 255 });
        BeginMode3D(cam);
        drawScene(scene(), draw, orbit.target, selected);
        if (grab.active && world().isValid(grab.body)) {
            const Vec3 anchor = world().transform(grab.body).apply(grab.local);
            drawLine(anchor, grab.target, WHITE);
            DrawSphereEx(rl(grab.target), 0.05f, 4, 6, WHITE);
        }
        EndMode3D();
        rlImGuiBegin();
        if (labels) { bodyLabels(cam); }
        controlPanel();
        statsPanel();
        bodiesPanel();
        inspector();
        rlImGuiEnd();
        if (screenshot) { saveScreenshot(screenshot); }
        EndDrawing();
    }
};

}

// Usage: twofold_viewer [--scene N] [--frames N --screenshot file.png] [--contacts]
// With --screenshot the viewer runs the scene for the given frames, saves one image to the
// working directory and exits.
int main(int argc, char** argv) {
    int sceneIndex = 0, frames = 0;
    bool contacts = false;
    const char* screenshot = nullptr;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const bool hasValue = i + 1 < argc;
        if (arg == "--scene" && hasValue) {
            sceneIndex = std::atoi(argv[++i]) - 1;
        } else if (arg == "--frames" && hasValue) {
            frames = std::atoi(argv[++i]);
        } else if (arg == "--screenshot" && hasValue) {
            screenshot = argv[++i];
        } else if (arg == "--contacts") {
            contacts = true;
        }
    }

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE | FLAG_VSYNC_HINT);
    InitWindow(1440, 860, "Twofold viewer");
    SetExitKey(KEY_NULL);
    SetTargetFPS(60);
    rlImGuiSetup(true);
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::GetStyle().FontSizeBase = 15.0f;

    App app;
    app.baseStyle = ImGui::GetStyle();
    app.updateScale();
    const int monitor = GetCurrentMonitor();
    const Vector2 corner = GetMonitorPosition(monitor);
    SetWindowPosition((int)corner.x + (GetMonitorWidth(monitor) - GetScreenWidth()) / 2,
                      (int)corner.y + (GetMonitorHeight(monitor) - GetScreenHeight()) / 2);
    app.load(sceneIndex, true);
    app.draw.contacts = contacts;
    int frame = 0;
    while (!WindowShouldClose()) {
        app.updateScale();
        app.orbit.update();
        const Camera3D cam = app.orbit.camera();
        if (screenshot) {
            app.stepForward();
            const bool last = ++frame >= frames;
            app.render(cam, last ? screenshot : nullptr);
            if (last) { break; }
        } else {
            app.input(cam);
            app.advance();
            app.render(cam, nullptr);
        }
    }
    rlImGuiShutdown();
    CloseWindow();
    return 0;
}
