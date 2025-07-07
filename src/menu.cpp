#include "menu.h"
#include "globals.h"
#include "renderer.h"
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <mutex>

Menu::Menu() : window_(nullptr) {
}

Menu::~Menu() {
}

bool Menu::init(SDL_Window* window, Simulation& sim) {
    window_ = window;
    // ImGui is initialized in renderer.cpp, so no need to reinitialize here
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.5f; // Increase font size by 1.5x for better visibility
    return true;
}

void Menu::render(Simulation& sim) {
    std::lock_guard<std::mutex> lock(global_state.sim_mutex);

    ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::SetWindowPos(ImVec2(10, 200)); // Match original menu position

    // Row 1: Display modes
    ImGui::Text("Display Modes");
    if (ImGui::Button("1##Points")) {
        global_state.current_display_mode = global_state.DisplayMode::POINTS;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Points mode");
    ImGui::SameLine();
    if (ImGui::Button("2##Isosurface")) {
        global_state.current_display_mode = global_state.DisplayMode::ISOSURFACE;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Isosurface mode");
    ImGui::SameLine();
    if (ImGui::Button("3##Wireframe")) {
        global_state.current_display_mode = global_state.DisplayMode::WIREFRAME;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Wireframe mode");
    ImGui::SameLine();
    if (ImGui::Button("4##Particles")) {
        global_state.current_display_mode = global_state.DisplayMode::PARTICLES;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Particles mode");
    ImGui::SameLine();
    if (ImGui::Button("5##Hybrid")) {
        global_state.current_display_mode = global_state.DisplayMode::HYBRID;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Hybrid mode");
    ImGui::SameLine();
    if (ImGui::Button("6##Surface")) {
        global_state.current_display_mode = global_state.DisplayMode::SURFACE;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Surface mode");
    ImGui::SameLine();
    if (ImGui::Button("7##Sphere")) {
        global_state.current_display_mode = global_state.DisplayMode::SPHERE_POINTS;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Sphere mode");
    ImGui::SameLine();
    if (ImGui::Button("0##Center")) {
        global_state.camera_angle = 0.0f;
        global_state.camera_tilt = 0.0f;
        global_state.camera_zoom = 0.3f;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Center camera");

    // Row 2: Camera navigation
    ImGui::Text("Camera Navigation");
    if (ImGui::Button("Up##ZoomIn")) {
        global_state.camera_zoom *= 0.95f;
        global_state.camera_zoom = std::max(global_state.camera_zoom, 0.001f);
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Zoom in");
    ImGui::SameLine();
    if (ImGui::Button("Left##PrevEq")) {
        sim.prev_equation();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Previous equation");
    ImGui::SameLine();
    if (ImGui::Button("Right##NextEq")) {
        sim.next_equation();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Next equation");
    ImGui::SameLine();
    if (ImGui::Button("Down##ZoomOut")) {
        global_state.camera_zoom *= 1.05f;
        global_state.camera_zoom = std::min(global_state.camera_zoom, 5.0f);
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Zoom out");

    // Row 3: Simulation controls
    ImGui::Text("Simulation Controls");
    if (ImGui::Button("Space##Pause")) {
        global_state.is_paused = !global_state.is_paused;
        global_state.rotate_camera = !global_state.is_paused;
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Toggle pause");
    ImGui::SameLine();
    if (ImGui::Button("Enter##Step")) {
        if (global_state.is_paused) {
            sim.step();
            global_state.needs_render = true;
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Step simulation (if paused)");

    // Row 4: Camera pan/tilt (if paused)
    ImGui::Text("Camera Pan/Tilt (Paused)");
    if (ImGui::Button("C##PanLeft")) {
        if (global_state.is_paused) {
            global_state.camera_angle -= 0.01f;
            global_state.needs_render = true;
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Pan camera left (if paused)");
    ImGui::SameLine();
    if (ImGui::Button("V##PanRight")) {
        if (global_state.is_paused) {
            global_state.camera_angle += 0.01f;
            global_state.needs_render = true;
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Pan camera right (if paused)");
    ImGui::SameLine();
    if (ImGui::Button("B##TiltUp")) {
        if (global_state.is_paused) {
            global_state.camera_tilt += 0.01f;
            global_state.needs_render = true;
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Tilt camera up (if paused)");
    ImGui::SameLine();
    if (ImGui::Button("G##TiltDown")) {
        if (global_state.is_paused) {
            global_state.camera_tilt -= 0.01f;
            global_state.needs_render = true;
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Tilt camera down (if paused)");

    // Row 5: Visualization adjustments
    ImGui::Text("Visualization Adjustments");
    if (ImGui::Button("A##DecScale")) {
        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                double current_scale = eq->getVisScale();
                current_scale /= 1.1;
                eq->setVisScale(current_scale);
                sim.compute();
                global_state.needs_render = true;
            }
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Decrease vis scale by 10%");
    ImGui::SameLine();
    if (ImGui::Button("S##IncScale")) {
        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                double current_scale = eq->getVisScale();
                current_scale *= 1.1;
                eq->setVisScale(current_scale);
                sim.compute();
                global_state.needs_render = true;
            }
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Increase vis scale by 10%");
    ImGui::SameLine();
    if (ImGui::Button("Q##DecIntensity")) {
        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                double current_intensity = eq->getVisColorIntensity();
                current_intensity /= 1.1;
                eq->setVisColorIntensity(current_intensity);
                sim.compute();
                global_state.needs_render = true;
            }
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Decrease color intensity by 10%");
    ImGui::SameLine();
    if (ImGui::Button("W##IncIntensity")) {
        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                double current_intensity = eq->getVisColorIntensity();
                current_intensity *= 1.1;
                eq->setVisColorIntensity(current_intensity);
                sim.compute();
                global_state.needs_render = true;
            }
        }
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Increase color intensity by 10%");

    // Row 6: Grid size and reset
    ImGui::Text("Grid and Reset");
    if (ImGui::Button("X##IncGrid")) {
        sim.size = std::min(sim.size + 1, 50);
        sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
        sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
        sim.computed_scalar.resize(sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L)));
        sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<long double>>>(
            sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L))));
        sim.setup();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Increase grid size");
    ImGui::SameLine();
    if (ImGui::Button("Z##DecGrid")) {
        sim.size = std::max(sim.size - 1, 5);
        sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
        sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
        sim.computed_scalar.resize(sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L)));
        sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<long double>>>(
            sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L))));
        sim.setup();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Decrease grid size");
    ImGui::SameLine();
    if (ImGui::Button("R##ResetSim")) {
        global_state.camera_angle = 0.0f;
        global_state.camera_tilt = 0.0f;
        global_state.camera_zoom = 0.3f;
        sim.t = 0.0;
        sim.setup();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Reset simulation");
    ImGui::SameLine();
    if (ImGui::Button("P##ResetParticles")) {
        sim.particles.clear();
        sim.setup();
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Reset particles");

    // Row 7: Fullscreen and quit
    ImGui::Text("Window and Exit");
    if (ImGui::Button("F##Fullscreen")) {
        global_state.is_fullscreen = !global_state.is_fullscreen;
        if (global_state.is_fullscreen) {
            if (SDL_SetWindowFullscreen(window_, SDL_WINDOW_FULLSCREEN) < 0) {
                log_message("Failed to switch to fullscreen: " + std::string(SDL_GetError()));
                global_state.is_fullscreen = false;
            }
        } else {
            if (SDL_SetWindowFullscreen(window_, 0) < 0) {
                log_message("Failed to switch to windowed mode: " + std::string(SDL_GetError()));
                global_state.is_fullscreen = true;
            } else {
                SDL_SetWindowSize(window_, 800, 600); // 4:3 resolution
                SDL_SetWindowPosition(window_, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
            }
        }
        global_state.needs_render = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Toggle fullscreen");
    ImGui::SameLine();
    if (ImGui::Button("Esc##Quit")) {
        SDL_Event quit_event;
        quit_event.type = SDL_EVENT_QUIT;
        SDL_PushEvent(&quit_event);
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Quit application");

    ImGui::End();
}

void Menu::handle_event(const SDL_Event& e, Simulation& sim) {
    // Pass events to ImGui
    ImGui_ImplSDL3_ProcessEvent(&e);
}