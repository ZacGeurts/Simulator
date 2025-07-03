#include "globals.h"
#include "renderer.h"
#include "menu.h"
#include <SDL3/SDL.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <mutex>

std::mutex sim_mutex;
Simulation sim;
bool rotate_camera = true;
bool is_paused = false;
bool needs_render = true;
bool is_fullscreen = true;
SDL_TimerID timer_id = 0;

void cleanup_ttf() {
    // Placeholder: No TTF used in ImGui menu
}

void log_message(const std::string& message) {
    static std::set<std::string> logged_messages;
    std::ofstream log("log.txt", std::ios::app);
    if (log.is_open() && logged_messages.find(message) == logged_messages.end()) {
        std::time_t now = std::time(nullptr);
        log << "[" << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "] " << message << "\n";
        logged_messages.insert(message);
        log.close();
    }
}

Uint32 timer(Uint32 interval, void* param) {
    {
        std::lock_guard<std::mutex> lock(sim_mutex);
        if (!is_paused) {
            sim.step();
            if (rotate_camera) {
                camera_angle += 0.05f;
                needs_render = true;
            }
        }
    }
    return 16;
}

int main() {
    log_message("Program started");
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0) {
        log_message("SDL initialization failed: " + std::string(SDL_GetError()));
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("Physics Visualization", 800, 600, SDL_WINDOW_VULKAN | SDL_WINDOW_FULLSCREEN_DESKTOP | SDL_WINDOW_RESIZABLE);
    if (!window) {
        log_message("Window creation failed: " + std::string(SDL_GetError()));
        SDL_Quit();
        return 1;
    }
    log_message("Window created successfully");

    VulkanContext ctx;
    {
        std::lock_guard<std::mutex> lock(sim_mutex);
        sim.setup();
        log_message("Simulation setup completed");
        if (!init_renderer(window, ctx, sim)) {
            log_message("Renderer initialization failed");
            SDL_DestroyWindow(window);
            SDL_Quit();
            return 1;
        }
    }
    log_message("Renderer initialized successfully");

    Menu menu;
    if (!menu.init(window, sim)) {
        log_message("Menu initialization failed");
        cleanup_renderer(ctx);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    log_message("Menu initialized successfully");

    timer_id = SDL_AddTimer(16, timer, nullptr);
    if (timer_id == 0) {
        log_message("Failed to create timer: " + std::string(SDL_GetError()));
        cleanup_renderer(ctx);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_Event e;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&e)) {
            menu.handle_event(e, sim);
            if (e.type == SDL_EVENT_QUIT) {
                quit = true;
            } else if (e.type == SDL_EVENT_KEY_DOWN) {
                std::lock_guard<std::mutex> lock(sim_mutex);
                switch (e.key.keysym.scancode) {
                    case SDL_SCANCODE_ESCAPE:
                        quit = true;
                        break;
                    case SDL_SCANCODE_1:
                        current_display_mode = POINTS;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_2:
                        current_display_mode = ISOSURFACE;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_3:
                        current_display_mode = WIREFRAME;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_4:
                        current_display_mode = PARTICLES;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_5:
                        current_display_mode = HYBRID;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_6:
                        current_display_mode = SURFACE;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_7:
                        current_display_mode = SPHERE_POINTS;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_LEFT:
                        sim.prev_equation();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_RIGHT:
                        sim.next_equation();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_R:
                        camera_angle = 0.0f;
                        camera_tilt = 0.0f;
                        camera_zoom = 0.3f;
                        sim.t = 0.0;
                        sim.setup();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_P:
                        sim.particles.clear();
                        sim.setup();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_SPACE:
                        is_paused = !is_paused;
                        rotate_camera = !is_paused;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_RETURN:
                        if (is_paused) {
                            sim.step();
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_C:
                        if (is_paused) {
                            camera_angle -= 0.01f;
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_V:
                        if (is_paused) {
                            camera_angle += 0.01f;
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_B:
                        if (is_paused) {
                            camera_tilt += 0.01f;
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_G:
                        if (is_paused) {
                            camera_tilt -= 0.01f;
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_0:
                        camera_angle = 0.0f;
                        camera_tilt = 0.0f;
                        camera_zoom = 0.3f;
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_UP:
                        camera_zoom *= 0.95f;
                        camera_zoom = std::max(camera_zoom, 0.001f);
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_DOWN:
                        camera_zoom *= 1.05f;
                        camera_zoom = std::min(camera_zoom, 5.0f);
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_X:
                        sim.size = std::min(sim.size + 1, 50);
                        sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                        sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                        sim.computed_scalar.resize(sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L)));
                        sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<long double>>>(
                            sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L))));
                        sim.setup();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_Z:
                        sim.size = std::max(sim.size - 1, 5);
                        sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                        sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                        sim.computed_scalar.resize(sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L)));
                        sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<long double>>>(
                            sim.size, std::vector<std::vector<long double>>(sim.size, std::vector<long double>(sim.size, 0.0L))));
                        sim.setup();
                        needs_render = true;
                        break;
                    case SDL_SCANCODE_S:
                        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                                double current_scale = eq->getVisScale();
                                current_scale *= 1.1;
                                eq->setVisScale(current_scale);
                                sim.compute();
                                needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_A:
                        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                                double current_scale = eq->getVisScale();
                                current_scale /= 1.1;
                                eq->setVisScale(current_scale);
                                sim.compute();
                                needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_W:
                        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                                double current_intensity = eq->getVisColorIntensity();
                                current_intensity *= 1.1;
                                eq->setVisColorIntensity(current_intensity);
                                sim.compute();
                                needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_Q:
                        if (!sim.equations.empty() && sim.current_equation < sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())) {
                                double current_intensity = eq->getVisColorIntensity();
                                current_intensity /= 1.1;
                                eq->setVisColorIntensity(current_intensity);
                                sim.compute();
                                needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_F:
                        is_fullscreen = !is_fullscreen;
                        if (is_fullscreen) {
                            if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP) < 0) {
                                log_message("Failed to switch to fullscreen: " + std::string(SDL_GetError()));
                                is_fullscreen = false;
                            }
                        } else {
                            if (SDL_SetWindowFullscreen(window, 0) < 0) {
                                log_message("Failed to switch to windowed mode: " + std::string(SDL_GetError()));
                                is_fullscreen = true;
                            } else {
                                SDL_SetWindowSize(window, 800, 600);
                                SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
                            }
                        }
                        needs_render = true;
                        break;
                    default:
                        break;
                }
            }
        }

        {
            std::lock_guard<std::mutex> lock(sim_mutex);
            render(sim, ctx, menu);
        }
    }

    if (timer_id != 0) {
        SDL_RemoveTimer(timer_id);
        log_message("Timer removed");
    }
    cleanup_renderer(ctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    log_message("Program exited");
    return 0;
}