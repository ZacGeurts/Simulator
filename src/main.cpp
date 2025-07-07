#include "globals.h"
#include "renderer.h"
#include "menu.h"
#include <SDL3/SDL.h>
#include <iostream>

int main() {
    log_message("Program started");
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        log_message("SDL initialization failed: " + std::string(SDL_GetError()));
        return 1;
    }

    // Create window with 4:3 resolution (800x600), resizable
    SDL_Window* window = SDL_CreateWindow("Physics Visualization", 800, 600, SDL_WINDOW_VULKAN | SDL_WINDOW_RESIZABLE);
    if (!window) {
        log_message("Window creation failed: " + std::string(SDL_GetError()));
        SDL_Quit();
        return 1;
    }
    log_message("Window created successfully");

    VulkanContext ctx;
    {
        std::lock_guard<std::mutex> lock(global_state.sim_mutex);
        global_state.sim.setup();
        log_message("Simulation setup completed");
        if (!init_renderer(window, ctx, global_state.sim)) {
            log_message("Renderer initialization failed");
            SDL_DestroyWindow(window);
            SDL_Quit();
            return 1;
        }
    }
    log_message("Renderer initialized successfully");

    Menu menu;
    if (!menu.init(window, global_state.sim)) {
        log_message("Menu initialization failed");
        cleanup_renderer(ctx);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }
    log_message("Menu initialized successfully");

    global_state.timer_id = SDL_AddTimer(16, timer, nullptr);
    if (global_state.timer_id == 0) {
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
            menu.handle_event(e, global_state.sim);
            if (e.type == SDL_EVENT_QUIT) {
                quit = true;
            } else if (e.type == SDL_EVENT_WINDOW_RESIZED || e.type == SDL_EVENT_WINDOW_PIXEL_SIZE_CHANGED) {
                int width, height;
                SDL_GetWindowSize(window, &width, &height);
                log_message("Window resized to " + std::to_string(width) + "x" + std::to_string(height));
                create_swapchain(window, ctx); // Changed to create_swapchain
                global_state.needs_render = true;
            } else if (e.type == SDL_EVENT_KEY_DOWN) {
                std::lock_guard<std::mutex> lock(global_state.sim_mutex);
                switch (e.key.scancode) {
                    case SDL_SCANCODE_ESCAPE:
                        quit = true;
                        break;
                    case SDL_SCANCODE_1:
                        global_state.current_display_mode = global_state.DisplayMode::POINTS;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_2:
                        global_state.current_display_mode = global_state.DisplayMode::ISOSURFACE;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_3:
                        global_state.current_display_mode = global_state.DisplayMode::WIREFRAME;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_4:
                        global_state.current_display_mode = global_state.DisplayMode::PARTICLES;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_5:
                        global_state.current_display_mode = global_state.DisplayMode::HYBRID;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_6:
                        global_state.current_display_mode = global_state.DisplayMode::SURFACE;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_7:
                        global_state.current_display_mode = global_state.DisplayMode::SPHERE_POINTS;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_LEFT:
                        global_state.sim.prev_equation();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_RIGHT:
                        global_state.sim.next_equation();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_R:
                        global_state.camera_angle = 0.0f;
                        global_state.camera_tilt = 0.0f;
                        global_state.camera_zoom = 0.3f;
                        global_state.sim.t = 0.0;
                        global_state.sim.setup();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_P:
                        global_state.sim.particles.clear();
                        global_state.sim.setup();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_SPACE:
                        global_state.is_paused = !global_state.is_paused;
                        global_state.rotate_camera = !global_state.is_paused;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_RETURN:
                        if (global_state.is_paused) {
                            global_state.sim.step();
                            global_state.needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_C:
                        if (global_state.is_paused) {
                            global_state.camera_angle -= 0.01f;
                            global_state.needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_V:
                        if (global_state.is_paused) {
                            global_state.camera_angle += 0.01f;
                            global_state.needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_B:
                        if (global_state.is_paused) {
                            global_state.camera_tilt += 0.01f;
                            global_state.needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_G:
                        if (global_state.is_paused) {
                            global_state.camera_tilt -= 0.01f;
                            global_state.needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_0:
                        global_state.camera_angle = 0.0f;
                        global_state.camera_tilt = 0.0f;
                        global_state.camera_zoom = 0.3f;
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_UP:
                        global_state.camera_zoom *= 0.95f;
                        global_state.camera_zoom = std::max(global_state.camera_zoom, 0.001f);
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_DOWN:
                        global_state.camera_zoom *= 1.05f;
                        global_state.camera_zoom = std::min(global_state.camera_zoom, 5.0f);
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_X:
                        global_state.sim.size = std::min(global_state.sim.size + 1, 50);
                        global_state.sim.ricci.resize(global_state.sim.size, std::vector<std::vector<Tensor>>(global_state.sim.size, std::vector<Tensor>(global_state.sim.size)));
                        global_state.sim.divergence.resize(global_state.sim.size, std::vector<std::vector<Tensor>>(global_state.sim.size, std::vector<Tensor>(global_state.sim.size)));
                        global_state.sim.computed_scalar.resize(global_state.sim.size, std::vector<std::vector<long double>>(global_state.sim.size, std::vector<long double>(global_state.sim.size, 0.0L)));
                        global_state.sim.scalar_4d.resize(global_state.sim.size, std::vector<std::vector<std::vector<long double>>>(
                            global_state.sim.size, std::vector<std::vector<long double>>(global_state.sim.size, std::vector<long double>(global_state.sim.size, 0.0L))));
                        global_state.sim.setup();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_Z:
                        global_state.sim.size = std::max(global_state.sim.size - 1, 5);
                        global_state.sim.ricci.resize(global_state.sim.size, std::vector<std::vector<Tensor>>(global_state.sim.size, std::vector<Tensor>(global_state.sim.size)));
                        global_state.sim.divergence.resize(global_state.sim.size, std::vector<std::vector<Tensor>>(global_state.sim.size, std::vector<Tensor>(global_state.sim.size)));
                        global_state.sim.computed_scalar.resize(global_state.sim.size, std::vector<std::vector<long double>>(global_state.sim.size, std::vector<long double>(global_state.sim.size, 0.0L)));
                        global_state.sim.scalar_4d.resize(global_state.sim.size, std::vector<std::vector<std::vector<long double>>>(
                            global_state.sim.size, std::vector<std::vector<long double>>(global_state.sim.size, std::vector<long double>(global_state.sim.size, 0.0L))));
                        global_state.sim.setup();
                        global_state.needs_render = true;
                        break;
                    case SDL_SCANCODE_S:
                        if (!global_state.sim.equations.empty() && global_state.sim.current_equation < global_state.sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(global_state.sim.equations[global_state.sim.current_equation].get())) {
                                double current_scale = eq->getVisScale();
                                current_scale *= 1.1;
                                eq->setVisScale(current_scale);
                                global_state.sim.compute();
                                global_state.needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_A:
                        if (!global_state.sim.equations.empty() && global_state.sim.current_equation < global_state.sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(global_state.sim.equations[global_state.sim.current_equation].get())) {
                                double current_scale = eq->getVisScale();
                                current_scale /= 1.1;
                                eq->setVisScale(current_scale);
                                global_state.sim.compute();
                                global_state.needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_W:
                        if (!global_state.sim.equations.empty() && global_state.sim.current_equation < global_state.sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(global_state.sim.equations[global_state.sim.current_equation].get())) {
                                double current_intensity = eq->getVisColorIntensity();
                                current_intensity *= 1.1;
                                eq->setVisColorIntensity(current_intensity);
                                global_state.sim.compute();
                                global_state.needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_Q:
                        if (!global_state.sim.equations.empty() && global_state.sim.current_equation < global_state.sim.equations.size()) {
                            if (auto* eq = dynamic_cast<CustomEquation*>(global_state.sim.equations[global_state.sim.current_equation].get())) {
                                double current_intensity = eq->getVisColorIntensity();
                                current_intensity /= 1.1;
                                eq->setVisColorIntensity(current_intensity);
                                global_state.sim.compute();
                                global_state.needs_render = true;
                            }
                        }
                        break;
                    case SDL_SCANCODE_F:
                        global_state.is_fullscreen = !global_state.is_fullscreen;
                        if (global_state.is_fullscreen) {
                            if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN) < 0) {
                                log_message("Failed to switch to fullscreen: " + std::string(SDL_GetError()));
                                global_state.is_fullscreen = false;
                            }
                        } else {
                            if (SDL_SetWindowFullscreen(window, 0) < 0) {
                                log_message("Failed to switch to windowed mode: " + std::string(SDL_GetError()));
                                global_state.is_fullscreen = true;
                            } else {
                                SDL_SetWindowSize(window, 800, 600); // Reset to 4:3
                                SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
                            }
                        }
                        create_swapchain(window, ctx); // Changed to create_swapchain
                        global_state.needs_render = true;
                        break;
                    default:
                        break;
                }
            }
        }

        {
            std::lock_guard<std::mutex> lock(global_state.sim_mutex);
            render(global_state.sim, ctx, menu); // Ensure render is defined in renderer.h
        }
    }

    if (global_state.timer_id != 0) {
        SDL_RemoveTimer(global_state.timer_id);
        log_message("Timer removed");
    }
    cleanup_renderer(ctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    log_message("Program exited");
    return 0;
}