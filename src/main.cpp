#include "renderer.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mutex>

std::mutex sim_mutex;

// Utility function to log messages
void log_message(const std::string& message) {
    std::ofstream log("log.txt", std::ios::app);
    if (log.is_open()) {
        std::time_t now = std::time(nullptr);
        log << "[" << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "] " << message << "\n";
        log.close();
    }
}

Simulation sim;
DisplayMode current_display_mode = POINTS;
float camera_angle = 0.0f;
float camera_zoom = 0.3f;
bool rotate_camera = true;
bool needs_render = false;
SDL_TimerID timer_id = 0;

// Timer function to update simulation and camera
Uint32 timer(Uint32 interval, void* param) {
    {
        std::lock_guard<std::mutex> lock(sim_mutex);
        sim.step();
        if (rotate_camera) {
            camera_angle += 0.05f;
        }
        needs_render = true;
    }
    return 16;
}

// Declare cleanup_ttf (defined in renderer.cpp)
void cleanup_ttf();

int main() {
    log_message("Program started");
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0) {
        log_message("SDL initialization failed: " + std::string(SDL_GetError()));
        std::cout << "SDL initialization failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    // Set OpenGL attributes
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    SDL_Window* window = SDL_CreateWindow(
        "Physics Visualization",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        800, 600,
        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN
    );
    if (!window) {
        log_message("Window creation failed: " + std::string(SDL_GetError()));
        std::cout << "Window creation failed: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    if (!gl_context) {
        log_message("GL context creation failed: " + std::string(SDL_GetError()));
        std::cout << "GL context creation failed: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    if (SDL_GL_MakeCurrent(window, gl_context) < 0) {
        log_message("Failed to make GL context current: " + std::string(SDL_GetError()));
        std::cout << "Failed to make GL context current: " << SDL_GetError() << std::endl;
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    SDL_GL_SetSwapInterval(1);

    {
        std::lock_guard<std::mutex> lock(sim_mutex);
        sim.setup();
        log_message("Simulation setup completed");
    }

    timer_id = SDL_AddTimer(16, timer, nullptr);
    if (timer_id == 0) {
        log_message("Failed to create timer: " + std::string(SDL_GetError()));
        std::cout << "Failed to create timer: " << SDL_GetError() << std::endl;
        SDL_GL_DeleteContext(gl_context);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_Event e;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
            if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                    case SDLK_ESCAPE:
                        quit = true;
                        break;
                    case SDLK_1:
                        current_display_mode = POINTS;
                        log_message("Display mode set to POINTS");
                        needs_render = true;
                        break;
                    case SDLK_2:
                        current_display_mode = ISOSURFACE;
                        log_message("Display mode set to ISOSURFACE");
                        needs_render = true;
                        break;
                    case SDLK_3:
                        current_display_mode = WIREFRAME;
                        log_message("Display mode set to WIREFRAME");
                        needs_render = true;
                        break;
                    case SDLK_4:
                        current_display_mode = PARTICLES;
                        log_message("Display mode set to PARTICLES");
                        needs_render = true;
                        break;
                    case SDLK_5:
                        current_display_mode = HYBRID;
                        log_message("Display mode set to HYBRID");
                        needs_render = true;
                        break;
                    case SDLK_LEFT:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.prev_equation();
                            log_message("Switched to previous equation");
                            needs_render = true;
                        }
                        break;
                    case SDLK_RIGHT:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.next_equation();
                            log_message("Switched to next equation");
                            needs_render = true;
                        }
                        break;
                    case SDLK_r:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_angle = 0.0f;
                            camera_zoom = 0.3f;
                            sim.t = 0.0;
                            sim.setup();
                            log_message("Simulation reset");
                            needs_render = true;
                        }
                        break;
                    case SDLK_p:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.particles.clear();
                            sim.setup();
                            log_message("Particles cleared and simulation reset");
                            needs_render = true;
                        }
                        break;
                    case SDLK_SPACE:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            rotate_camera = !rotate_camera;
                            log_message(rotate_camera ? "Camera rotation resumed" : "Camera rotation paused");
                            needs_render = true;
                        }
                        break;
                    case SDLK_UP:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_zoom *= 0.95f;
                            camera_zoom = std::max(camera_zoom, 0.05f);
                            log_message("Zoomed in to " + std::to_string(camera_zoom));
                            needs_render = true;
                        }
                        break;
                    case SDLK_DOWN:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_zoom *= 1.05f;
                            camera_zoom = std::min(camera_zoom, 5.0f);
                            log_message("Zoomed out to " + std::to_string(camera_zoom));
                            needs_render = true;
                        }
                        break;
                    case SDLK_z:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.size = std::min(sim.size + 1, 50);
                            sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                            sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                            sim.scalar.resize(sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0)));
                            sim.computed_scalar.resize(sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0)));
                            sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<double>>>(
                                sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0))));
                            sim.setup();
                            log_message("Simulation resized to " + std::to_string(sim.size));
                            needs_render = true;
                        }
                        break;
                    case SDLK_x:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.size = std::max(sim.size - 1, 5);
                            sim.ricci.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                            sim.divergence.resize(sim.size, std::vector<std::vector<Tensor>>(sim.size, std::vector<Tensor>(sim.size)));
                            sim.scalar.resize(sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0)));
                            sim.computed_scalar.resize(sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0)));
                            sim.scalar_4d.resize(sim.size, std::vector<std::vector<std::vector<double>>>(
                                sim.size, std::vector<std::vector<double>>(sim.size, std::vector<double>(sim.size, 0.0))));
                            sim.setup();
                            log_message("Simulation resized to " + std::to_string(sim.size));
                            needs_render = true;
                        }
                        break;
                    case SDLK_c:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_angle = 0.0f;
                            camera_zoom = 0.3f;
                            log_message("Camera centered");
                            needs_render = true;
                        }
                        break;
                    default:
                        break;
                }
            }
        }

        // Render if needed
        if (needs_render) {
            {
                std::lock_guard<std::mutex> lock(sim_mutex);
                render(sim);
            }
            SDL_GL_SwapWindow(window);
            needs_render = false;
        }
    }

    // Cleanup
    if (timer_id != 0) {
        SDL_RemoveTimer(timer_id);
        log_message("Timer removed");
    }
    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    cleanup_ttf();
    SDL_Quit();
    log_message("Program exited");
    return 0;
}