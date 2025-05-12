#include "renderer.h"
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mutex>
#include <set>

std::mutex sim_mutex;

// Utility function to log messages, ensuring each unique message is logged only once
void log_message(const std::string& message) {
    static std::set<std::string> logged_messages;
    std::ofstream log("log.txt", std::ios::app);
    if (log.is_open()) {
        if (logged_messages.find(message) == logged_messages.end()) {
            std::time_t now = std::time(nullptr);
            log << "[" << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "] " << message << "\n";
            logged_messages.insert(message);
        }
        log.close();
    }
}

Simulation sim;
DisplayMode current_display_mode = POINTS;
float camera_angle = 0.0f;
float camera_tilt = 0.0f;
float camera_zoom = 0.3f;
bool rotate_camera = true;
bool is_paused = false;
bool needs_render = false;
bool is_fullscreen = false;
SDL_TimerID timer_id = 0;

// Timer function to update simulation and camera
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
        SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE
    );
    if (!window) {
        log_message("Window creation failed: " + std::string(SDL_GetError()));
        std::cout << "Window creation failed: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }
    log_message("Window created successfully");

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
    					current_display_mode = SPHERE_POINTS; // or 7
						needs_render = true;
                        break;
                    case SDL_SCANCODE_LEFT:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.prev_equation();
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_RIGHT:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.next_equation();
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_R:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_angle = 0.0f;
                            camera_tilt = 0.0f;
                            camera_zoom = 0.3f;
                            sim.t = 0.0;
                            sim.setup();
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_P:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            sim.particles.clear();
                            sim.setup();
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_SPACE:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            is_paused = !is_paused;
                            rotate_camera = !is_paused;
                            //log_message(is_paused ? "Simulation and camera rotation paused" : "Simulation and camera rotation unpaused");
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_RETURN:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (is_paused) {
                                sim.step();
                                needs_render = true;
                                //log_message("Simulation advanced one frame");
                            }
                        }
                        break;
                    case SDL_SCANCODE_C:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (is_paused) {
                                camera_angle -= 0.01f;
                                needs_render = true;
                                //log_message("Camera rotated left");
                            }
                        }
                        break;
                    case SDL_SCANCODE_V:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (is_paused) {
                                camera_angle += 0.01f;
                                needs_render = true;
                                //log_message("Camera rotated right");
                            }
                        }
                        break;
                    case SDL_SCANCODE_B:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (is_paused) {
                                camera_tilt += 0.01f;
                                needs_render = true;
                                //log_message("Camera tilted up");
                            }
                        }
                        break;
                    case SDL_SCANCODE_G:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (is_paused) {
                                camera_tilt -= 0.01f;
                                needs_render = true;
                                //log_message("Camera tilted down");
                            }
                        }
                        break;
                    case SDL_SCANCODE_0:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_angle = 0.0f;
                            camera_tilt = 0.0f;
                            camera_zoom = 0.3f;
                            needs_render = true;
                            //log_message("Camera reset");
                        }
                        break;
                    case SDL_SCANCODE_UP:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_zoom *= 0.95f;
                            camera_zoom = std::max(camera_zoom, 0.001f);
                            needs_render = true;
                            //log_message("Zoomed in, camera_zoom = " + std::to_string(camera_zoom));
                        }
                        break;
                    case SDL_SCANCODE_DOWN:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            camera_zoom *= 1.05f;
                            camera_zoom = std::min(camera_zoom, 5.0f);
                            needs_render = true;
                            //log_message("Zoomed out, camera_zoom = " + std::to_string(camera_zoom));
                        }
                        break;
                    case SDL_SCANCODE_X:
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
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_Z:
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
                            needs_render = true;
                        }
                        break;
                    case SDL_SCANCODE_S:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (!sim.equations.empty()) {
                                if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
                                    double current_scale = eq->getVisScale();
                                    current_scale *= 1.1;
                                    eq->setVisScale(current_scale);
                                    sim.compute();
                                    needs_render = true;
                                }
                            }
                        }
                        break;
                    case SDL_SCANCODE_A:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (!sim.equations.empty()) {
                                if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
                                    double current_scale = eq->getVisScale();
                                    current_scale /= 1.1;
                                    eq->setVisScale(current_scale);
                                    sim.compute();
                                    needs_render = true;
                                }
                            }
                        }
                        break;
                    case SDL_SCANCODE_W:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (!sim.equations.empty()) {
                                if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
                                    double current_intensity = eq->getVisColorIntensity();
                                    current_intensity *= 1.1;
                                    eq->setVisColorIntensity(current_intensity);
                                    sim.compute();
                                    needs_render = true;
                                }
                            }
                        }
                        break;
                    case SDL_SCANCODE_Q:
                        {
                            std::lock_guard<std::mutex> lock(sim_mutex);
                            if (!sim.equations.empty()) {
                                if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
                                    double current_intensity = eq->getVisColorIntensity();
                                    current_intensity /= 1.1;
                                    eq->setVisColorIntensity(current_intensity);
                                    sim.compute();
                                    needs_render = true;
                                }
                            }
                        }
                        break;
                    case SDL_SCANCODE_F:
                    {
                        is_fullscreen = !is_fullscreen;
                        const int DEFAULT_WIDTH = 800;
                        const int DEFAULT_HEIGHT = 600;
                        const float TARGET_ASPECT = 4.0f / 3.0f;

                        if (is_fullscreen) {
                            if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP) < 0) {
                                log_message("Failed to switch to fullscreen: " + std::string(SDL_GetError()));
                                is_fullscreen = false;
                                break;
                            }

                            int drawable_w, drawable_h;
                            SDL_GL_GetDrawableSize(window, &drawable_w, &drawable_h);
                            float window_aspect = (float)drawable_w / drawable_h;

                            int vp_width, vp_height, vp_x, vp_y;
                            if (window_aspect > TARGET_ASPECT) {
                                vp_height = drawable_h;
                                vp_width = (int)(vp_height * TARGET_ASPECT);
                                vp_x = (drawable_w - vp_width) / 2;
                                vp_y = 0;
                            } else {
                                vp_width = drawable_w;
                                vp_height = (int)(vp_width / TARGET_ASPECT);
                                vp_x = 0;
                                vp_y = (drawable_h - vp_height) / 2;
                            }

                            glViewport(vp_x, vp_y, vp_width, vp_height);
                        } else {
                            if (SDL_SetWindowFullscreen(window, 0) < 0) {
                                log_message("Failed to switch to windowed mode: " + std::string(SDL_GetError()));
                                is_fullscreen = true;
                                break;
                            }

                            SDL_SetWindowSize(window, DEFAULT_WIDTH, DEFAULT_HEIGHT);
                            SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
                            glViewport(0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT);
                        }
                        needs_render = true;
                    }
                    break;
                    default:
                        break;
                }
            }
        }

        if (needs_render) {
            const int DEFAULT_WIDTH = 800;
            const int DEFAULT_HEIGHT = 600;
            const float TARGET_ASPECT = 4.0f / 3.0f;

            if (is_fullscreen) {
                int drawable_w, drawable_h;
                SDL_GL_GetDrawableSize(window, &drawable_w, &drawable_h);
                float window_aspect = (float)drawable_w / drawable_h;

                int vp_width, vp_height, vp_x, vp_y;
                if (window_aspect > TARGET_ASPECT) {
                    vp_height = drawable_h;
                    vp_width = (int)(vp_height * TARGET_ASPECT);
                    vp_x = (drawable_w - vp_width) / 2;
                    vp_y = 0;
                } else {
                    vp_width = drawable_w;
                    vp_height = (int)(vp_width / TARGET_ASPECT);
                    vp_x = 0;
                    vp_y = (drawable_h - vp_height) / 2;
                }
                glViewport(vp_x, vp_y, vp_width, vp_height);
            } else {
                glViewport(0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT);
            }

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            {
                std::lock_guard<std::mutex> lock(sim_mutex);
                render(sim);
            }
            SDL_GL_SwapWindow(window);
            needs_render = false;
        }
    }

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