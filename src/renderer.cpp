#include "renderer.h"
#include "equations.h"
#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>
#include <SDL2/SDL_ttf.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>
#include "tables.h"

// Font and TTF management
static TTF_Font* font = nullptr;
static bool ttf_initialized = false;

// Check if file exists and is readable
bool file_exists(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    bool exists = file.good();
    if (!exists) {
        std::string warn = "Font file does not exist or is inaccessible: " + path;
        log_message(warn);
        std::cerr << warn << std::endl;
    } else {
        std::string info = "Font file exists and is readable: " + path;
        log_message(info);
        std::cerr << info << std::endl;
    }
    return exists;
}

// Initialize SDL_ttf
bool init_ttf() {
    std::string info = "Starting SDL_ttf initialization";
    log_message(info);
    std::cerr << info << std::endl;

    if (TTF_Init() == -1) {
        std::string error = "SDL_ttf initialization failed: " + std::string(TTF_GetError());
        log_message(error);
        std::cerr << error << std::endl;
        return false;
    }
    ttf_initialized = true;
    log_message("SDL_ttf initialized successfully");
    std::cerr << "SDL_ttf initialized successfully" << std::endl;

    const char* font_paths[] = {
        "../DejaVuSans.ttf",
        "DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "C:/Windows/Fonts/DejaVuSans.ttf"
    };

    for (const char* path : font_paths) {
        std::string attempt = "Attempting to load font at: " + std::string(path);
        log_message(attempt);
        std::cerr << attempt << std::endl;

        if (!file_exists(path)) {
            continue;
        }

        font = TTF_OpenFont(path, 24);
        if (font) {
            std::string success = "Successfully loaded font: " + std::string(path);
            log_message(success);
            std::cerr << success << std::endl;
            return true;
        }
        std::string error = "Failed to load font at " + std::string(path) + ": " + (TTF_GetError() ? TTF_GetError() : "Unknown error");
        log_message(error);
        std::cerr << error << std::endl;
    }

    std::string error = "Failed to load DejaVuSans.ttf from any path";
    log_message(error);
    std::cerr << error << std::endl;
    TTF_Quit();
    ttf_initialized = false;
    return false;
}

// Clean up SDL_ttf
void cleanup_ttf() {
    if (font) {
        TTF_CloseFont(font);
        font = nullptr;
    }
    if (ttf_initialized) {
        TTF_Quit();
        ttf_initialized = false;
    }
}

// Render text to texture
GLuint render_text_to_texture(const std::string& text, SDL_Color color, int& w, int& h) {
    if (!font || text.empty()) {
        std::string error = std::string("Cannot render text: ") + (text.empty() ? "Text is empty" : "Font not loaded");
        log_message(error);
        std::cerr << error << std::endl;
        w = 0;
        h = 0;
        return 0;
    }

    SDL_Surface* surface = TTF_RenderText_Blended(font, text.c_str(), color);
    if (!surface) {
        std::string error = "Failed to render text '" + text + "': " + (TTF_GetError() ? TTF_GetError() : "Unknown error");
        log_message(error);
        std::cerr << error << std::endl;
        w = 0;
        h = 0;
        return 0;
    }

    SDL_Surface* converted = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGBA32, 0);
    SDL_FreeSurface(surface);
    if (!converted) {
        std::string error = "Failed to convert surface to RGBA: " + std::string(SDL_GetError());
        log_message(error);
        std::cerr << error << std::endl;
        w = 0;
        h = 0;
        return 0;
    }

    w = converted->w;
    h = converted->h;

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, converted->w, converted->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, converted->pixels);

    SDL_FreeSurface(converted);
    glBindTexture(GL_TEXTURE_2D, 0);
    return texture;
}

void render(const Simulation& sim) {
    static bool ttf_init_done = false;
    if (!ttf_init_done) {
        if (!init_ttf()) {
            log_message("Failed to initialize TTF. Text rendering disabled.");
            std::cerr << "Failed to initialize TTF. Text rendering disabled." << std::endl;
        }
        ttf_init_done = true;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 800.0 / 600.0, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    float cam_dist = 5.0f / camera_zoom;
    float cam_x = cam_dist * std::cos(camera_angle);
    float cam_y = cam_dist * std::sin(camera_angle);
    float cam_z = cam_dist * 0.5f;
    float tilted_y = cam_y * std::cos(camera_tilt) - cam_z * std::sin(camera_tilt);
    float tilted_z = cam_y * std::sin(camera_tilt) + cam_z * std::cos(camera_tilt);
    gluLookAt(cam_x, tilted_y, tilted_z, 0, 0, 0, 0, 0, 1);

    double max_scalar = 1e-10, min_scalar = 1e10;
    double max_4d = 1e-10, min_4d = 1e10;
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                max_scalar = std::max(max_scalar, sim.computed_scalar[i][j][k]);
                min_scalar = std::min(min_scalar, sim.computed_scalar[i][j][k]);
                for (int w = 0; w < sim.size; ++w) {
                    max_4d = std::max(max_4d, sim.scalar_4d[i][j][k][w]);
                    min_4d = std::min(min_4d, sim.scalar_4d[i][j][k][w]);
                }
            }
        }
    }
    double range = max_scalar - min_scalar;
    if (range == 0) range = 1.0;
    double range_4d = max_4d - min_4d;
    if (range_4d == 0) range_4d = 1.0;
    double iso_level = min_scalar + range * (0.3 + 0.2 * std::sin(sim.t));

    int w_slice = static_cast<int>(sim.t * 2) % sim.size;

    float vis_scale = 1.0f;
    float vis_color_intensity = 1.0f;
    if (!sim.equations.empty()) {
        if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
            vis_scale = static_cast<float>(eq->getVisScale());
            vis_color_intensity = static_cast<float>(eq->getVisColorIntensity());
        }
    }

    if (current_display_mode == PARTICLES || current_display_mode == HYBRID) {
        glPointSize(8.0);
        glBegin(GL_POINTS);
        for (const auto& p : sim.particles) {
            float speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
            float r = std::min(speed / 5.0f, 1.0f) * vis_color_intensity;
            float g = 0.5f * vis_color_intensity;
            float b = (1.0f - r) * vis_color_intensity;
            glColor3f(r, g, b);
            glVertex3f(p.x * vis_scale, p.y * vis_scale, p.z * vis_scale);
        }
        glEnd();
        glBegin(GL_LINES);
        glColor3f(0.5f, 0.5f, 1.0f);
        for (const auto& p : sim.particles) {
            glVertex3f(p.x * vis_scale, p.y * vis_scale, p.z * vis_scale);
            glVertex3f((p.x - p.vx * sim.dt * 10) * vis_scale, (p.y - p.vy * sim.dt * 10) * vis_scale, (p.z - p.vz * sim.dt * 10) * vis_scale);
        }
        glEnd();
    }

    if (current_display_mode != PARTICLES && sim.size > 0) {
        if (current_display_mode == POINTS || current_display_mode == HYBRID) {
            glPointSize(6.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        float x = static_cast<float>((i - sim.size / 2) * sim.dx);
                        float y = static_cast<float>((j - sim.size / 2) * sim.dx);
                        float z = static_cast<float>((k - sim.size / 2) * sim.dx);
                        float value = (sim.computed_scalar[i][j][k] - min_scalar) / range;
                        float r = value * vis_color_intensity;
                        float g = (1.0f - value) * vis_color_intensity;
                        float b = 0.5f * vis_color_intensity;
                        glColor3f(r, g, b);
                        glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                    }
                }
            }
            glEnd();
        }

        if (current_display_mode == ISOSURFACE) {
            static bool glew_initialized = false;
            if (!glew_initialized) {
                GLenum err = glewInit();
                if (err != GLEW_OK) {
                    log_message("GLEW initialization failed: " + std::string(reinterpret_cast<const char*>(glewGetErrorString(err))));
                    return;
                }
                glew_initialized = true;
            }

            if (sim.size <= 0 || w_slice < 0 || w_slice >= sim.size) {
                log_message("Invalid simulation size or w_slice: size=" + std::to_string(sim.size) + ", w_slice=" + std::to_string(w_slice));
                return;
            }

            // Precompute colors and normals with boundary handling
            std::vector<std::vector<std::vector<float>>> colors(sim.size, std::vector<std::vector<float>>(sim.size, std::vector<float>(sim.size)));
            std::vector<std::vector<std::vector<std::array<float, 3>>>> normals(sim.size, std::vector<std::vector<std::array<float, 3>>>(sim.size, std::vector<std::array<float, 3>>(sim.size)));
            #pragma omp parallel for
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        colors[i][j][k] = (sim.scalar_4d[i][j][k][w_slice] - min_4d) / range_4d;
                        float nx = 0, ny = 0, nz = 0;
                        // Use one-sided differences at boundaries
                        if (i == 0)
                            nx = (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i][j][k]) / sim.dx;
                        else if (i == sim.size - 1)
                            nx = (sim.computed_scalar[i][j][k] - sim.computed_scalar[i - 1][j][k]) / sim.dx;
                        else
                            nx = (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx);
                        if (j == 0)
                            ny = (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j][k]) / sim.dx;
                        else if (j == sim.size - 1)
                            ny = (sim.computed_scalar[i][j][k] - sim.computed_scalar[i][j - 1][k]) / sim.dx;
                        else
                            ny = (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx);
                        if (k == 0)
                            nz = (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k]) / sim.dx;
                        else if (k == sim.size - 1)
                            nz = (sim.computed_scalar[i][j][k] - sim.computed_scalar[i][j][k - 1]) / sim.dx;
                        else
                            nz = (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx);
                        float norm = std::sqrt(nx * nx + ny * ny + nz * nz);
                        normals[i][j][k] = {norm > 0 ? nx / norm : 0, norm > 0 ? ny / norm : 0, norm > 0 ? nz / norm : 0};
                    }
                }
            }

            std::vector<std::vector<float>> thread_vertices(omp_get_max_threads());
            static int cube_error_count = 0;
            static int edge_error_count = 0;
            static int degenerate_triangle_count = 0;
            const int max_errors_per_frame = 5;

            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < sim.size - 1; ++i) {
                int thread_id = omp_get_thread_num();
                thread_vertices[thread_id].reserve(1000 * 9);
                for (int j = 0; j < sim.size - 1; ++j) {
                    for (int k = 0; k < sim.size - 1; ++k) {
                        float values[8] = {
                            static_cast<float>(sim.computed_scalar[i][j][k]),
                            static_cast<float>(sim.computed_scalar[i + 1][j][k]),
                            static_cast<float>(sim.computed_scalar[i + 1][j + 1][k]),
                            static_cast<float>(sim.computed_scalar[i][j + 1][k]),
                            static_cast<float>(sim.computed_scalar[i][j][k + 1]),
                            static_cast<float>(sim.computed_scalar[i + 1][j][k + 1]),
                            static_cast<float>(sim.computed_scalar[i + 1][j + 1][k + 1]),
                            static_cast<float>(sim.computed_scalar[i][j + 1][k + 1])
                        };
                        float vertices[8][3] = {
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)}
                        };
                        int vi[8] = {i, i + 1, i + 1, i, i, i + 1, i + 1, i};
                        int vj[8] = {j, j, j + 1, j + 1, j, j, j + 1, j + 1};
                        int vk[8] = {k, k, k, k, k + 1, k + 1, k + 1, k + 1};

                        int cube_index = 0;
                        for (int v = 0; v < 8; ++v) {
                            if (values[v] >= iso_level) cube_index |= (1 << v);
                        }
                        if (cube_index < 0 || cube_index > 255) {
                            #pragma omp critical
                            {
                                if (cube_error_count < max_errors_per_frame) {
                                    log_message("Invalid cube_index: " + std::to_string(cube_index) + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k));
                                    cube_error_count++;
                                }
                            }
                            continue;
                        }

                        if (edgeMasks[cube_index] == 0) {
                            continue;
                        }

                        // Store triangle vertices to check for degeneracy
                        std::vector<std::array<float, 3>> triangle_vertices;
                        for (int t = 0; t < 16 && triTable[cube_index][t] != -1; t += 3) {
                            triangle_vertices.clear();
                            for (int v = 0; v < 3; ++v) {
                                int edge = triTable[cube_index][t + v];
                                if (edge < 0 || edge >= 12) {
                                    #pragma omp critical
                                    {
                                        if (edge_error_count < max_errors_per_frame) {
                                            log_message("Invalid edge index: " + std::to_string(edge) + " at cube_index=" + std::to_string(cube_index) + ", t=" + std::to_string(t) + ", i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k));
                                            edge_error_count++;
                                        }
                                    }
                                    triangle_vertices.clear();
                                    break;
                                }
                                int v0 = edgeTable[edge][0];
                                int v1 = edgeTable[edge][1];
                                float t0 = 0.5f;
                                // Add epsilon to avoid division by near-zero
                                float denom = values[v1] - values[v0];
                                if (std::abs(denom) > 1e-6) {
                                    t0 = (iso_level - values[v0]) / denom;
                                    t0 = std::max(0.0f, std::min(1.0f, t0));
                                }
                                float vert[3], norm[3];
                                float color = colors[vi[v0]][vj[v0]][vk[v0]] + t0 * (colors[vi[v1]][vj[v1]][vk[v1]] - colors[vi[v0]][vj[v0]][vk[v0]]);
                                for (int d = 0; d < 3; ++d) {
                                    vert[d] = vertices[v0][d] + t0 * (vertices[v1][d] - vertices[v0][d]);
                                    norm[d] = normals[vi[v0]][vj[v0]][vk[v0]][d] + t0 * (normals[vi[v1]][vj[v1]][vk[v1]][d] - normals[vi[v0]][vj[v0]][vk[v0]][d]);
                                }
                                float norm_len = std::sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
                                if (norm_len > 0) {
                                    norm[0] /= norm_len;
                                    norm[1] /= norm_len;
                                    norm[2] /= norm_len;
                                } else {
                                    norm[0] = norm[1] = norm[2] = 0.0f; // Fallback for zero normals
                                }
                                triangle_vertices.push_back({vert[0], vert[1], vert[2]});
                                thread_vertices[thread_id].insert(thread_vertices[thread_id].end(), {
                                    vert[0], vert[1], vert[2],
                                    norm[0], norm[1], norm[2],
                                    std::max(0.2f, color * vis_color_intensity), 0.5f, std::max(0.2f, (1.0f - color) * vis_color_intensity)
                                });
                            }
                            // Check for degenerate triangles
                            if (triangle_vertices.size() == 3) {
                                float v0[3] = {triangle_vertices[0][0], triangle_vertices[0][1], triangle_vertices[0][2]};
                                float v1[3] = {triangle_vertices[1][0], triangle_vertices[1][1], triangle_vertices[1][2]};
                                float v2[3] = {triangle_vertices[2][0], triangle_vertices[2][1], triangle_vertices[2][2]};
                                float edge1[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
                                float edge2[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
                                float cross[3] = {
                                    edge1[1] * edge2[2] - edge1[2] * edge2[1],
                                    edge1[2] * edge2[0] - edge1[0] * edge2[2],
                                    edge1[0] * edge2[1] - edge1[1] * edge2[0]
                                };
                                float area = 0.5f * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
                                if (area < 1e-6) {
                                    #pragma omp critical
                                    {
                                        if (degenerate_triangle_count < max_errors_per_frame) {
                                            log_message("Degenerate triangle at cube_index=" + std::to_string(cube_index) + ", i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k) + ", area=" + std::to_string(area));
                                            degenerate_triangle_count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            cube_error_count = 0;
            edge_error_count = 0;
            degenerate_triangle_count = 0;

            std::vector<float> vertex_data;
            size_t total_size = 0;
            for (const auto& tv : thread_vertices) {
                total_size += tv.size();
            }
            vertex_data.reserve(total_size);
            for (const auto& tv : thread_vertices) {
                vertex_data.insert(vertex_data.end(), tv.begin(), tv.end());
            }

            static GLuint vbo = 0;
            if (vbo == 0) {
                glGenBuffers(1, &vbo);
            }
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, vertex_data.size() * sizeof(float), vertex_data.data(), GL_STATIC_DRAW);

            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            GLfloat light_pos[] = {cam_x, cam_y, cam_z, 1.0f};
            GLfloat diffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
            GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
            glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glVertexPointer(3, GL_FLOAT, 9 * sizeof(float), (void*)0);
            glNormalPointer(GL_FLOAT, 9 * sizeof(float), (void*)(3 * sizeof(float)));
            glColorPointer(3, GL_FLOAT, 9 * sizeof(float), (void*)(6 * sizeof(float)));
            glDrawArrays(GL_TRIANGLES, 0, vertex_data.size() / 9);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            glDisable(GL_LIGHTING);
            glDisable(GL_COLOR_MATERIAL);
        }

        if (current_display_mode == SURFACE) {
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            GLfloat light_pos[] = {cam_x, cam_y, cam_z, 1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

            int k = sim.size / 2;
            glBegin(GL_QUADS);
            for (int i = 0; i < sim.size - 1; ++i) {
                for (int j = 0; j < sim.size - 1; ++j) {
                    float x0 = static_cast<float>((i - sim.size / 2) * sim.dx);
                    float y0 = static_cast<float>((j - sim.size / 2) * sim.dx);
                    float x1 = static_cast<float>((i + 1 - sim.size / 2) * sim.dx);
                    float y1 = static_cast<float>((j + 1 - sim.size / 2) * sim.dx);

                    float values[4] = {
                        static_cast<float>((sim.computed_scalar[i][j][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i + 1][j][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i + 1][j + 1][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i][j + 1][k] - min_scalar) / range)
                    };
                    float z[4] = {
                        static_cast<float>(values[0] * sim.dx * vis_scale),
                        static_cast<float>(values[1] * sim.dx * vis_scale),
                        static_cast<float>(values[2] * sim.dx * vis_scale),
                        static_cast<float>(values[3] * sim.dx * vis_scale)
                    };

                    float vec1[3] = {x1 - x0, y0 - y0, z[1] - z[0]};
                    float vec2[3] = {x0 - x0, y1 - y0, z[3] - z[0]};
                    float normal[3] = {
                        vec1[1] * vec2[2] - vec1[2] * vec2[1],
                        vec1[2] * vec2[0] - vec1[0] * vec2[2],
                        vec1[0] * vec2[1] - vec1[1] * vec2[0]
                    };
                    float norm_len = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
                    if (norm_len > 0) {
                        normal[0] /= norm_len;
                        normal[1] /= norm_len;
                        normal[2] /= norm_len;
                    }

                    glNormal3fv(normal);

                    glColor3f(values[0] * vis_color_intensity, 0.5f, (1.0f - values[0]) * vis_color_intensity);
                    glVertex3f(x0 * vis_scale, y0 * vis_scale, z[0]);

                    glColor3f(values[1] * vis_color_intensity, 0.5f, (1.0f - values[1]) * vis_color_intensity);
                    glVertex3f(x1 * vis_scale, y0 * vis_scale, z[1]);

                    glColor3f(values[2] * vis_color_intensity, 0.5f, (1.0f - values[2]) * vis_color_intensity);
                    glVertex3f(x1 * vis_scale, y1 * vis_scale, z[2]);

                    glColor3f(values[3] * vis_color_intensity, 0.5f, (1.0f - values[3]) * vis_color_intensity);
                    glVertex3f(x0 * vis_scale, y1 * vis_scale, z[3]);
                }
            }
            glEnd();

            glDisable(GL_LIGHTING);
            glDisable(GL_COLOR_MATERIAL);
        }

        if (current_display_mode == WIREFRAME) {
            glBegin(GL_LINES);
            glColor3f(0.5f, 0.5f, 1.0f);
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        float x = static_cast<float>((i - sim.size / 2) * sim.dx);
                        float y = static_cast<float>((j - sim.size / 2) * sim.dx);
                        float z = static_cast<float>((k - sim.size / 2) * sim.dx);
                        if (i < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f((x + sim.dx) * vis_scale, y * vis_scale, z * vis_scale);
                        }
                        if (j < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f(x * vis_scale, (y + sim.dx) * vis_scale, z * vis_scale);
                        }
                        if (k < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f(x * vis_scale, y * vis_scale, (z + sim.dx) * vis_scale);
                        }
                    }
                }
            }
            glEnd();
        }

        if (current_display_mode == SPHERE_POINTS) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glPointSize(12.0); // Larger points to reduce gaps
            glBegin(GL_POINTS);
            float radius = static_cast<float>(sim.size / 2) * sim.dx; // Sphere radius
            float radius_sq = radius * radius;
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        float x = static_cast<float>((i - sim.size / 2) * sim.dx);
                        float y = static_cast<float>((j - sim.size / 2) * sim.dx);
                        float z = static_cast<float>((k - sim.size / 2) * sim.dx);
                        // Check if point is within sphere
                        if (x * x + y * y + z * z <= radius_sq) {
                            float value = (sim.computed_scalar[i][j][k] - min_scalar) / range;
                            float r = value * vis_color_intensity;
                            float g = (1.0f - value) * vis_color_intensity;
                            float b = 0.5f * vis_color_intensity;
                            glColor4f(r, g, b, 0.5f); // Semi-transparent (alpha = 0.5)
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                        }
                    }
                }
            }
            glEnd();
            glDisable(GL_BLEND);
            glPointSize(1.0); // Reset point size
        }
    }

    glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(-1.0f * vis_scale, 1.0f * vis_scale, -1.0f * vis_scale);
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, 1.0f * vis_scale);
    glEnd();

    if (font && !sim.equations.empty()) {
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, 800, 0, 600, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glEnable(GL_TEXTURE_2D);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        std::string text = sim.equations[sim.current_equation]->name() + " | t = " + std::to_string(sim.t);
        SDL_Color color = {255, 255, 255, 255};
        int w, h;
        GLuint texture = render_text_to_texture(text, color, w, h);
        if (texture) {
            glBindTexture(GL_TEXTURE_2D, texture);
            glBegin(GL_QUADS);
            glTexCoord2f(0, 1); glVertex2f(10, 590 - h);
            glTexCoord2f(1, 1); glVertex2f(10 + w, 590 - h);
            glTexCoord2f(1, 0); glVertex2f(10 + w, 590);
            glTexCoord2f(0, 0); glVertex2f(10, 590);
            glEnd();
            glDeleteTextures(1, &texture);
        }

        glDisable(GL_TEXTURE_2D);
        glDisable(GL_BLEND);
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
}