#include "menu.h"
#include "globals.h"
#include <SDL2/SDL_opengl.h>
#include <iostream>

Menu::Menu() : hover_index_(-1), tooltip_texture_(0), tooltip_w_(0), tooltip_h_(0), window_(nullptr) {
}

Menu::~Menu() {
    for (auto& button : buttons_) {
        if (button.texture) {
            glDeleteTextures(1, &button.texture);
        }
    }
    if (tooltip_texture_) {
        glDeleteTextures(1, &tooltip_texture_);
    }
}

bool Menu::init(SDL_Window* window, Simulation& sim) {
    window_ = window;

    // Define buttons with scancode, label, tooltip, and action
    buttons_ = {
        {SDL_SCANCODE_1, "1", "Points mode", {}, [](Simulation& s) { current_display_mode = POINTS; needs_render = true; }},
        {SDL_SCANCODE_2, "2", "Isosurface mode", {}, [](Simulation& s) { current_display_mode = ISOSURFACE; needs_render = true; }},
        {SDL_SCANCODE_3, "3", "Wireframe mode", {}, [](Simulation& s) { current_display_mode = WIREFRAME; needs_render = true; }},
        {SDL_SCANCODE_4, "4", "Particles mode", {}, [](Simulation& s) { current_display_mode = PARTICLES; needs_render = true; }},
        {SDL_SCANCODE_5, "5", "Hybrid mode", {}, [](Simulation& s) { current_display_mode = HYBRID; needs_render = true; }},
        {SDL_SCANCODE_6, "6", "Surface mode", {}, [](Simulation& s) { current_display_mode = SURFACE; needs_render = true; }},
        {SDL_SCANCODE_7, "7", "Sphere mode", {}, [](Simulation& s) { current_display_mode = SPHERE_POINTS; needs_render = true; }},
        {SDL_SCANCODE_0, "0", "Center camera", {}, [](Simulation& s) {
            camera_angle = 0.0f; camera_tilt = 0.0f; camera_zoom = 0.3f; needs_render = true;
        }},
        {SDL_SCANCODE_UP, "Up", "Zoom in", {}, [](Simulation& s) {
            camera_zoom *= 0.95f; camera_zoom = std::max(camera_zoom, 0.001f); needs_render = true;
        }},
        {SDL_SCANCODE_LEFT, "Left", "Previous equation", {}, [](Simulation& s) { s.prev_equation(); needs_render = true; }},
        {SDL_SCANCODE_RIGHT, "Right", "Next equation", {}, [](Simulation& s) { s.next_equation(); needs_render = true; }},
        {SDL_SCANCODE_DOWN, "Down", "Zoom out", {}, [](Simulation& s) {
            camera_zoom *= 1.05f; camera_zoom = std::min(camera_zoom, 5.0f); needs_render = true;
        }},
        {SDL_SCANCODE_SPACE, "Space", "Toggle pause", {}, [](Simulation& s) {
            is_paused = !is_paused; rotate_camera = !is_paused; needs_render = true;
        }},
        {SDL_SCANCODE_RETURN, "Enter", "Step simulation (if paused)", {}, [](Simulation& s) {
            if (is_paused) { s.step(); needs_render = true; }
        }},
        {SDL_SCANCODE_C, "C", "Pan camera left (if paused)", {}, [](Simulation& s) {
            if (is_paused) { camera_angle -= 0.01f; needs_render = true; }
        }},
        {SDL_SCANCODE_V, "V", "Pan camera right (if paused)", {}, [](Simulation& s) {
            if (is_paused) { camera_angle += 0.01f; needs_render = true; }
        }},
        {SDL_SCANCODE_B, "B", "Tilt camera up (if paused)", {}, [](Simulation& s) {
            if (is_paused) { camera_tilt += 0.01f; needs_render = true; }
        }},
        {SDL_SCANCODE_G, "G", "Tilt camera down (if paused)", {}, [](Simulation& s) {
            if (is_paused) { camera_tilt -= 0.01f; needs_render = true; }
        }},
        {SDL_SCANCODE_A, "A", "Decrease vis scale by 10%", {}, [](Simulation& s) {
            if (!s.equations.empty()) {
                if (auto* eq = dynamic_cast<CustomEquation*>(s.equations[s.current_equation])) {
                    double current_scale = eq->getVisScale();
                    current_scale /= 1.1;
                    eq->setVisScale(current_scale);
                    s.compute();
                    needs_render = true;
                }
            }
        }},
        {SDL_SCANCODE_S, "S", "Increase vis scale by 10%", {}, [](Simulation& s) {
            if (!s.equations.empty()) {
                if (auto* eq = dynamic_cast<CustomEquation*>(s.equations[s.current_equation])) {
                    double current_scale = eq->getVisScale();
                    current_scale *= 1.1;
                    eq->setVisScale(current_scale);
                    s.compute();
                    needs_render = true;
                }
            }
        }},
        {SDL_SCANCODE_Q, "Q", "Decrease color intensity by 10%", {}, [](Simulation& s) {
            if (!s.equations.empty()) {
                if (auto* eq = dynamic_cast<CustomEquation*>(s.equations[s.current_equation])) {
                    double current_intensity = eq->getVisColorIntensity();
                    current_intensity /= 1.1;
                    eq->setVisColorIntensity(current_intensity);
                    s.compute();
                    needs_render = true;
                }
            }
        }},
        {SDL_SCANCODE_W, "W", "Increase color intensity by 10%", {}, [](Simulation& s) {
            if (!s.equations.empty()) {
                if (auto* eq = dynamic_cast<CustomEquation*>(s.equations[s.current_equation])) {
                    double current_intensity = eq->getVisColorIntensity();
                    current_intensity *= 1.1;
                    eq->setVisColorIntensity(current_intensity);
                    s.compute();
                    needs_render = true;
                }
            }
        }},
        {SDL_SCANCODE_X, "X", "Increase grid size", {}, [](Simulation& s) {
            s.size = std::min(s.size + 1, 50);
            s.ricci.resize(s.size, std::vector<std::vector<Tensor>>(s.size, std::vector<Tensor>(s.size)));
            s.divergence.resize(s.size, std::vector<std::vector<Tensor>>(s.size, std::vector<Tensor>(s.size)));
            s.scalar.resize(s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0)));
            s.computed_scalar.resize(s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0)));
            s.scalar_4d.resize(s.size, std::vector<std::vector<std::vector<double>>>(
                s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0))));
            s.setup();
            needs_render = true;
        }},
        {SDL_SCANCODE_Z, "Z", "Decrease grid size", {}, [](Simulation& s) {
            s.size = std::max(s.size - 1, 5);
            s.ricci.resize(s.size, std::vector<std::vector<Tensor>>(s.size, std::vector<Tensor>(s.size)));
            s.divergence.resize(s.size, std::vector<std::vector<Tensor>>(s.size, std::vector<Tensor>(s.size)));
            s.scalar.resize(s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0)));
            s.computed_scalar.resize(s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0)));
            s.scalar_4d.resize(s.size, std::vector<std::vector<std::vector<double>>>(
                s.size, std::vector<std::vector<double>>(s.size, std::vector<double>(s.size, 0.0))));
            s.setup();
            needs_render = true;
        }},
        {SDL_SCANCODE_R, "R", "Reset simulation", {}, [](Simulation& s) {
            camera_angle = 0.0f; camera_tilt = 0.0f; camera_zoom = 0.3f; s.t = 0.0; s.setup(); needs_render = true;
        }},
        {SDL_SCANCODE_P, "P", "Reset particles", {}, [](Simulation& s) { s.particles.clear(); s.setup(); needs_render = true; }},
        {SDL_SCANCODE_F, "F", "Toggle fullscreen", {}, [](Simulation& s) {
            is_fullscreen = !is_fullscreen;
            const int DEFAULT_WIDTH = 800;
            const int DEFAULT_HEIGHT = 600;
            const float TARGET_ASPECT = 4.0f / 3.0f;

            SDL_Window* window = SDL_GL_GetCurrentWindow();
            if (is_fullscreen) {
                if (SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN_DESKTOP) < 0) {
                    log_message("Failed to switch to fullscreen: " + std::string(SDL_GetError()));
                    is_fullscreen = false;
                    return;
                }
                int drawable_w, drawable_h, window_w, window_h;
                SDL_GL_GetDrawableSize(window, &drawable_w, &drawable_h);
                SDL_GetWindowSize(window, &window_w, &window_h);
                float scale_x = (float)drawable_w / window_w;
                float scale_y = (float)drawable_h / window_h;
                float window_aspect = (float)window_w / window_h;
                int vp_width, vp_height, vp_x, vp_y;
                if (window_aspect > TARGET_ASPECT) {
                    vp_height = window_h;
                    vp_width = (int)(vp_height * TARGET_ASPECT);
                    vp_x = (window_w - vp_width) / 2;
                    vp_y = 0;
                } else {
                    vp_width = window_w;
                    vp_height = (int)(vp_width / TARGET_ASPECT);
                    vp_x = 0;
                    vp_y = (window_h - vp_height) / 2;
                }
                glViewport(vp_x * scale_x, vp_y * scale_y, vp_width * scale_x, vp_height * scale_y);
            } else {
                if (SDL_SetWindowFullscreen(window, 0) < 0) {
                    log_message("Failed to switch to windowed mode: " + std::string(SDL_GetError()));
                    is_fullscreen = true;
                    return;
                }
                SDL_SetWindowSize(window, DEFAULT_WIDTH, DEFAULT_HEIGHT);
                SDL_SetWindowPosition(window, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED);
                glViewport(0, 0, DEFAULT_WIDTH, DEFAULT_HEIGHT);
            }
            needs_render = true;
        }},
        {SDL_SCANCODE_ESCAPE, "Esc", "Quit application", {}, [](Simulation& s) {
            SDL_Event quit_event;
            quit_event.type = SDL_QUIT;
            SDL_PushEvent(&quit_event);
        }}
    };

    // Set up menu box and button positions
    const int button_w = 50, button_h = 30, padding = 5;
    const int menu_x = 10, menu_y = 200; // Top-left, 200 pixels down
    int current_y = menu_y + padding;

    // Row 1: 1, 2, 3, 4, 5, 6, 7, 0
    for (size_t i = 0; i < 8; ++i) {
        int x = menu_x + padding + static_cast<int>(i) * (button_w + padding); // Avoid narrowing
        buttons_[i].rect = {x, current_y, button_w, button_h};
    }
    current_y += button_h + padding;

    // Arrow group: Up, Left, Right, Down
    // Up centered above Left and Right
    buttons_[8].rect = {menu_x + padding + button_w / 2 + padding / 2, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Left and Right side-by-side
    buttons_[9].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[10].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Down centered below
    buttons_[11].rect = {menu_x + padding + button_w / 2 + padding / 2, current_y, button_w, button_h};
    current_y += button_h + padding;

    // Row 3: Space, Enter
    buttons_[12].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[13].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;

    // Function group
    // Row 1: C, V (pan left/right)
    buttons_[14].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[15].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Row 2: B, G (tilt up/down)
    buttons_[16].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[17].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Row 3: A, S (vis scale)
    buttons_[18].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[19].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Row 4: Q, W (color intensity)
    buttons_[20].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[21].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Row 5: X, Z (grid size)
    buttons_[22].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[23].rect = {menu_x + padding + button_w + padding, current_y, button_w, button_h};
    current_y += button_h + padding;
    // Row 6: R, P
    buttons_[24].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[25].rect = {menu_x + padding + (button_w + padding), current_y, button_w, button_h};
	current_y += button_h + padding;
	// Row 7: F, Esc
    buttons_[26].rect = {menu_x + padding, current_y, button_w, button_h};
    buttons_[27].rect = {menu_x + padding + (button_w + padding), current_y, button_w, button_h};
    current_y += button_h + padding;

    // Set menu rectangle to encompass all buttons
    menu_rect_ = {menu_x, menu_y, 8 * button_w + 9 * padding, current_y - menu_y};

    // Create textures for button labels
    return create_button_textures();
}

bool Menu::create_button_textures() {
    SDL_Color white = {255, 255, 255, 255};
    for (auto& button : buttons_) {
        button.texture = render_text_to_texture(button.label, white, button.texture_w, button.texture_h);
        if (!button.texture) {
            std::cerr << "Failed to create texture for button: " << button.label << std::endl;
            return false;
        }
    }
    return true;
}

void Menu::update_tooltip(int button_index) {
    if (tooltip_texture_) {
        glDeleteTextures(1, &tooltip_texture_);
        tooltip_texture_ = 0;
    }
    if (button_index >= 0 && button_index < static_cast<int>(buttons_.size())) {
        SDL_Color yellow = {255, 255, 0, 255}; // Bright yellow
        tooltip_texture_ = render_text_to_texture(buttons_[button_index].tooltip, yellow, tooltip_w_, tooltip_h_);
    }
}

bool Menu::point_in_rect(int x, int y, const SDL_Rect& rect, int vp_x, int vp_y) const {
    // Adjust mouse coordinates to viewport
    int adjusted_x = x - vp_x;
    int adjusted_y = y - vp_y;
    return adjusted_x >= rect.x && adjusted_x < rect.x + rect.w &&
           adjusted_y >= rect.y && adjusted_y < rect.y + rect.h;
}

void Menu::handle_event(const SDL_Event& e, Simulation& sim) {
    // Get window and drawable sizes
    int window_w, window_h, drawable_w, drawable_h;
    SDL_GetWindowSize(window_, &window_w, &window_h);
    SDL_GL_GetDrawableSize(window_, &drawable_w, &drawable_h);
    float scale_x = (float)drawable_w / window_w;
    float scale_y = (float)drawable_h / window_h;

    // Compute viewport for 4:3 aspect ratio
    const float TARGET_ASPECT = 4.0f / 3.0f;
    float window_aspect = (float)window_w / window_h;
    int vp_width, vp_height, vp_x, vp_y;
    if (window_aspect > TARGET_ASPECT) {
        vp_height = window_h;
        vp_width = (int)(vp_height * TARGET_ASPECT);
        vp_x = (window_w - vp_width) / 2;
        vp_y = 0;
    } else {
        vp_width = window_w;
        vp_height = (int)(vp_width / TARGET_ASPECT);
        vp_x = 0;
        vp_y = (window_h - vp_height) / 2;
    }

    if (e.type == SDL_MOUSEMOTION) {
        // Scale mouse coordinates to drawable size
        int x = (int)(e.motion.x * scale_x);
        int y = (int)(e.motion.y * scale_y);
        int new_hover = -1;
        for (size_t i = 0; i < buttons_.size(); ++i) {
            if (point_in_rect(x, y, buttons_[i].rect, vp_x * scale_x, vp_y * scale_y)) {
                new_hover = i;
                break;
            }
        }
        if (new_hover != hover_index_) {
            hover_index_ = new_hover;
            update_tooltip(hover_index_);
        }
    } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
        // Scale mouse coordinates to drawable size
        int x = (int)(e.button.x * scale_x);
        int y = (int)(e.button.y * scale_y);
        for (size_t i = 0; i < buttons_.size(); ++i) {
            if (point_in_rect(x, y, buttons_[i].rect, vp_x * scale_x, vp_y * scale_y)) {
                std::lock_guard<std::mutex> lock(sim_mutex);
                buttons_[i].action(sim);
                break;
            }
        }
    }
}

void Menu::render() const {
    // Get window and drawable sizes
    int window_w, window_h, drawable_w, drawable_h;
    SDL_GetWindowSize(window_, &window_w, &window_h);
    SDL_GL_GetDrawableSize(window_, &drawable_w, &drawable_h);
    float scale_x = (float)drawable_w / window_w;
    float scale_y = (float)drawable_h / window_h;

    // Compute viewport for 4:3 aspect ratio
    const float TARGET_ASPECT = 4.0f / 3.0f;
    float window_aspect = (float)window_w / window_h;
    int vp_width, vp_height, vp_x, vp_y;
    if (window_aspect > TARGET_ASPECT) {
        vp_height = window_h;
        vp_width = (int)(vp_height * TARGET_ASPECT);
        vp_x = (window_w - vp_width) / 2;
        vp_y = 0;
    } else {
        vp_width = window_w;
        vp_height = (int)(vp_width / TARGET_ASPECT);
        vp_x = 0;
        vp_y = (window_h - vp_height) / 2;
    }

    // Set up orthographic projection for 2D rendering in logical coordinates
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, window_w, window_h, 0, -1, 1); // Logical coordinates, Y flipped
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_TEXTURE_2D);

    // Draw semi-transparent menu background
    glColor4f(0.1f, 0.1f, 0.1f, 0.7f);
    glBegin(GL_QUADS);
    glVertex2i(menu_rect_.x, menu_rect_.y);
    glVertex2i(menu_rect_.x + menu_rect_.w, menu_rect_.y);
    glVertex2i(menu_rect_.x + menu_rect_.w, menu_rect_.y + menu_rect_.h);
    glVertex2i(menu_rect_.x, menu_rect_.y + menu_rect_.h);
    glEnd();

    // Draw buttons
    for (size_t i = 0; i < buttons_.size(); ++i) {
        const auto& button = buttons_[i];
        // Highlight if hovered
        glColor4f(i == static_cast<size_t>(hover_index_) ? 0.5f : 0.3f, 0.3f, 0.3f, 0.8f);
        glBegin(GL_QUADS);
        glVertex2i(button.rect.x, button.rect.y);
        glVertex2i(button.rect.x + button.rect.w, button.rect.y);
        glVertex2i(button.rect.x + button.rect.w, button.rect.y + button.rect.h);
        glVertex2i(button.rect.x, button.rect.y + button.rect.h);
        glEnd();

        // Draw button label
        glBindTexture(GL_TEXTURE_2D, button.texture);
        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2i(button.rect.x, button.rect.y);
        glTexCoord2f(1, 0); glVertex2i(button.rect.x + button.rect.w, button.rect.y);
        glTexCoord2f(1, 1); glVertex2i(button.rect.x + button.rect.w, button.rect.y + button.rect.h);
        glTexCoord2f(0, 1); glVertex2i(button.rect.x, button.rect.y + button.rect.h);
        glEnd();
    }

    // Draw tooltip if hovering
    if (hover_index_ >= 0 && tooltip_texture_) {
        int mouse_x, mouse_y;
        SDL_GetMouseState(&mouse_x, &mouse_y);
        // Scale mouse position to logical coordinates
        mouse_x = (int)(mouse_x * scale_x / scale_x); // Effectively mouse_x, but kept for clarity
        mouse_y = (int)(mouse_y * scale_y / scale_y); // Effectively mouse_y
        glBindTexture(GL_TEXTURE_2D, tooltip_texture_);
        glColor4f(0.0f, 0.0f, 0.0f, 0.8f); // Black background for tooltip
        glBegin(GL_QUADS);
        glVertex2i(mouse_x, mouse_y);
        glVertex2i(mouse_x + tooltip_w_, mouse_y);
        glVertex2i(mouse_x + tooltip_w_, mouse_y + tooltip_h_);
        glVertex2i(mouse_x, mouse_y + tooltip_h_);
        glEnd();
        glColor4f(1.0f, 1.0f, 0.0f, 1.0f); // Yellow text
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2i(mouse_x, mouse_y);
        glTexCoord2f(1, 0); glVertex2i(mouse_x + tooltip_w_, mouse_y);
        glTexCoord2f(1, 1); glVertex2i(mouse_x + tooltip_w_, mouse_y + tooltip_h_);
        glTexCoord2f(0, 1); glVertex2i(mouse_x, mouse_y + tooltip_h_);
        glEnd();
    }

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    // Restore projection and modelview matrices
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}