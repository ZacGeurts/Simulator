#ifndef MENU_H
#define MENU_H

#include "renderer.h"
#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <string>
#include <vector>
#include <functional>

class Simulation;

struct MenuButton {
    SDL_Scancode scancode;
    std::string label;
    std::string tooltip;
    SDL_Rect rect;
    std::function<void(Simulation&)> action;
    GLuint texture;
    int texture_w, texture_h;
};

class Menu {
public:
    Menu();
    ~Menu();
    bool init(SDL_Window* window, Simulation& sim);
    void render() const;
    void handle_event(const SDL_Event& e, Simulation& sim);

private:
    std::vector<MenuButton> buttons_;
    SDL_Rect menu_rect_;
    int hover_index_;
    GLuint tooltip_texture_;
    int tooltip_w_, tooltip_h_;
    SDL_Window* window_;
    bool create_button_textures();
    void update_tooltip(int button_index);
    bool point_in_rect(int x, int y, const SDL_Rect& rect, int vp_x, int vp_y) const;
};

#endif // MENU_H