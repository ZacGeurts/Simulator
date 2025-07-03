#ifndef MENU_H
#define MENU_H

#include <SDL3/SDL.h>
#include "equations.h"

class Menu {
public:
    Menu();
    ~Menu();
    bool init(SDL_Window* window, Simulation& sim);
    void render(Simulation& sim);
    void handle_event(const SDL_Event& e, Simulation& sim);

private:
    SDL_Window* window_;
};

#endif // MENU_H