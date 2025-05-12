#ifndef RENDERER_H
#define RENDERER_H

#include <SDL2/SDL_ttf.h>
#include <string>
#include "equations.h"

extern float camera_angle;
extern float camera_zoom;
extern bool rotate_camera;
extern float camera_tilt;

enum DisplayMode { POINTS, ISOSURFACE, WIREFRAME, PARTICLES, HYBRID, SURFACE, SPHERE_POINTS };

extern DisplayMode current_display_mode;

void render(const Simulation& sim);
void cleanup_ttf();
uint render_text_to_texture(const std::string& text, SDL_Color color, int& w, int& h);

#endif