#ifndef RENDERER_H
#define RENDERER_H

#include "equations.h"

extern float camera_angle;
extern float camera_zoom;
extern bool rotate_camera;
extern float camera_tilt;

enum DisplayMode { POINTS, ISOSURFACE, WIREFRAME, PARTICLES, HYBRID, SURFACE };

extern DisplayMode current_display_mode;

void render(const Simulation& sim);
void cleanup_ttf();

#endif