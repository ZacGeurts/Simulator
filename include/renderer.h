#ifndef RENDERER_H
#define RENDERER_H

#include "equations.h"

enum DisplayMode {
    POINTS,
    ISOSURFACE,
    WIREFRAME,
    PARTICLES,
    HYBRID
};

extern DisplayMode current_display_mode;
extern float camera_angle;
extern float camera_zoom;

void render(const Simulation& sim);
void cleanup_ttf();

#endif