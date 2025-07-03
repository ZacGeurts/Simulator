#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL3/SDL.h>
#include <mutex>
#include "renderer.h"

extern std::mutex sim_mutex;
extern Simulation sim;
extern DisplayMode current_display_mode;
extern float camera_angle;
extern float camera_tilt;
extern float camera_zoom;
extern bool rotate_camera;
extern bool is_paused;
extern bool needs_render;
extern bool is_fullscreen;
extern SDL_TimerID timer_id;

void log_message(const std::string& message);
Uint32 timer(Uint32 interval, void* param);
void cleanup_ttf();

#endif // GLOBALS_H