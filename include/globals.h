#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL3/SDL.h>
#include <mutex>
#include <string>
#include "equations.h" // Include for Simulation

// Forward declaration of DisplayMode to avoid including renderer.h
enum class DisplayMode; // Use enum class for stronger typing

// GlobalState struct to encapsulate global variables
struct GlobalState {
    std::mutex sim_mutex;
    Simulation sim;
    DisplayMode current_display_mode;
    float camera_angle = 0.0f;
    float camera_tilt = 0.0f;
    float camera_zoom = 0.3f;
    bool rotate_camera = true;
    bool is_paused = false;
    bool needs_render = true;
    bool is_fullscreen = true;
    SDL_TimerID timer_id = 0;
};

// Single global instance
extern GlobalState global_state;

void log_message(const std::string& message);
Uint32 timer(void* param, SDL_TimerID timerID, Uint32 interval); // Updated for SDL3
void cleanup_ttf();

#endif // GLOBALS_H