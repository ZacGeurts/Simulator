#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL3/SDL.h>
#include <mutex>
#include <string>
#include "equations.h" // Include for Simulation

// GlobalState struct to encapsulate global variables
struct GlobalState {
    std::mutex sim_mutex;
    Simulation sim;
    enum class DisplayMode { // Define enum class with all required members
        POINTS,
        ISOSURFACE,
        WIREFRAME,
        PARTICLES,
        HYBRID,
        SURFACE,
        SPHERE_POINTS
    } current_display_mode = DisplayMode::POINTS; // Default to POINTS
    float camera_angle = 0.0f;
    float camera_tilt = 0.0f;
    float camera_zoom = 0.3f;
    bool rotate_camera = true;
    bool is_paused = false;
    bool needs_render = true;
    bool is_fullscreen = false; // Changed to false for 4:3 windowed start
    SDL_TimerID timer_id = 0;
};

// Single global instance
extern GlobalState global_state;

void log_message(const std::string& message);
Uint32 timer(void* param, SDL_TimerID timerID, Uint32 interval); // Updated for SDL3
void cleanup_ttf();

#endif // GLOBALS_H