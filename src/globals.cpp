#include "globals.h"
#include "renderer.h"

GlobalState global_state;

Uint32 timer(void* param, SDL_TimerID timerID, Uint32 interval) {
    {
        std::lock_guard<std::mutex> lock(global_state.sim_mutex);
        if (!global_state.is_paused) {
            global_state.sim.step();
            if (global_state.rotate_camera) {
                global_state.camera_angle += 0.05f;
                global_state.needs_render = true;
            }
        }
    }
    return 16;
}

void log_message(const std::string& message) {
    static std::set<std::string> logged_messages;
    std::ofstream log("log.txt", std::ios::app);
    if (log.is_open() && logged_messages.find(message) == logged_messages.end()) {
        std::time_t now = std::time(nullptr);
        log << "[" << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "] " << message << "\n";
        logged_messages.insert(message);
        log.close();
    }
}