cmake_minimum_required(VERSION 3.10)
project(PhysicsVisualization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# set(CMAKE_BUILD_TYPE Debug)
# for symbols with gdb --args simulation
# bt full
# info registers
# I do not recommend enabling the symbols.

# Set OpenGL preference to GLVND
set(OpenGL_GL_PREFERENCE GLVND)

# Find SDL2
find_package(SDL2 REQUIRED)
if (NOT SDL2_FOUND)
    message(FATAL_ERROR "SDL2 not found. Please install SDL2 development libraries (e.g., 'sudo apt install libsdl2-dev' on Ubuntu).")
endif()

# Find GLEW
find_package(GLEW REQUIRED)
if (NOT GLEW_FOUND)
    message(FATAL_ERROR "GLEW not found. Please install GLEW development libraries (e.g., 'sudo apt install libglew-dev' on Ubuntu).")
endif()

# Find SDL2_ttf
find_package(SDL2_ttf)
if (NOT SDL2_TTF_FOUND)
    # Fallback: Manually find SDL2_ttf library and headers
    find_library(SDL2_TTF_LIBRARY NAMES SDL2_ttf sdl2_ttf)
    find_path(SDL2_TTF_INCLUDE_DIR NAMES SDL_ttf.h PATH_SUFFIXES SDL2)
    if (SDL2_TTF_LIBRARY AND SDL2_TTF_INCLUDE_DIR)
        set(SDL2_TTF_FOUND TRUE)
        message(STATUS "Manually found SDL2_ttf: ${SDL2_TTF_LIBRARY}")
        message(STATUS "SDL2_ttf include dir: ${SDL2_TTF_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "SDL2_ttf not found. Please install SDL2_ttf development libraries (e.g., 'sudo apt install libsdl2-ttf-dev' on Ubuntu).")
    endif()
endif()

# Find OpenGL and GLU
find_package(OpenGL REQUIRED)
if (NOT OPENGL_FOUND)
    message(FATAL_ERROR "OpenGL not found.")
endif()

find_library(GLU_LIBRARY NAMES GLU glu)
if (NOT GLU_LIBRARY)
    message(FATAL_ERROR "GLU library not found. Please install GLU (e.g., 'sudo apt install libglu1-mesa-dev' on Ubuntu).")
endif()

# Find OpenMP for parallelization
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Print found libraries for debugging
message(STATUS "SDL2 Libraries: ${SDL2_LIBRARIES}")
message(STATUS "SDL2 Include Dirs: ${SDL2_INCLUDE_DIRS}")
message(STATUS "SDL2_ttf Library: ${SDL2_TTF_LIBRARY}")
message(STATUS "SDL2_ttf Include Dir: ${SDL2_TTF_INCLUDE_DIR}")
message(STATUS "OpenGL Libraries: ${OPENGL_LIBRARIES}")
message(STATUS "GLU Library: ${GLU_LIBRARY}")
message(STATUS "GLEW Library: ${GLEW_LIBRARIES}")

# Source files
set(SOURCES
    src/main.cpp
    src/equations.cpp
    src/renderer.cpp
    src/menu.cpp
)

# Create output directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

# Executable
add_executable(simulation ${SOURCES})

# Include directories for the simulation target
target_include_directories(simulation PRIVATE
    ${SDL2_INCLUDE_DIRS}
    ${SDL2_TTF_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)

# Link libraries
target_link_libraries(simulation PRIVATE
    ${SDL2_LIBRARIES}
    ${SDL2_TTF_LIBRARY}
    ${OPENGL_LIBRARIES}
    ${GLU_LIBRARY}
	${GLEW_LIBRARIES}
)

if (OPENMP_FOUND)
    target_link_libraries(simulation PRIVATE OpenMP::OpenMP_CXX)
endif()

# Platform-specific settings
if (ANDROID)
    # Use OpenGL ES for Android
    target_link_libraries(simulation PRIVATE GLESv2)
elseif (UNIX AND NOT APPLE)
    # Linux: Use GLU for desktop OpenGL
    target_link_libraries(simulation PRIVATE ${SDL2_LIBRARIES} ${SDL2_TTF_LIBRARY} ${OPENGL_LIBRARIES} ${GLU_LIBRARY})
elseif (APPLE)
    # macOS: Use OpenGL and GLU frameworks
    target_link_libraries(simulation PRIVATE ${SDL2_LIBRARIES} ${SDL2_TTF_LIBRARY} ${OPENGL_LIBRARIES} ${GLU_LIBRARY})
elseif (WIN32)
    # Windows: Link against OpenGL and GLU
    target_link_libraries(simulation PRIVATE ${SDL2_LIBRARIES} ${SDL2_TTF_LIBRARY} ${OPENGL_LIBRARIES} glu32)
else()
    message(FATAL_ERROR "Unsupported platform. This CMakeLists.txt supports Android, Linux, macOS, and Windows.")
endif()

# Copy equations.txt to build directory
configure_file(${CMAKE_SOURCE_DIR}/equations.txt ${CMAKE_BINARY_DIR}/equations.txt COPYONLY)

# Install target (optional)
install(TARGETS simulation DESTINATION bin)
