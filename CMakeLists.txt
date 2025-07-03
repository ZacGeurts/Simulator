cmake_minimum_required(VERSION 3.10)
project(PhysicsVisualization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find SDL3
find_package(SDL3 REQUIRED)
if (NOT SDL3_FOUND)
    message(FATAL_ERROR "SDL3 not found. Install SDL3 development libraries (e.g., 'sudo apt install libsdl3-dev' on Ubuntu or build from https://github.com/libsdl-org/SDL).")
endif()

# Find Vulkan
find_package(Vulkan REQUIRED COMPONENTS glslc glslangValidator)
if (NOT Vulkan_FOUND)
    message(FATAL_ERROR "Vulkan not found. Install Vulkan SDK (e.g., 'sudo apt install vulkan-sdk' on Ubuntu or download from https://vulkan.lunarg.com/).")
endif()

# Find GLM
find_package(glm REQUIRED)
if (NOT glm_FOUND)
    message(FATAL_ERROR "GLM not found. Install GLM development libraries (e.g., 'sudo apt install libglm-dev' on Ubuntu).")
endif()

# Check for glslangValidator
find_program(GLSLANG_VALIDATOR glslangValidator HINTS ${Vulkan_GLSLANG_VALIDATOR_EXECUTABLE})
if (NOT GLSLANG_VALIDATOR)
    message(FATAL_ERROR "glslangValidator not found. Install glslang-tools (e.g., 'sudo apt install glslang-tools' on Ubuntu).")
endif()

# ImGui: Add as source files (assumes imgui/ directory in project root)
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/imgui)
add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/backends/imgui_impl_sdl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_vulkan.cpp
)
target_include_directories(imgui PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)
target_link_libraries(imgui PUBLIC SDL3::SDL3 Vulkan::Vulkan)

# Vulkan Memory Allocator (header-only)
set(VMA_DIR ${CMAKE_SOURCE_DIR}/VulkanMemoryAllocator)
# Check if VMA header exists
if (NOT EXISTS ${VMA_DIR}/include/vk_mem_alloc.h)
    message(FATAL_ERROR "VulkanMemoryAllocator header not found at ${VMA_DIR}/include/vk_mem_alloc.h. "
                        "Please clone VMA (e.g., 'git clone https://github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator.git' in ${CMAKE_SOURCE_DIR}).")
endif()
add_library(vma INTERFACE)
target_include_directories(vma INTERFACE ${VMA_DIR}/include)
target_link_libraries(vma INTERFACE Vulkan::Vulkan)

# Source files
set(SOURCES
    src/main.cpp
    src/renderer.cpp
    src/menu.cpp
    src/equations.cpp
)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Executable
add_executable(simulation ${SOURCES})
target_include_directories(simulation PRIVATE
    ${SDL3_INCLUDE_DIRS}
    ${Vulkan_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(simulation PRIVATE
    SDL3::SDL3
    Vulkan::Vulkan
    glm::glm
    imgui
    vma
)

# Shader compilation (assumes shaders in shaders/ directory)
set(SHADER_DIR ${CMAKE_SOURCE_DIR}/shaders)
file(GLOB SHADERS "${SHADER_DIR}/*.vert" "${SHADER_DIR}/*.frag" "${SHADER_DIR}/*.comp" "${SHADER_DIR}/*.glsl")
foreach(SHADER ${SHADERS})
    get_filename_component(SHADER_NAME ${SHADER} NAME)
    add_custom_command(
        OUTPUT ${CMAKE_BINARY_DIR}/${SHADER_NAME}.spv
        COMMAND ${GLSLANG_VALIDATOR} -V ${SHADER} -o ${CMAKE_BINARY_DIR}/${SHADER_NAME}.spv
        DEPENDS ${SHADER}
        COMMENT "Compiling ${SHADER_NAME} to SPIR-V"
    )
    list(APPEND SHADER_SPVS ${CMAKE_BINARY_DIR}/${SHADER_NAME}.spv)
endforeach()

# Custom target for shaders
add_custom_target(shaders DEPENDS ${SHADER_SPVS})
add_dependencies(simulation shaders)

# Copy SPIR-V shaders and binary to parent directory (next to equations.txt)
add_custom_command(TARGET simulation POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/*.spv ${CMAKE_SOURCE_DIR}/../
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:simulation> ${CMAKE_SOURCE_DIR}/../simulation
    COMMENT "Copying simulation binary and SPIR-V shaders to ${CMAKE_SOURCE_DIR}/../"
)

# Platform-specific settings
if (UNIX AND NOT APPLE)
    message(STATUS "Configuring for Linux")
    target_link_libraries(simulation PRIVATE dl)
elseif (APPLE)
    message(WARNING "macOS support requires SDL3 and Vulkan setup. Ensure MoltenVK is installed.")
    target_link_libraries(simulation PRIVATE "-framework Cocoa")
elseif (WIN32)
    message(STATUS "Configuring for Windows")
    target_link_libraries(simulation PRIVATE gdi32)
else()
    message(FATAL_ERROR "Unsupported platform. This CMakeLists.txt supports Linux, macOS, and Windows.")
endif()