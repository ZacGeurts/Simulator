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
find_package(Vulkan REQUIRED COMPONENTS glslc)
if (NOT Vulkan_FOUND)
    message(FATAL_ERROR "Vulkan not found. Install Vulkan SDK (e.g., 'sudo apt install vulkan-sdk' on Ubuntu or download from https://vulkan.lunarg.com/).")
endif()

# Find GLM
find_package(glm REQUIRED)
if (NOT glm_FOUND)
    message(FATAL_ERROR "GLM not found. Install GLM development libraries (e.g., 'sudo apt install libglm-dev' on Ubuntu).")
endif()

# Check for glslc
find_program(GLSLC glslc HINTS ${Vulkan_GLSLC_EXECUTABLE})
if (NOT GLSLC)
    message(FATAL_ERROR "glslc not found. Install glslc (e.g., 'sudo apt install glslc' on Ubuntu).")
endif()

# ImGui: Add as source files
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

# Vulkan Memory Allocator
set(VMA_DIR ${CMAKE_SOURCE_DIR}/VulkanMemoryAllocator)
if (NOT EXISTS ${VMA_DIR}/include/vk_mem_alloc.h)
    message(FATAL_ERROR "VulkanMemoryAllocator header not found at ${VMA_DIR}/include/vk_mem_alloc.h. "
                        "Please clone VMA (e.g., 'git clone https://github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator.git' in ${CMAKE_SOURCE_DIR}).")
endif()

# Check if vk_mem_alloc.cpp exists for non-header-only VMA
set(VMA_SRC ${VMA_DIR}/src/vk_mem_alloc.cpp)
if (EXISTS ${VMA_SRC})
    # Build VMA as a static library if vk_mem_alloc.cpp exists
    add_library(vma STATIC ${VMA_SRC})
    target_include_directories(vma PUBLIC ${VMA_DIR}/include)
    target_link_libraries(vma PUBLIC Vulkan::Vulkan)
else()
    # Use header-only VMA
    add_library(vma INTERFACE)
    target_include_directories(vma INTERFACE ${VMA_DIR}/include)
    target_link_libraries(vma INTERFACE Vulkan::Vulkan)
    target_compile_definitions(vma INTERFACE VMA_STATIC_VULKAN_FUNCTIONS=0 VMA_DYNAMIC_VULKAN_FUNCTIONS=1)
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/renderer.cpp
    src/menu.cpp
    src/equations.cpp
    src/globals.cpp
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
# Check if shader directory exists
if (NOT EXISTS ${SHADER_DIR})
    message(WARNING "Shader directory ${SHADER_DIR} not found. Skipping shader compilation.")
else()
    file(GLOB SHADERS "${SHADER_DIR}/*.vert" "${SHADER_DIR}/*.frag" "${SHADER_DIR}/*.comp" "${SHADER_DIR}/frag.glsl" "${SHADER_DIR}/vert.glsl")
    if (NOT SHADERS)
        message(WARNING "No shader files found in ${SHADER_DIR}. Skipping shader compilation.")
    else()
        set(SHADER_SPVS)
        foreach(SHADER ${SHADERS})
            get_filename_component(SHADER_NAME ${SHADER} NAME)
            get_filename_component(SHADER_EXT ${SHADER} EXT)
            
            # Skip include files like tables.glsl
            if (SHADER_NAME MATCHES "tables.glsl")
                continue()
            endif()

            set(SHADER_OUTPUT ${CMAKE_BINARY_DIR}/${SHADER_NAME}.spv)
            # Determine dependencies (e.g., tables.glsl for compute.comp)
            set(SHADER_DEPS ${SHADER})
            if (SHADER_NAME MATCHES "compute.comp")
                if (EXISTS ${SHADER_DIR}/tables.glsl)
                    list(APPEND SHADER_DEPS ${SHADER_DIR}/tables.glsl)
                else()
                    message(WARNING "Dependency tables.glsl not found for compute.comp")
                endif()
            endif()

            # Determine shader stage for .glsl files
            set(SHADER_STAGE_FLAG "")
            if (SHADER_NAME MATCHES "frag.glsl")
                set(SHADER_STAGE_FLAG "-fshader-stage=frag")
            elseif (SHADER_NAME MATCHES "vert.glsl")
                set(SHADER_STAGE_FLAG "-fshader-stage=vert")
            endif()

            add_custom_command(
                OUTPUT ${SHADER_OUTPUT}
                COMMAND ${GLSLC} --target-env=vulkan1.2 ${SHADER_STAGE_FLAG} ${SHADER} -o ${SHADER_OUTPUT}
                DEPENDS ${SHADER_DEPS}
                COMMENT "Compiling ${SHADER_NAME} to SPIR-V"
                VERBATIM
            )
            list(APPEND SHADER_SPVS ${SHADER_OUTPUT})
        endforeach()

        # Custom target for shaders
        add_custom_target(shaders DEPENDS ${SHADER_SPVS})
        add_dependencies(simulation shaders)
    endif()
endif()

# Copy SPIR-V shaders and binary to project root (if shaders exist)
if (SHADER_SPVS)
    add_custom_command(TARGET simulation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/shaders
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SHADER_SPVS} ${CMAKE_SOURCE_DIR}/shaders/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:simulation> ${CMAKE_SOURCE_DIR}/simulation
        COMMENT "Copying simulation binary and SPIR-V shaders to ${CMAKE_SOURCE_DIR}/shaders/"
        VERBATIM
    )
else()
    add_custom_command(TARGET simulation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:simulation> ${CMAKE_SOURCE_DIR}/simulation
        COMMENT "Copying simulation binary to ${CMAKE_SOURCE_DIR}/ (no shaders to copy)"
        VERBATIM
    )
endif()

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