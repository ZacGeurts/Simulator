cmake_minimum_required(VERSION 3.20)
project(Simulation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find dependencies
find_package(Vulkan REQUIRED)
find_package(SDL3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(glm REQUIRED)

# ImGui
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/imgui)
add_library(imgui STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/backends/imgui_impl_sdl3.cpp
    ${IMGUI_DIR}/backends/imgui_impl_vulkan.cpp
)
target_include_directories(imgui PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends ${Vulkan_INCLUDE_DIRS} ${SDL3_INCLUDE_DIRS})

# Vulkan Memory Allocator (Header-Only)
set(VMA_DIR ${CMAKE_SOURCE_DIR}/include)
add_library(vma INTERFACE)
target_include_directories(vma INTERFACE ${VMA_DIR} ${Vulkan_INCLUDE_DIRS})
target_compile_definitions(vma INTERFACE VMA_STATIC_VULKAN_FUNCTIONS=0 VMA_DYNAMIC_VULKAN_FUNCTIONS=1)

# Source files
set(SOURCES
    src/main.cpp
    src/renderer.cpp
    src/menu.cpp
    src/equations.cpp
    src/globals.cpp
)

# Shader compilation
file(GLOB SHADERS shaders/*.vert shaders/*.frag)
foreach(SHADER ${SHADERS})
    get_filename_component(SHADER_NAME ${SHADER} NAME)
    set(SHADER_OUTPUT ${CMAKE_BINARY_DIR}/shaders/${SHADER_NAME}.spv)
    add_custom_command(
        OUTPUT ${SHADER_OUTPUT}
        COMMAND ${Vulkan_GLSLC_EXECUTABLE} ${SHADER} -o ${SHADER_OUTPUT}
        DEPENDS ${SHADER}
        COMMENT "Compiling ${SHADER_NAME}"
    )
    list(APPEND SHADER_OUTPUTS ${SHADER_OUTPUT})
endforeach()
add_custom_target(shaders DEPENDS ${SHADER_OUTPUTS})

# Executable
add_executable(${PROJECT_NAME} ${SOURCES})
add_dependencies(${PROJECT_NAME} shaders)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${Vulkan_INCLUDE_DIRS}
    ${SDL3_INCLUDE_DIRS}
    ${GLM_INCLUDE_DIRS}
)
target_link_libraries(${PROJECT_NAME} PRIVATE
    Vulkan::Vulkan
    SDL3::SDL3
    OpenMP::OpenMP_CXX
    glm::glm
    imgui
    vma
)