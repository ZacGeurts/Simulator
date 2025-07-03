#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include <vulkan/vulkan.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <VulkanMemoryAllocator/include/vk_mem_alloc.h>
#include "equations.h"

enum DisplayMode {
    PARTICLES, POINTS, ISOSURFACE, SURFACE, WIREFRAME, SPHERE_POINTS, HYBRID
};

extern float camera_zoom, camera_angle, camera_tilt;
extern DisplayMode current_display_mode;

struct VulkanContext {
    VkInstance instance;
    VkPhysicalDevice physicalDevice;
    VkDevice device;
    VkQueue graphicsQueue, computeQueue;
    VkSwapchainKHR swapchain;
    std::vector<VkImage> swapchainImages;
    std::vector<VkImageView> swapchainImageViews;
    std::vector<VkFramebuffer> framebuffers;
    VkRenderPass renderPass;
    VkPipelineLayout pipelineLayout;
    VkPipeline graphicsPipeline, computePipeline;
    VkCommandPool commandPool;
    std::vector<VkCommandBuffer> commandBuffers;
    VkDescriptorSetLayout descriptorSetLayout;
    VkDescriptorPool descriptorPool;
    std::vector<VkDescriptorSet> descriptorSets;
    VkBuffer uniformBuffer, scalarBuffer, scalar4DBuffer, vertexBuffer, particleBuffer;
    VmaAllocation uniformAlloc, scalarAlloc, scalar4DAlloc, vertexAlloc, particleAlloc;
    VmaAllocator allocator;
    VkSemaphore imageAvailableSemaphore, renderFinishedSemaphore;
    VkFence inFlightFence;
};

bool init_renderer(SDL_Window* window, VulkanContext& ctx, const Simulation& sim);
void render(const Simulation& sim, VulkanContext& ctx, class Menu& menu);
void cleanup_renderer(VulkanContext& ctx);

#endif