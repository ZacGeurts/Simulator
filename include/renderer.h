// Copyright (c) 2025, Zachary Geurts

#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include <SDL3/SDL_vulkan.h>
#include <vulkan/vulkan.h>
#include <vk_mem_alloc.h>
#include <vector>

// Forward declarations
class Simulation;
class Menu;

struct VulkanContext {
    VkInstance instance = VK_NULL_HANDLE;
    VkPhysicalDevice physicalDevice = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    VkQueue graphicsQueue = VK_NULL_HANDLE;
    VkQueue computeQueue = VK_NULL_HANDLE;
    VkQueue presentQueue = VK_NULL_HANDLE;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VmaAllocator allocator = VK_NULL_HANDLE;
    std::vector<VkImage> swapchainImages;
    VkFormat swapchainFormat;
    VkExtent2D swapchainExtent;
    VkSwapchainKHR swapchain = VK_NULL_HANDLE;
    VkRenderPass renderPass = VK_NULL_HANDLE;
    VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
    VkPipeline graphicsPipeline = VK_NULL_HANDLE;
    VkPipeline computePipeline = VK_NULL_HANDLE;
    VkDescriptorSetLayout descriptorSetLayout = VK_NULL_HANDLE;
    VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
    std::vector<VkDescriptorSet> descriptorSets;
    VkCommandPool commandPool = VK_NULL_HANDLE;
    std::vector<VkCommandBuffer> commandBuffers;
    VkSemaphore imageAvailableSemaphore = VK_NULL_HANDLE;
    VkSemaphore renderFinishedSemaphore = VK_NULL_HANDLE;
    VkFence inFlightFence = VK_NULL_HANDLE;
    VkBuffer uniformBuffer = VK_NULL_HANDLE;
    VmaAllocation uniformAlloc = VK_NULL_HANDLE;
    VkBuffer scalarBuffer = VK_NULL_HANDLE;
    VmaAllocation scalarAlloc = VK_NULL_HANDLE;
    VkBuffer scalar4DBuffer = VK_NULL_HANDLE;
    VmaAllocation scalar4DAlloc = VK_NULL_HANDLE;
    VkBuffer vertexBuffer = VK_NULL_HANDLE;
    VmaAllocation vertexAlloc = VK_NULL_HANDLE;
    VkBuffer particleBuffer = VK_NULL_HANDLE;
    VmaAllocation particleAlloc = VK_NULL_HANDLE;
    VkImage fontImage = VK_NULL_HANDLE;
    VmaAllocation fontImageAlloc = VK_NULL_HANDLE;
    VkImageView fontImageView = VK_NULL_HANDLE;
};

bool init_renderer(SDL_Window* window, VulkanContext& ctx, const Simulation& sim);
bool create_swapchain(SDL_Window* window, VulkanContext& ctx);
VkCommandBuffer begin_single_time_commands(VulkanContext& ctx);
void end_single_time_commands(VulkanContext& ctx, VkCommandBuffer commandBuffer);
void cleanup_renderer(VulkanContext& ctx);
void draw_frame(VulkanContext& ctx, const Simulation& sim, const Menu& menu);
void render(const Simulation& sim, VulkanContext& ctx, const Menu& menu);

#endif