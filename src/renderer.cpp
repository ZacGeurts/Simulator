// renderer.cpp
// Copyright (c) 2025, Zachary Geurts

#define VMA_IMPLEMENTATION
#include <vk_mem_alloc.h>
#include "renderer.h"
#include "equations.h" // Assuming Simulation is defined here
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_vulkan.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <fstream>
#include <set>

// Forward declaration for log_message
void log_message(const std::string& message);

// Utility function to read SPIR-V files
std::vector<char> read_spirv_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::ate | std::ios::binary);
    if (!file.is_open()) {
        log_message("Failed to open SPIR-V file: " + filename);
        return std::vector<char>();
    }

    size_t fileSize = static_cast<size_t>(file.tellg());
    std::vector<char> buffer(fileSize);
    file.seekg(0);
    file.read(buffer.data(), fileSize);
    file.close();
    log_message("Successfully read SPIR-V file: " + filename);
    return buffer;
}

// Utility function to create a Vulkan shader module
VkShaderModule create_shader_module(VkDevice device, const std::vector<char>& code) {
    VkShaderModuleCreateInfo createInfo = {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());

    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        log_message("Failed to create shader module");
        return VK_NULL_HANDLE;
    }
    log_message("Shader module created successfully");
    return shaderModule;
}

// Implementation of begin_single_time_commands
VkCommandBuffer begin_single_time_commands(VulkanContext& ctx) {
    VkCommandBufferAllocateInfo allocInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandPool = ctx.commandPool;
    allocInfo.commandBufferCount = 1;

    VkCommandBuffer commandBuffer;
    vkAllocateCommandBuffers(ctx.device, &allocInfo, &commandBuffer);

    VkCommandBufferBeginInfo beginInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(commandBuffer, &beginInfo);

    return commandBuffer;
}

// Implementation of end_single_time_commands
void end_single_time_commands(VulkanContext& ctx, VkCommandBuffer commandBuffer) {
    vkEndCommandBuffer(commandBuffer);

    VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;

    vkQueueSubmit(ctx.graphicsQueue, 1, &submitInfo, VK_NULL_HANDLE);
    vkQueueWaitIdle(ctx.graphicsQueue);

    vkFreeCommandBuffers(ctx.device, ctx.commandPool, 1, &commandBuffer);
}

// Implementation of create_swapchain
bool create_swapchain(SDL_Window* window, VulkanContext& ctx) {
    log_message("Creating swapchain");

    // Get surface capabilities
    VkSurfaceCapabilitiesKHR capabilities;
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(ctx.physicalDevice, ctx.surface, &capabilities);

    // Get surface formats
    uint32_t formatCount;
    vkGetPhysicalDeviceSurfaceFormatsKHR(ctx.physicalDevice, ctx.surface, &formatCount, nullptr);
    std::vector<VkSurfaceFormatKHR> formats(formatCount);
    vkGetPhysicalDeviceSurfaceFormatsKHR(ctx.physicalDevice, ctx.surface, &formatCount, formats.data());

    // Select surface format
    VkSurfaceFormatKHR surfaceFormat = formats[0];
    for (const auto& format : formats) {
        if (format.format == VK_FORMAT_B8G8R8A8_SRGB && format.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR) {
            surfaceFormat = format;
            break;
        }
    }
    ctx.swapchainFormat = surfaceFormat.format;

    // Get present modes
    uint32_t presentModeCount;
    vkGetPhysicalDeviceSurfacePresentModesKHR(ctx.physicalDevice, ctx.surface, &presentModeCount, nullptr);
    std::vector<VkPresentModeKHR> presentModes(presentModeCount);
    vkGetPhysicalDeviceSurfacePresentModesKHR(ctx.physicalDevice, ctx.surface, &presentModeCount, presentModes.data());

    // Select present mode (prefer mailbox, fallback to FIFO)
    VkPresentModeKHR presentMode = VK_PRESENT_MODE_FIFO_KHR;
    for (const auto& mode : presentModes) {
        if (mode == VK_PRESENT_MODE_MAILBOX_KHR) {
            presentMode = mode;
            break;
        }
    }

    // Determine swapchain extent
    VkExtent2D extent = capabilities.currentExtent;
    if (extent.width == UINT32_MAX) {
        int width, height;
        SDL_GetWindowSize(window, &width, &height);
        extent.width = std::clamp(static_cast<uint32_t>(width), capabilities.minImageExtent.width, capabilities.maxImageExtent.width);
        extent.height = std::clamp(static_cast<uint32_t>(height), capabilities.minImageExtent.height, capabilities.maxImageExtent.height);
    }
    ctx.swapchainExtent = extent;

    // Determine image count
    uint32_t imageCount = capabilities.minImageCount + 1;
    if (capabilities.maxImageCount > 0) {
        imageCount = std::min(imageCount, capabilities.maxImageCount);
    }

    // Create swapchain
    VkSwapchainCreateInfoKHR createInfo = {VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    createInfo.surface = ctx.surface;
    createInfo.minImageCount = imageCount;
    createInfo.imageFormat = surfaceFormat.format;
    createInfo.imageColorSpace = surfaceFormat.colorSpace;
    createInfo.imageExtent = extent;
    createInfo.imageArrayLayers = 1;
    createInfo.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    createInfo.preTransform = capabilities.currentTransform;
    createInfo.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    createInfo.presentMode = presentMode;
    createInfo.clipped = VK_TRUE;
    createInfo.oldSwapchain = VK_NULL_HANDLE;

    uint32_t queueFamilyIndices[] = {0, 0}; // Placeholder: Update with actual indices
    if (queueFamilyIndices[0] != queueFamilyIndices[1]) {
        createInfo.imageSharingMode = VK_SHARING_MODE_CONCURRENT;
        createInfo.queueFamilyIndexCount = 2;
        createInfo.pQueueFamilyIndices = queueFamilyIndices;
    } else {
        createInfo.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    }

    if (vkCreateSwapchainKHR(ctx.device, &createInfo, nullptr, &ctx.swapchain) != VK_SUCCESS) {
        log_message("Failed to create swapchain");
        return false;
    }

    // Get swapchain images
    vkGetSwapchainImagesKHR(ctx.device, ctx.swapchain, &imageCount, nullptr);
    ctx.swapchainImages.resize(imageCount);
    vkGetSwapchainImagesKHR(ctx.device, ctx.swapchain, &imageCount, ctx.swapchainImages.data());

    log_message("Swapchain created successfully");
    return true;
}

// Implementation of cleanup_renderer
void cleanup_renderer(VulkanContext& ctx) {
    log_message("Cleaning up renderer");

    if (ctx.device != VK_NULL_HANDLE) {
        vkDeviceWaitIdle(ctx.device);

        ImGui_ImplVulkan_Shutdown();
        ImGui_ImplSDL3_Shutdown();
        ImGui::DestroyContext();

        vkDestroySemaphore(ctx.device, ctx.renderFinishedSemaphore, nullptr);
        vkDestroySemaphore(ctx.device, ctx.imageAvailableSemaphore, nullptr);
        vkDestroyFence(ctx.device, ctx.inFlightFence, nullptr);
        vkDestroyCommandPool(ctx.device, ctx.commandPool, nullptr);
        vkDestroyDescriptorPool(ctx.device, ctx.descriptorPool, nullptr);
        vkDestroyDescriptorSetLayout(ctx.device, ctx.descriptorSetLayout, nullptr);
        vkDestroyRenderPass(ctx.device, ctx.renderPass, nullptr);
        vkDestroyPipeline(ctx.device, ctx.computePipeline, nullptr);
        vkDestroyPipeline(ctx.device, ctx.graphicsPipeline, nullptr);
        vkDestroyPipelineLayout(ctx.device, ctx.pipelineLayout, nullptr);
        vkDestroySwapchainKHR(ctx.device, ctx.swapchain, nullptr);
        vkDestroyImageView(ctx.device, ctx.fontImageView, nullptr);
        if (ctx.fontImage != VK_NULL_HANDLE) {
            vmaDestroyImage(ctx.allocator, ctx.fontImage, ctx.fontImageAlloc);
        }
        if (ctx.particleBuffer != VK_NULL_HANDLE) {
            vmaDestroyBuffer(ctx.allocator, ctx.particleBuffer, ctx.particleAlloc);
        }
        if (ctx.vertexBuffer != VK_NULL_HANDLE) {
            vmaDestroyBuffer(ctx.allocator, ctx.vertexBuffer, ctx.vertexAlloc);
        }
        if (ctx.scalar4DBuffer != VK_NULL_HANDLE) {
            vmaDestroyBuffer(ctx.allocator, ctx.scalar4DBuffer, ctx.scalar4DAlloc);
        }
        if (ctx.scalarBuffer != VK_NULL_HANDLE) {
            vmaDestroyBuffer(ctx.allocator, ctx.scalarBuffer, ctx.scalarAlloc);
        }
        if (ctx.uniformBuffer != VK_NULL_HANDLE) {
            vmaDestroyBuffer(ctx.allocator, ctx.uniformBuffer, ctx.uniformAlloc);
        }

        vkDestroyDevice(ctx.device, nullptr);
    }
    if (ctx.surface != VK_NULL_HANDLE) {
        vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
    }
    if (ctx.instance != VK_NULL_HANDLE) {
        vkDestroyInstance(ctx.instance, nullptr);
    }
    if (ctx.allocator != VK_NULL_HANDLE) {
        vmaDestroyAllocator(ctx.allocator);
    }
}

// Implementation of draw_frame
void draw_frame(VulkanContext& ctx, const Simulation& sim, const Menu& menu) {
    vkWaitForFences(ctx.device, 1, &ctx.inFlightFence, VK_TRUE, UINT64_MAX);
    vkResetFences(ctx.device, 1, &ctx.inFlightFence);

    uint32_t imageIndex;
    vkAcquireNextImageKHR(ctx.device, ctx.swapchain, UINT64_MAX, ctx.imageAvailableSemaphore, VK_NULL_HANDLE, &imageIndex);

    vkResetCommandBuffer(ctx.commandBuffers[imageIndex], 0);
    VkCommandBufferBeginInfo beginInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    vkBeginCommandBuffer(ctx.commandBuffers[imageIndex], &beginInfo);

    // Note: Framebuffer creation is missing. This is a placeholder.
    VkRenderPassBeginInfo renderPassInfo = {VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    renderPassInfo.renderPass = ctx.renderPass;
    renderPassInfo.framebuffer = VK_NULL_HANDLE; // To be fixed
    renderPassInfo.renderArea.offset = {0, 0};
    renderPassInfo.renderArea.extent = ctx.swapchainExtent;
    VkClearValue clearColor = {{{0.0f, 0.0f, 0.0f, 1.0f}}};
    renderPassInfo.clearValueCount = 1;
    renderPassInfo.pClearValues = &clearColor;

    vkCmdBeginRenderPass(ctx.commandBuffers[imageIndex], &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

    VkViewport viewport = {0.0f, 0.0f, (float)ctx.swapchainExtent.width, (float)ctx.swapchainExtent.height, 0.0f, 1.0f};
    vkCmdSetViewport(ctx.commandBuffers[imageIndex], 0, 1, &viewport);
    VkRect2D scissor = {{0, 0}, ctx.swapchainExtent};
    vkCmdSetScissor(ctx.commandBuffers[imageIndex], 0, 1, &scissor);

    // Placeholder: Draw a simple triangle
    vkCmdBindPipeline(ctx.commandBuffers[imageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, ctx.graphicsPipeline);
    VkDeviceSize offset = 0;
    vkCmdDraw(ctx.commandBuffers[imageIndex], 3, 1, 0, 0);

    // Render ImGui
    ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), ctx.commandBuffers[imageIndex]);

    vkCmdEndRenderPass(ctx.commandBuffers[imageIndex]);
    vkEndCommandBuffer(ctx.commandBuffers[imageIndex]);

    VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    VkSemaphore waitSemaphores[] = {ctx.imageAvailableSemaphore};
    VkPipelineStageFlags waitStages[] = {VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT};
    submitInfo.waitSemaphoreCount = 1;
    submitInfo.pWaitSemaphores = waitSemaphores;
    submitInfo.pWaitDstStageMask = waitStages;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &ctx.commandBuffers[imageIndex];
    VkSemaphore signalSemaphores[] = {ctx.renderFinishedSemaphore};
    submitInfo.signalSemaphoreCount = 1;
    submitInfo.pSignalSemaphores = signalSemaphores;

    vkQueueSubmit(ctx.graphicsQueue, 1, &submitInfo, ctx.inFlightFence);

    VkPresentInfoKHR presentInfo = {VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    presentInfo.waitSemaphoreCount = 1;
    presentInfo.pWaitSemaphores = signalSemaphores;
    presentInfo.swapchainCount = 1;
    presentInfo.pSwapchains = &ctx.swapchain;
    presentInfo.pImageIndices = &imageIndex;
    vkQueuePresentKHR(ctx.presentQueue, &presentInfo);
}

// Implementation of render
void render(const Simulation& sim, VulkanContext& ctx, const Menu& menu) {
    log_message("Rendering frame");
    // Assuming Menu has a non-const render method or modify Simulation state elsewhere
    draw_frame(ctx, sim, menu);
}

bool init_renderer(SDL_Window* window, VulkanContext& ctx, const Simulation& sim) {
    log_message("Starting init_renderer");

    // Enable validation layers in debug mode
    std::vector<const char*> layers;
#ifdef _DEBUG
    layers.push_back("VK_LAYER_KHRONOS_validation");
#endif

    // Get required instance extensions from SDL
    uint32_t extensionCount = 0;
    if (!SDL_Vulkan_GetInstanceExtensions(&extensionCount)) {
        log_message("Failed to get Vulkan instance extension count: " + std::string(SDL_GetError()));
        return false;
    }
    std::vector<const char*> extensionList(extensionCount);
#ifdef _DEBUG
    extensionList.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
#endif

    // Create Vulkan instance
    VkApplicationInfo appInfo = {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    appInfo.pApplicationName = "Physics Simulator";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_3;

    VkInstanceCreateInfo createInfo = {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    createInfo.pApplicationInfo = &appInfo;
    createInfo.enabledLayerCount = static_cast<uint32_t>(layers.size());
    createInfo.ppEnabledLayerNames = layers.empty() ? nullptr : layers.data();
    createInfo.enabledExtensionCount = static_cast<uint32_t>(extensionList.size());
    createInfo.ppEnabledExtensionNames = extensionList.data();

    log_message("Creating Vulkan instance");
    if (vkCreateInstance(&createInfo, nullptr, &ctx.instance) != VK_SUCCESS) {
        log_message("Failed to create Vulkan instance");
        return false;
    }
    log_message("Vulkan instance created");

    // Create Vulkan surface
    if (!SDL_Vulkan_CreateSurface(window, ctx.instance, nullptr, &ctx.surface)) {
        log_message("Failed to create Vulkan surface: " + std::string(SDL_GetError()));
        vkDestroyInstance(ctx.instance, nullptr);
        return false;
    }
    log_message("Vulkan surface created");

    // Select physical device
    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, nullptr);
    if (deviceCount == 0) {
        log_message("No Vulkan-capable devices found");
        vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
        vkDestroyInstance(ctx.instance, nullptr);
        return false;
    }
    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, devices.data());
    ctx.physicalDevice = devices[0]; // Simplified: pick first device

    // Find queue families
    uint32_t queueFamilyCount = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(ctx.physicalDevice, &queueFamilyCount, nullptr);
    std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
    vkGetPhysicalDeviceQueueFamilyProperties(ctx.physicalDevice, &queueFamilyCount, queueFamilies.data());
    uint32_t graphicsFamily = UINT32_MAX, computeFamily = UINT32_MAX, presentFamily = UINT32_MAX;
    for (uint32_t i = 0; i < queueFamilyCount; ++i) {
        if (queueFamilies[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) graphicsFamily = i;
        if (queueFamilies[i].queueFlags & VK_QUEUE_COMPUTE_BIT) computeFamily = i;
        VkBool32 presentSupport = false;
        vkGetPhysicalDeviceSurfaceSupportKHR(ctx.physicalDevice, i, ctx.surface, &presentSupport);
        if (presentSupport) presentFamily = i;
    }
    if (graphicsFamily == UINT32_MAX || computeFamily == UINT32_MAX || presentFamily == UINT32_MAX) {
        log_message("Failed to find required queue families");
        vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
        vkDestroyInstance(ctx.instance, nullptr);
        return false;
    }

    // Create logical device
    float queuePriority = 1.0f;
    std::vector<VkDeviceQueueCreateInfo> queueCreateInfos;
    std::set<uint32_t> uniqueFamilies = {graphicsFamily, computeFamily, presentFamily};
    for (uint32_t family : uniqueFamilies) {
        VkDeviceQueueCreateInfo queueCreateInfo = {VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        queueCreateInfo.queueFamilyIndex = family;
        queueCreateInfo.queueCount = 1;
        queueCreateInfo.pQueuePriorities = &queuePriority;
        queueCreateInfos.push_back(queueCreateInfo);
    }
    VkPhysicalDeviceFeatures deviceFeatures = {};
    const char* deviceExtensions[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    VkDeviceCreateInfo deviceCreateInfo = {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceCreateInfo.queueCreateInfoCount = static_cast<uint32_t>(queueCreateInfos.size());
    deviceCreateInfo.pQueueCreateInfos = queueCreateInfos.data();
    deviceCreateInfo.pEnabledFeatures = &deviceFeatures;
    deviceCreateInfo.enabledExtensionCount = 1;
    deviceCreateInfo.ppEnabledExtensionNames = deviceExtensions;
    if (vkCreateDevice(ctx.physicalDevice, &deviceCreateInfo, nullptr, &ctx.device) != VK_SUCCESS) {
        log_message("Failed to create Vulkan device");
        vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
        vkDestroyInstance(ctx.instance, nullptr);
        return false;
    }
    vkGetDeviceQueue(ctx.device, graphicsFamily, 0, &ctx.graphicsQueue);
    vkGetDeviceQueue(ctx.device, computeFamily, 0, &ctx.computeQueue);
    vkGetDeviceQueue(ctx.device, presentFamily, 0, &ctx.presentQueue);
    log_message("Vulkan device and queues created");

    // Create VMA allocator
    VmaAllocatorCreateInfo allocatorInfo = {};
    allocatorInfo.vulkanApiVersion = VK_API_VERSION_1_3;
    allocatorInfo.physicalDevice = ctx.physicalDevice;
    allocatorInfo.device = ctx.device;
    allocatorInfo.instance = ctx.instance;
    if (vmaCreateAllocator(&allocatorInfo, &ctx.allocator) != VK_SUCCESS) {
        log_message("Failed to create VMA allocator");
        vkDestroyDevice(ctx.device, nullptr);
        vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
        vkDestroyInstance(ctx.instance, nullptr);
        return false;
    }

    // Create swapchain
    if (!create_swapchain(window, ctx)) {
        log_message("Failed to create swapchain");
        cleanup_renderer(ctx);
        return false;
    }

    // Create render pass
    VkAttachmentDescription colorAttachment = {};
    colorAttachment.format = VK_FORMAT_B8G8R8A8_SRGB;
    colorAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
    colorAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    colorAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    colorAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    colorAttachment.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    VkAttachmentReference colorAttachmentRef = {0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkSubpassDescription subpass = {};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &colorAttachmentRef;
    VkRenderPassCreateInfo renderPassInfo = {VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    renderPassInfo.attachmentCount = 1;
    renderPassInfo.pAttachments = &colorAttachment;
    renderPassInfo.subpassCount = 1;
    renderPassInfo.pSubpasses = &subpass;
    if (vkCreateRenderPass(ctx.device, &renderPassInfo, nullptr, &ctx.renderPass) != VK_SUCCESS) {
        log_message("Failed to create render pass");
        cleanup_renderer(ctx);
        return false;
    }

    // Create descriptor set layout
    VkDescriptorSetLayoutBinding bindings[4] = {
        {0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {3, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT | VK_SHADER_STAGE_VERTEX_BIT, nullptr}
    };
    VkDescriptorSetLayoutCreateInfo layoutInfo = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    layoutInfo.bindingCount = 4;
    layoutInfo.pBindings = bindings;
    if (vkCreateDescriptorSetLayout(ctx.device, &layoutInfo, nullptr, &ctx.descriptorSetLayout) != VK_SUCCESS) {
        log_message("Failed to create descriptor set layout");
        cleanup_renderer(ctx);
        return false;
    }

    // Create descriptor pool
    uint32_t imageCount = static_cast<uint32_t>(ctx.swapchainImages.size());
    VkDescriptorPoolSize poolSizes[] = {
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, imageCount},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, imageCount * 3},
        {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, imageCount} // For ImGui
    };
    VkDescriptorPoolCreateInfo descriptorPoolInfo = {VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    descriptorPoolInfo.poolSizeCount = 3;
    descriptorPoolInfo.pPoolSizes = poolSizes;
    descriptorPoolInfo.maxSets = imageCount + 1; // +1 for ImGui
    descriptorPoolInfo.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    if (vkCreateDescriptorPool(ctx.device, &descriptorPoolInfo, nullptr, &ctx.descriptorPool) != VK_SUCCESS) {
        log_message("Failed to create descriptor pool");
        cleanup_renderer(ctx);
        return false;
    }

    // Allocate descriptor sets
    ctx.descriptorSets.resize(imageCount);
    std::vector<VkDescriptorSetLayout> layouts(imageCount, ctx.descriptorSetLayout);
    VkDescriptorSetAllocateInfo descriptorAllocInfo = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    descriptorAllocInfo.descriptorPool = ctx.descriptorPool;
    descriptorAllocInfo.descriptorSetCount = imageCount;
    descriptorAllocInfo.pSetLayouts = layouts.data();
    if (vkAllocateDescriptorSets(ctx.device, &descriptorAllocInfo, ctx.descriptorSets.data()) != VK_SUCCESS) {
        log_message("Failed to allocate descriptor sets");
        cleanup_renderer(ctx);
        return false;
    }

    // Create pipelines
    std::vector<char> vertexShaderCode = read_spirv_file("../shaders/vert.spv");
    std::vector<char> fragmentShaderCode = read_spirv_file("../shaders/frag.spv");
    std::vector<char> computeShaderCode = read_spirv_file("../shaders/compute.comp.spv");
    VkShaderModule vertexShaderModule = create_shader_module(ctx.device, vertexShaderCode);
    VkShaderModule fragmentShaderModule = create_shader_module(ctx.device, fragmentShaderCode);
    VkShaderModule computeShaderModule = create_shader_module(ctx.device, computeShaderCode);

    if (vertexShaderModule == VK_NULL_HANDLE || fragmentShaderModule == VK_NULL_HANDLE || computeShaderModule == VK_NULL_HANDLE) {
        log_message("Failed to create one or more shader modules");
        cleanup_renderer(ctx);
        return false;
    }

    VkPipelineShaderStageCreateInfo shaderStages[2] = {
        {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, vertexShaderModule, "main"},
        {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, fragmentShaderModule, "main"}
    };
    VkPipelineVertexInputStateCreateInfo vertexInputInfo = {VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    VkVertexInputAttributeDescription attributes[2] = {
        {0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0}, // Position
        {1, 0, VK_FORMAT_R32G32B32_SFLOAT, 12} // Color
    };
    VkVertexInputBindingDescription bindingDescription = {0, 24, VK_VERTEX_INPUT_RATE_VERTEX};
    vertexInputInfo.vertexBindingDescriptionCount = 1;
    vertexInputInfo.pVertexBindingDescriptions = &bindingDescription;
    vertexInputInfo.vertexAttributeDescriptionCount = 2;
    vertexInputInfo.pVertexAttributeDescriptions = attributes;
    VkPipelineInputAssemblyStateCreateInfo inputAssembly = {VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    inputAssembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    VkPipelineViewportStateCreateInfo viewportState = {VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewportState.viewportCount = 1;
    viewportState.scissorCount = 1;
    VkPipelineRasterizationStateCreateInfo rasterizer = {VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    rasterizer.lineWidth = 1.0f;
    rasterizer.cullMode = VK_CULL_MODE_NONE;
    VkPipelineMultisampleStateCreateInfo multisampling = {VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisampling.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;
    VkPipelineColorBlendAttachmentState colorBlendAttachment = {};
    colorBlendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    colorBlendAttachment.blendEnable = VK_FALSE;
    VkPipelineColorBlendStateCreateInfo colorBlending = {VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    colorBlending.attachmentCount = 1;
    colorBlending.pAttachments = &colorBlendAttachment;
    VkPipelineLayoutCreateInfo pipelineLayoutInfo = {VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pipelineLayoutInfo.setLayoutCount = 1;
    pipelineLayoutInfo.pSetLayouts = &ctx.descriptorSetLayout;
    VkPushConstantRange pushConstantRange = {VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(float) * 6 + sizeof(int) * 2};
    pipelineLayoutInfo.pushConstantRangeCount = 1;
    pipelineLayoutInfo.pPushConstantRanges = &pushConstantRange;
    if (vkCreatePipelineLayout(ctx.device, &pipelineLayoutInfo, nullptr, &ctx.pipelineLayout) != VK_SUCCESS) {
        log_message("Failed to create pipeline layout");
        cleanup_renderer(ctx);
        return false;
    }
    VkPipelineDynamicStateCreateInfo dynamicState = {VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    VkDynamicState dynamicStates[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    dynamicState.dynamicStateCount = 2;
    dynamicState.pDynamicStates = dynamicStates;
    VkGraphicsPipelineCreateInfo pipelineInfo = {VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    pipelineInfo.stageCount = 2;
    pipelineInfo.pStages = shaderStages;
    pipelineInfo.pVertexInputState = &vertexInputInfo;
    pipelineInfo.pInputAssemblyState = &inputAssembly;
    pipelineInfo.pViewportState = &viewportState;
    pipelineInfo.pRasterizationState = &rasterizer;
    pipelineInfo.pMultisampleState = &multisampling;
    pipelineInfo.pColorBlendState = &colorBlending;
    pipelineInfo.pDynamicState = &dynamicState;
    pipelineInfo.renderPass = ctx.renderPass;
    pipelineInfo.layout = ctx.pipelineLayout;
    if (vkCreateGraphicsPipelines(ctx.device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &ctx.graphicsPipeline) != VK_SUCCESS) {
        log_message("Failed to create graphics pipeline");
        cleanup_renderer(ctx);
        return false;
    }
    VkComputePipelineCreateInfo computePipelineInfo = {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    computePipelineInfo.stage = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_COMPUTE_BIT, computeShaderModule, "main"};
    computePipelineInfo.layout = ctx.pipelineLayout;
    if (vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &computePipelineInfo, nullptr, &ctx.computePipeline) != VK_SUCCESS) {
        log_message("Failed to create compute pipeline");
        cleanup_renderer(ctx);
        return false;
    }
    vkDestroyShaderModule(ctx.device, vertexShaderModule, nullptr);
    vkDestroyShaderModule(ctx.device, fragmentShaderModule, nullptr);
    vkDestroyShaderModule(ctx.device, computeShaderModule, nullptr);

    // Create buffers
    VkBufferCreateInfo bufferCreateInfo = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    VmaAllocationCreateInfo allocCreateInfo = {};
    allocCreateInfo.usage = VMA_MEMORY_USAGE_AUTO;
    allocCreateInfo.flags = VMA_ALLOCATION_CREATE_HOST_ACCESS_SEQUENTIAL_WRITE_BIT;

    bufferCreateInfo.size = sizeof(glm::mat4) + sizeof(float) * 2 + sizeof(int);
    bufferCreateInfo.usage = VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    if (vmaCreateBuffer(ctx.allocator, &bufferCreateInfo, &allocCreateInfo, &ctx.uniformBuffer, &ctx.uniformAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create uniform buffer");
        cleanup_renderer(ctx);
        return false;
    }

    bufferCreateInfo.size = sim.size * sim.size * sim.size * sizeof(float);
    bufferCreateInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    if (vmaCreateBuffer(ctx.allocator, &bufferCreateInfo, &allocCreateInfo, &ctx.scalarBuffer, &ctx.scalarAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create scalar buffer");
        cleanup_renderer(ctx);
        return false;
    }

    bufferCreateInfo.size = sim.size * sim.size * sim.size * sim.size * sizeof(float);
    if (vmaCreateBuffer(ctx.allocator, &bufferCreateInfo, &allocCreateInfo, &ctx.scalar4DBuffer, &ctx.scalar4DAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create scalar_4d buffer");
        cleanup_renderer(ctx);
        return false;
    }

    bufferCreateInfo.size = sim.size * sim.size * sim.size * 15 * sizeof(float) * 6;
    bufferCreateInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    if (vmaCreateBuffer(ctx.allocator, &bufferCreateInfo, &allocCreateInfo, &ctx.vertexBuffer, &ctx.vertexAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create vertex buffer");
        cleanup_renderer(ctx);
        return false;
    }

    bufferCreateInfo.size = sim.particles.size() * 6 * sizeof(float);
    if (vmaCreateBuffer(ctx.allocator, &bufferCreateInfo, &allocCreateInfo, &ctx.particleBuffer, &ctx.particleAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create particle buffer");
        cleanup_renderer(ctx);
        return false;
    }

    // Update descriptor sets
    for (uint32_t i = 0; i < imageCount; ++i) {
        VkDescriptorBufferInfo uniformInfo = {ctx.uniformBuffer, 0, sizeof(glm::mat4) + sizeof(float) * 2 + sizeof(int)};
        VkDescriptorBufferInfo scalarInfo = {ctx.scalarBuffer, 0, sim.size * sim.size * sim.size * sizeof(float)};
        VkDescriptorBufferInfo scalar4DInfo = {ctx.scalar4DBuffer, 0, sim.size * sim.size * sim.size * sim.size * sizeof(float)};
        VkDescriptorBufferInfo vertexInfo = {ctx.vertexBuffer, 0, sim.size * sim.size * sim.size * 15 * sizeof(float) * 6};
        VkWriteDescriptorSet descriptorWrites[4] = {
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, ctx.descriptorSets[i], 0, 0, 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, nullptr, &uniformInfo, nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, ctx.descriptorSets[i], 1, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &scalarInfo, nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, ctx.descriptorSets[i], 2, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &scalar4DInfo, nullptr},
            {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET, nullptr, ctx.descriptorSets[i], 3, 0, 1, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, nullptr, &vertexInfo, nullptr}
        };
        vkUpdateDescriptorSets(ctx.device, 4, descriptorWrites, 0, nullptr);
    }

    // Create command pool and buffers
    VkCommandPoolCreateInfo commandPoolInfo = {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    commandPoolInfo.queueFamilyIndex = graphicsFamily;
    commandPoolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    if (vkCreateCommandPool(ctx.device, &commandPoolInfo, nullptr, &ctx.commandPool) != VK_SUCCESS) {
        log_message("Failed to create command pool");
        cleanup_renderer(ctx);
        return false;
    }
    ctx.commandBuffers.resize(imageCount);
    VkCommandBufferAllocateInfo cmdAllocInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cmdAllocInfo.commandPool = ctx.commandPool;
    cmdAllocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cmdAllocInfo.commandBufferCount = imageCount;
    if (vkAllocateCommandBuffers(ctx.device, &cmdAllocInfo, ctx.commandBuffers.data()) != VK_SUCCESS) {
        log_message("Failed to allocate command buffers");
        cleanup_renderer(ctx);
        return false;
    }

    // Create synchronization objects
    VkSemaphoreCreateInfo semaphoreInfo = {VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    if (vkCreateSemaphore(ctx.device, &semaphoreInfo, nullptr, &ctx.imageAvailableSemaphore) != VK_SUCCESS ||
        vkCreateSemaphore(ctx.device, &semaphoreInfo, nullptr, &ctx.renderFinishedSemaphore) != VK_SUCCESS) {
        log_message("Failed to create semaphores");
        cleanup_renderer(ctx);
        return false;
    }
    VkFenceCreateInfo fenceInfo = {VK_STRUCTURE_TYPE_FENCE_CREATE_INFO, nullptr, VK_FENCE_CREATE_SIGNALED_BIT};
    if (vkCreateFence(ctx.device, &fenceInfo, nullptr, &ctx.inFlightFence) != VK_SUCCESS) {
        log_message("Failed to create fence");
        cleanup_renderer(ctx);
        return false;
    }

    // Initialize ImGui
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.5f;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    if (!ImGui_ImplSDL3_InitForVulkan(window)) {
        log_message("Failed to initialize ImGui SDL backend");
        cleanup_renderer(ctx);
        return false;
    }
    ImGui_ImplVulkan_InitInfo init_info = {};
    init_info.Instance = ctx.instance;
    init_info.PhysicalDevice = ctx.physicalDevice;
    init_info.Device = ctx.device;
    init_info.QueueFamily = graphicsFamily;
    init_info.Queue = ctx.graphicsQueue;
    init_info.PipelineCache = VK_NULL_HANDLE;
    init_info.DescriptorPool = ctx.descriptorPool;
    init_info.Subpass = 0;
    init_info.MinImageCount = imageCount;
    init_info.ImageCount = imageCount;
    init_info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;
    init_info.Allocator = nullptr;
    init_info.CheckVkResultFn = nullptr;
    init_info.RenderPass = ctx.renderPass;
    if (!ImGui_ImplVulkan_Init(&init_info)) {
        log_message("Failed to initialize ImGui Vulkan backend");
        cleanup_renderer(ctx);
        return false;
    }
    VkCommandBuffer commandBuffer = begin_single_time_commands(ctx);
    // Upload ImGui font atlas
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    VkDeviceSize imageSize = width * height * 4;

    // Create staging buffer
    VkBuffer stagingBuffer;
    VmaAllocation stagingAlloc;
    VkBufferCreateInfo stagingBufferInfo = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    stagingBufferInfo.size = imageSize;
    stagingBufferInfo.usage = VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
    VmaAllocationCreateInfo stagingAllocInfo = {};
    stagingAllocInfo.usage = VMA_MEMORY_USAGE_AUTO;
    stagingAllocInfo.flags = VMA_ALLOCATION_CREATE_HOST_ACCESS_SEQUENTIAL_WRITE_BIT;
    if (vmaCreateBuffer(ctx.allocator, &stagingBufferInfo, &stagingAllocInfo, &stagingBuffer, &stagingAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create staging buffer for ImGui fonts");
        end_single_time_commands(ctx, commandBuffer);
        cleanup_renderer(ctx);
        return false;
    }

    // Copy font data to staging buffer
    void* data;
    vmaMapMemory(ctx.allocator, stagingAlloc, &data);
    memcpy(data, pixels, imageSize);
    vmaUnmapMemory(ctx.allocator, stagingAlloc);

    // Create image for font atlas
    VkImage fontImage;
    VmaAllocation fontImageAlloc;
    VkImageCreateInfo fontImageInfo = {VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    fontImageInfo.imageType = VK_IMAGE_TYPE_2D;
    fontImageInfo.extent.width = width;
    fontImageInfo.extent.height = height;
    fontImageInfo.extent.depth = 1;
    fontImageInfo.mipLevels = 1;
    fontImageInfo.arrayLayers = 1;
    fontImageInfo.format = VK_FORMAT_R8G8B8A8_UNORM;
    fontImageInfo.tiling = VK_IMAGE_TILING_OPTIMAL;
    fontImageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    fontImageInfo.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
    fontImageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
    fontImageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    if (vmaCreateImage(ctx.allocator, &fontImageInfo, &stagingAllocInfo, &fontImage, &fontImageAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create font image");
        vmaDestroyBuffer(ctx.allocator, stagingBuffer, stagingAlloc);
        end_single_time_commands(ctx, commandBuffer);
        cleanup_renderer(ctx);
        return false;
    }

    // Transition image layout and copy data
    VkCommandBufferBeginInfo beginInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    vkBeginCommandBuffer(commandBuffer, &beginInfo);
    VkImageMemoryBarrier barrier = {VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    barrier.image = fontImage;
    barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    barrier.subresourceRange.levelCount = 1;
    barrier.subresourceRange.layerCount = 1;
    barrier.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    barrier.srcAccessMask = 0;
    barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    vkCmdPipelineBarrier(commandBuffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 0, nullptr, 0, nullptr, 1, &barrier);
    VkBufferImageCopy copyRegion = {};
    copyRegion.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copyRegion.imageSubresource.layerCount = 1;
    copyRegion.imageExtent = {static_cast<uint32_t>(width), static_cast<uint32_t>(height), 1};
    vkCmdCopyBufferToImage(commandBuffer, stagingBuffer, fontImage, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &copyRegion);
    barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    barrier.newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
    barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    vkCmdPipelineBarrier(commandBuffer, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 0, nullptr, 0, nullptr, 1, &barrier);
    vkEndCommandBuffer(commandBuffer);

    // Submit command buffer
    VkSubmitInfo submitInfo = {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;
    vkQueueSubmit(ctx.graphicsQueue, 1, &submitInfo, VK_NULL_HANDLE);
    vkQueueWaitIdle(ctx.graphicsQueue);

    // Create image view
    VkImageView fontImageView;
    VkImageViewCreateInfo viewInfo = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    viewInfo.image = fontImage;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = VK_FORMAT_R8G8B8A8_UNORM;
    viewInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    viewInfo.subresourceRange.levelCount = 1;
    viewInfo.subresourceRange.layerCount = 1;
    if (vkCreateImageView(ctx.device, &viewInfo, nullptr, &fontImageView) != VK_SUCCESS) {
        log_message("Failed to create font image view");
        vmaDestroyImage(ctx.allocator, fontImage, fontImageAlloc);
        vmaDestroyBuffer(ctx.allocator, stagingBuffer, stagingAlloc);
        end_single_time_commands(ctx, commandBuffer);
        cleanup_renderer(ctx);
        return false;
    }

    // Update ImGui font texture ID
    io.Fonts->SetTexID((ImTextureID)fontImageView);

    // Store for cleanup
    ctx.fontImage = fontImage;
    ctx.fontImageAlloc = fontImageAlloc;
    ctx.fontImageView = fontImageView;

    // Clean up
    vmaDestroyBuffer(ctx.allocator, stagingBuffer, stagingAlloc);
    end_single_time_commands(ctx, commandBuffer);

    return true;
}