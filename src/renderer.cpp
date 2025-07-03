#include "renderer.h"
#include "menu.h"
#include <SDL3/SDL_vulkan.h>
#include <vk_mem_alloc.h> // Updated include for VMA (header-only)
#include <imgui.h>
#include <imgui_impl_sdl3.h>
#include <imgui_impl_vulkan.h>
#include <stdexcept>
#include <array>
#include <cstring>
#include <fstream>
#include "tables.h"

// Global variables
float camera_zoom = 0.3f, camera_angle = 0.0f, camera_tilt = 0.0f;
DisplayMode current_display_mode = POINTS;

// Logging function
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

// Load SPIR-V shader from file
std::vector<char> read_spirv_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open shader file: " + filename);
    }
    size_t fileSize = file.tellg();
    std::vector<char> buffer(fileSize);
    file.seekg(0);
    file.read(buffer.data(), fileSize);
    file.close();
    return buffer;
}

VkShaderModule create_shader_module(VkDevice device, const std::vector<char>& code) {
    VkShaderModuleCreateInfo createInfo = {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t*>(code.data());
    VkShaderModule shaderModule;
    if (vkCreateShaderModule(device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create shader module");
    }
    return shaderModule;
}

bool init_renderer(SDL_Window* window, VulkanContext& ctx, const Simulation& sim) {
    VkApplicationInfo appInfo = {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    appInfo.pApplicationName = "Physics Simulator";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_3;
    uint32_t extensionCount;
    SDL_Vulkan_GetInstanceExtensions(window, &extensionCount, nullptr);
    std::vector<const char*> extensions(extensionCount);
    SDL_Vulkan_GetInstanceExtensions(window, &extensionCount, extensions.data());
    VkInstanceCreateInfo createInfo = {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    createInfo.pApplicationInfo = &appInfo;
    createInfo.enabledExtensionCount = extensionCount;
    createInfo.ppEnabledExtensionNames = extensions.data();
    if (vkCreateInstance(&createInfo, nullptr, &ctx.instance) != VK_SUCCESS) {
        log_message("Failed to create Vulkan instance");
        return false;
    }

    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, nullptr);
    if (deviceCount == 0) {
        log_message("No Vulkan-capable devices found");
        return false;
    }
    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(ctx.instance, &deviceCount, devices.data());
    ctx.physicalDevice = devices[0];
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
        return false;
    }

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
    VkDeviceCreateInfo deviceCreateInfo = {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceCreateInfo.queueCreateInfoCount = static_cast<uint32_t>(queueCreateInfos.size());
    deviceCreateInfo.pQueueCreateInfos = queueCreateInfos.data();
    deviceCreateInfo.pEnabledFeatures = &deviceFeatures;
    const char* deviceExtensions[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    deviceCreateInfo.enabledExtensionCount = 1;
    deviceCreateInfo.ppEnabledExtensionNames = deviceExtensions;
    if (vkCreateDevice(ctx.physicalDevice, &deviceCreateInfo, nullptr, &ctx.device) != VK_SUCCESS) {
        log_message("Failed to create Vulkan device");
        return false;
    }
    vkGetDeviceQueue(ctx.device, graphicsFamily, 0, &ctx.graphicsQueue);
    vkGetDeviceQueue(ctx.device, computeFamily, 0, &ctx.computeQueue);
    vkGetDeviceQueue(ctx.device, presentFamily, 0, &ctx.presentQueue);

    VkSurfaceKHR surface;
    if (!SDL_Vulkan_CreateSurface(window, ctx.instance, nullptr, &surface)) {
        log_message("Failed to create Vulkan surface");
        return false;
    }
    ctx.surface = surface;

    VkSwapchainCreateInfoKHR swapchainInfo = {VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    swapchainInfo.surface = ctx.surface;
    swapchainInfo.minImageCount = 2;
    swapchainInfo.imageFormat = VK_FORMAT_B8G8R8A8_SRGB;
    swapchainInfo.imageColorSpace = VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
    swapchainInfo.imageExtent = {800, 600};
    swapchainInfo.imageArrayLayers = 1;
    swapchainInfo.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    swapchainInfo.preTransform = VK_SURFACE_TRANSFORM_IDENTITY_BIT_KHR;
    swapchainInfo.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    swapchainInfo.presentMode = VK_PRESENT_MODE_FIFO_KHR;
    if (vkCreateSwapchainKHR(ctx.device, &swapchainInfo, nullptr, &ctx.swapchain) != VK_SUCCESS) {
        log_message("Failed to create swapchain");
        return false;
    }
    uint32_t imageCount;
    vkGetSwapchainImagesKHR(ctx.device, ctx.swapchain, &imageCount, nullptr);
    ctx.swapchainImages.resize(imageCount);
    vkGetSwapchainImagesKHR(ctx.device, ctx.swapchain, &imageCount, ctx.swapchainImages.data());
    ctx.swapchainImageViews.resize(imageCount);
    for (uint32_t i = 0; i < imageCount; ++i) {
        VkImageViewCreateInfo viewInfo = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        viewInfo.image = ctx.swapchainImages[i];
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = VK_FORMAT_B8G8R8A8_SRGB;
        viewInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        viewInfo.subresourceRange.levelCount = 1;
        viewInfo.subresourceRange.layerCount = 1;
        if (vkCreateImageView(ctx.device, &viewInfo, nullptr, &ctx.swapchainImageViews[i]) != VK_SUCCESS) {
            log_message("Failed to create image view");
            return false;
        }
    }

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
        return false;
    }

    ctx.framebuffers.resize(imageCount);
    for (uint32_t i = 0; i < imageCount; ++i) {
        VkFramebufferCreateInfo framebufferInfo = {VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        framebufferInfo.renderPass = ctx.renderPass;
        framebufferInfo.attachmentCount = 1;
        framebufferInfo.pAttachments = &ctx.swapchainImageViews[i];
        framebufferInfo.width = 800;
        framebufferInfo.height = 600;
        framebufferInfo.layers = 1;
        if (vkCreateFramebuffer(ctx.device, &framebufferInfo, nullptr, &ctx.framebuffers[i]) != VK_SUCCESS) {
            log_message("Failed to create framebuffer");
            return false;
        }
    }

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
        return false;
    }

    VkDescriptorPoolSize poolSizes[] = {
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, imageCount},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, imageCount * 3}
    };
    VkDescriptorPoolCreateInfo poolInfo = {VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    poolInfo.poolSizeCount = 2;
    poolInfo.pPoolSizes = poolSizes;
    poolInfo.maxSets = imageCount;
    if (vkCreateDescriptorPool(ctx.device, &poolInfo, nullptr, &ctx.descriptorPool) != VK_SUCCESS) {
        log_message("Failed to create descriptor pool");
        return false;
    }

    ctx.descriptorSets.resize(imageCount);
    std::vector<VkDescriptorSetLayout> layouts(imageCount, ctx.descriptorSetLayout);
    VkDescriptorSetAllocateInfo allocInfo = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    allocInfo.descriptorPool = ctx.descriptorPool;
    allocInfo.descriptorSetCount = imageCount;
    allocInfo.pSetLayouts = layouts.data();
    if (vkAllocateDescriptorSets(ctx.device, &allocInfo, ctx.descriptorSets.data()) != VK_SUCCESS) {
        log_message("Failed to allocate descriptor sets");
        return false;
    }

    // Load SPIR-V shaders
    std::vector<char> vertexShaderCode = read_spirv_file("../vert.spv");
    std::vector<char> fragmentShaderCode = read_spirv_file("../frag.spv");
    std::vector<char> computeShaderCode = read_spirv_file("../compute.comp.spv");

    VkPipelineShaderStageCreateInfo shaderStages[2] = {
        {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_VERTEX_BIT, create_shader_module(ctx.device, vertexShaderCode), "main"},
        {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_FRAGMENT_BIT, create_shader_module(ctx.device, fragmentShaderCode), "main"}
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
    VkViewport viewport = {0.0f, 0.0f, 800.0f, 600.0f, 0.0f, 1.0f};
    VkRect2D scissor = {{0, 0}, {800, 600}};
    VkPipelineViewportStateCreateInfo viewportState = {VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewportState.viewportCount = 1;
    viewportState.pViewports = &viewport;
    viewportState.scissorCount = 1;
    viewportState.pScissors = &scissor;
    VkPipelineRasterizationStateCreateInfo rasterizer = {VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    rasterizer.lineWidth = 1.0f;
    rasterizer.cullMode = VK_CULL_MODE_NONE; // Disable culling to avoid missing triangles
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
        return false;
    }
    VkGraphicsPipelineCreateInfo pipelineInfo = {VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    pipelineInfo.stageCount = 2;
    pipelineInfo.pStages = shaderStages;
    pipelineInfo.pVertexInputState = &vertexInputInfo;
    pipelineInfo.pInputAssemblyState = &inputAssembly;
    pipelineInfo.pViewportState = &viewportState;
    pipelineInfo.pRasterizationState = &rasterizer;
    pipelineInfo.pMultisampleState = &multisampling;
    pipelineInfo.pColorBlendState = &colorBlending;
    pipelineInfo.renderPass = ctx.renderPass;
    pipelineInfo.layout = ctx.pipelineLayout;
    if (vkCreateGraphicsPipelines(ctx.device, VK_NULL_HANDLE, 1, &pipelineInfo, nullptr, &ctx.graphicsPipeline) != VK_SUCCESS) {
        log_message("Failed to create graphics pipeline");
        return false;
    }
    VkComputePipelineCreateInfo computePipelineInfo = {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    computePipelineInfo.stage = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO, nullptr, 0, VK_SHADER_STAGE_COMPUTE_BIT, create_shader_module(ctx.device, computeShaderCode), "main"};
    computePipelineInfo.layout = ctx.pipelineLayout;
    if (vkCreateComputePipelines(ctx.device, VK_NULL_HANDLE, 1, &computePipelineInfo, nullptr, &ctx.computePipeline) != VK_SUCCESS) {
        log_message("Failed to create compute pipeline");
        return false;
    }

    VmaAllocatorCreateInfo allocatorInfo = {};
    allocatorInfo.vulkanApiVersion = VK_API_VERSION_1_3;
    allocatorInfo.physicalDevice = ctx.physicalDevice;
    allocatorInfo.device = ctx.device;
    allocatorInfo.instance = ctx.instance;
    if (vmaCreateAllocator(&allocatorInfo, &ctx.allocator) != VK_SUCCESS) {
        log_message("Failed to create VMA allocator");
        return false;
    }

    // Buffer creation with optimized sizes
    VkBufferCreateInfo bufferInfo = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bufferInfo.size = sizeof(glm::mat4) + sizeof(float) * 2 + sizeof(int);
    bufferInfo.usage = VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    VmaAllocationCreateInfo allocInfo = {0, VMA_MEMORY_USAGE_AUTO, VMA_ALLOCATION_CREATE_HOST_ACCESS_SEQUENTIAL_WRITE_BIT};
    if (vmaCreateBuffer(ctx.allocator, &bufferInfo, &allocInfo, &ctx.uniformBuffer, &ctx.uniformAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create uniform buffer");
        return false;
    }
    bufferInfo.size = sim.size * sim.size * sim.size * sizeof(float);
    bufferInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    if (vmaCreateBuffer(ctx.allocator, &bufferInfo, &allocInfo, &ctx.scalarBuffer, &ctx.scalarAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create scalar buffer");
        return false;
    }
    bufferInfo.size = sim.size * sim.size * sim.size * sim.size * sizeof(float);
    if (vmaCreateBuffer(ctx.allocator, &bufferInfo, &allocInfo, &ctx.scalar4DBuffer, &ctx.scalar4DAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create scalar_4d buffer");
        return false;
    }
    // Vertex buffer size: up to 15 vertices per cube (5 triangles * 3 vertices)
    bufferInfo.size = sim.size * sim.size * sim.size * 15 * sizeof(float) * 6; // pos (3) + color (3)
    bufferInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    if (vmaCreateBuffer(ctx.allocator, &bufferInfo, &allocInfo, &ctx.vertexBuffer, &ctx.vertexAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create vertex buffer");
        return false;
    }
    bufferInfo.size = sim.particles.size() * 6 * sizeof(float);
    if (vmaCreateBuffer(ctx.allocator, &bufferInfo, &allocInfo, &ctx.particleBuffer, &ctx.particleAlloc, nullptr) != VK_SUCCESS) {
        log_message("Failed to create particle buffer");
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

    VkCommandPoolCreateInfo poolInfo = {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    poolInfo.queueFamilyIndex = graphicsFamily;
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    if (vkCreateCommandPool(ctx.device, &poolInfo, nullptr, &ctx.commandPool) != VK_SUCCESS) {
        log_message("Failed to create command pool");
        return false;
    }
    ctx.commandBuffers.resize(imageCount);
    VkCommandBufferAllocateInfo cmdAllocInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cmdAllocInfo.commandPool = ctx.commandPool;
    cmdAllocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cmdAllocInfo.commandBufferCount = imageCount;
    if (vkAllocateCommandBuffers(ctx.device, &cmdAllocInfo, ctx.commandBuffers.data()) != VK_SUCCESS) {
        log_message("Failed to allocate command buffers");
        return false;
    }

    VkSemaphoreCreateInfo semaphoreInfo = {VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    vkCreateSemaphore(ctx.device, &semaphoreInfo, nullptr, &ctx.imageAvailableSemaphore);
    vkCreateSemaphore(ctx.device, &semaphoreInfo, nullptr, &ctx.renderFinishedSemaphore);
    VkFenceCreateInfo fenceInfo = {VK_STRUCTURE_TYPE_FENCE_CREATE_INFO, nullptr, VK_FENCE_CREATE_SIGNALED_BIT};
    vkCreateFence(ctx.device, &fenceInfo, nullptr, &ctx.inFlightFence);

    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.5f;
    ImGui_ImplSDL3_InitForVulkan(window);
    ImGui_ImplVulkan_InitInfo init_info = {};
    init_info.Instance = ctx.instance;
    init_info.PhysicalDevice = ctx.physicalDevice;
    init_info.Device = ctx.device;
    init_info.Queue = ctx.graphicsQueue;
    init_info.DescriptorPool = ctx.descriptorPool;
    init_info.MinImageCount = imageCount;
    init_info.ImageCount = imageCount;
    ImGui_ImplVulkan_Init(&init_info, ctx.renderPass);

    return true;
}

void render(const Simulation& sim, VulkanContext& ctx, Menu& menu) {
    vkWaitForFences(ctx.device, 1, &ctx.inFlightFence, VK_TRUE, UINT64_MAX);
    vkResetFences(ctx.device, 1, &ctx.inFlightFence);

    uint32_t imageIndex;
    VkResult result = vkAcquireNextImageKHR(ctx.device, ctx.swapchain, UINT64_MAX, ctx.imageAvailableSemaphore, VK_NULL_HANDLE, &imageIndex);
    if (result != VK_SUCCESS) {
        log_message("Failed to acquire next image");
        return;
    }

    struct Uniforms {
        glm::mat4 view_projection;
        float vis_scale, vis_color_intensity;
        int display_mode;
    } uniforms;
    uniforms.vis_scale = sim.equations.empty() ? 1.0f : dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())->getVisScale();
    uniforms.vis_color_intensity = sim.equations.empty() ? 1.0f : dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation].get())->getVisColorIntensity();
    uniforms.vis_scale *= camera_zoom; // Apply camera zoom
    uniforms.display_mode = current_display_mode;
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.01f, 100.0f);
    float cam_dist = 5.0f / camera_zoom;
    glm::vec3 eye(cam_dist * glm::cos(camera_angle), cam_dist * glm::sin(camera_angle), cam_dist * std::sin(camera_tilt));
    glm::mat4 view = glm::lookAt(eye, glm::vec3(0), glm::vec3(0, 0, 1));
    uniforms.view_projection = projection * view;
    void* data;
    vmaMapMemory(ctx.allocator, ctx.uniformAlloc, &data);
    memcpy(data, &uniforms, sizeof(Uniforms));
    vmaUnmapMemory(ctx.allocator, ctx.uniformAlloc);

    double max_scalar = 1e-10, min_scalar = 1e10, max_4d = 1e-10, min_4d = 1e10;
    for (int i = 0; i < sim.size; ++i)
        for (int j = 0; j < sim.size; ++j)
            for (int k = 0; k < sim.size; ++k) {
                max_scalar = std::max(max_scalar, (double)sim.computed_scalar[i][j][k]);
                min_scalar = std::min(min_scalar, (double)sim.computed_scalar[i][j][k]);
                for (int w = 0; w < sim.size; ++w) {
                    max_4d = std::max(max_4d, (double)sim.scalar_4d[i][j][k][w]);
                    min_4d = std::min(min_4d, (double)sim.scalar_4d[i][j][k][w]);
                }
            }
    double range = max_scalar - min_scalar ? max_scalar - min_scalar : 1.0;
    double range_4d = max_4d - min_4d ? max_4d - min_4d : 1.0;
    float iso_level = min_scalar + range * (0.3 + 0.2 * std::sin(sim.t));
    int w_slice = static_cast<int>(sim.t * 2) % sim.size;

    struct PushConstants {
        float iso_level, vis_scale, min_scalar, range, min_4d, range_4d;
        int size, w_slice;
    } pushConstants = {iso_level, uniforms.vis_scale, (float)min_scalar, (float)range, (float)min_4d, (float)range_4d, sim.size, w_slice};

    std::vector<float> scalar_data(sim.size * sim.size * sim.size);
    for (int i = 0; i < sim.size; ++i)
        for (int j = 0; j < sim.size; ++j)
            for (int k = 0; k < sim.size; ++k)
                scalar_data[i * sim.size * sim.size + j * sim.size + k] = sim.computed_scalar[i][j][k];
    vmaMapMemory(ctx.allocator, ctx.scalarAlloc, &data);
    memcpy(data, scalar_data.data(), scalar_data.size() * sizeof(float));
    vmaUnmapMemory(ctx.allocator, ctx.scalarAlloc);
    std::vector<float> scalar4d_data(sim.size * sim.size * sim.size * sim.size);
    for (int i = 0; i < sim.size; ++i)
        for (int j = 0; j < sim.size; ++j)
            for (int k = 0; k < sim.size; ++k)
                for (int w = 0; w < sim.size; ++w)
                    scalar4d_data[i * sim.size * sim.size * sim.size + j * sim.size * sim.size + k * sim.size + w] = sim.scalar_4d[i][j][k][w];
    vmaMapMemory(ctx.allocator, ctx.scalar4DAlloc, &data);
    memcpy(data, scalar4d_data.data(), scalar4d_data.size() * sizeof(float));
    vmaUnmapMemory(ctx.allocator, ctx.scalar4DAlloc);

    std::vector<float> particle_data;
    if (current_display_mode == PARTICLES || current_display_mode == HYBRID) {
        particle_data.reserve(sim.particles.size() * 6);
        for (const auto& p : sim.particles) {
            float speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
            float r = std::min(speed / 5.0f, 1.0f) * uniforms.vis_color_intensity;
            particle_data.insert(particle_data.end(), {
                p.x, p.y, p.z, r, 0.5f * uniforms.vis_color_intensity, (1.0f - r) * uniforms.vis_color_intensity
            });
        }
        vmaMapMemory(ctx.allocator, ctx.particleAlloc, &data);
        memcpy(data, particle_data.data(), particle_data.size() * sizeof(float));
        vmaUnmapMemory(ctx.allocator, ctx.particleAlloc);
    }

    VkCommandBufferBeginInfo beginInfo = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    vkBeginCommandBuffer(ctx.commandBuffers[imageIndex], &beginInfo);
    VkRenderPassBeginInfo renderPassInfo = {VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    renderPassInfo.renderPass = ctx.renderPass;
    renderPassInfo.framebuffer = ctx.framebuffers[imageIndex];
    renderPassInfo.renderArea.extent = {800, 600};
    VkClearValue clearColor = {{{0.0f, 0.0f, 0.0f, 1.0f}}};
    renderPassInfo.clearValueCount = 1;
    renderPassInfo.pClearValues = &clearColor;
    vkCmdBeginRenderPass(ctx.commandBuffers[imageIndex], &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

    // Compute pass
    if (current_display_mode == ISOSURFACE || current_display_mode == HYBRID) {
        vkCmdBindPipeline(ctx.commandBuffers[imageIndex], VK_PIPELINE_BIND_POINT_COMPUTE, ctx.computePipeline);
        vkCmdBindDescriptorSets(ctx.commandBuffers[imageIndex], VK_PIPELINE_BIND_POINT_COMPUTE, ctx.pipelineLayout, 0, 1, &ctx.descriptorSets[imageIndex], 0, nullptr);
        vkCmdPushConstants(ctx.commandBuffers[imageIndex], ctx.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants), &pushConstants);
        vkCmdDispatch(ctx.commandBuffers[imageIndex], (sim.size + 7) / 8, (sim.size + 7) / 8, (sim.size + 7) / 8);
        // Add memory barrier to ensure compute shader completes before vertex fetch
        VkMemoryBarrier barrier = {VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_VERTEX_ATTRIBUTE_READ_BIT;
        vkCmdPipelineBarrier(ctx.commandBuffers[imageIndex], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_VERTEX_INPUT_BIT, 0, 1, &barrier, 0, nullptr, 0, nullptr);
    }

    // Graphics pass
    vkCmdBindPipeline(ctx.commandBuffers[imageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, ctx.graphicsPipeline);
    vkCmdBindDescriptorSets(ctx.commandBuffers[imageIndex], VK_PIPELINE_BIND_POINT_GRAPHICS, ctx.pipelineLayout, 0, 1, &ctx.descriptorSets[imageIndex], 0, nullptr);
    VkDeviceSize offsets[] = {0};
    if (current_display_mode == PARTICLES || current_display_mode == HYBRID) {
        vkCmdBindVertexBuffers(ctx.commandBuffers[imageIndex], 0, 1, &ctx.particleBuffer, offsets);
        vkCmdDraw(ctx.commandBuffers[imageIndex], sim.particles.size(), 1, 0, 0);
    }
    if (current_display_mode == ISOSURFACE || current_display_mode == HYBRID) {
        vkCmdBindVertexBuffers(ctx.commandBuffers[imageIndex], 0, 1, &ctx.vertexBuffer, offsets);
        vkCmdDraw(ctx.commandBuffers[imageIndex], sim.size * sim.size * sim.size * 15, 1, 0, 0);
    }

    // ImGui rendering
    ImGui_ImplVulkan_NewFrame();
    ImGui_ImplSDL3_NewFrame();
    ImGui::NewFrame();
    if (!sim.equations.empty()) {
        ImGui::Begin("Info");
        ImGui::SetWindowPos(ImVec2(10, 10));
        ImGui::Text("%s | t = %.2f", sim.equations[sim.current_equation]->name().c_str(), sim.t);
        ImGui::End();
    }
    menu.render(sim);
    ImGui::Render();
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
    if (vkQueueSubmit(ctx.graphicsQueue, 1, &submitInfo, ctx.inFlightFence) != VK_SUCCESS) {
        log_message("Failed to submit draw command buffer");
        return;
    }
    VkPresentInfoKHR presentInfo = {VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    presentInfo.waitSemaphoreCount = 1;
    presentInfo.pWaitSemaphores = signalSemaphores;
    presentInfo.swapchainCount = 1;
    presentInfo.pSwapchains = &ctx.swapchain;
    presentInfo.pImageIndices = &imageIndex;
    vkQueuePresentKHR(ctx.presentQueue, &presentInfo);
}

void cleanup_renderer(VulkanContext& ctx) {
    vkDeviceWaitIdle(ctx.device);
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplSDL3_Shutdown();
    ImGui::DestroyContext();
    vmaDestroyBuffer(ctx.allocator, ctx.uniformBuffer, ctx.uniformAlloc);
    vmaDestroyBuffer(ctx.allocator, ctx.scalarBuffer, ctx.scalarAlloc);
    vmaDestroyBuffer(ctx.allocator, ctx.scalar4DBuffer, ctx.scalar4DAlloc);
    vmaDestroyBuffer(ctx.allocator, ctx.vertexBuffer, ctx.vertexAlloc);
    vmaDestroyBuffer(ctx.allocator, ctx.particleBuffer, ctx.particleAlloc);
    vmaDestroyAllocator(ctx.allocator);
    vkDestroySemaphore(ctx.device, ctx.imageAvailableSemaphore, nullptr);
    vkDestroySemaphore(ctx.device, ctx.renderFinishedSemaphore, nullptr);
    vkDestroyFence(ctx.device, ctx.inFlightFence, nullptr);
    vkDestroyCommandPool(ctx.device, ctx.commandPool, nullptr);
    vkDestroyPipeline(ctx.device, ctx.graphicsPipeline, nullptr);
    vkDestroyPipeline(ctx.device, ctx.computePipeline, nullptr);
    vkDestroyPipelineLayout(ctx.device, ctx.pipelineLayout, nullptr);
    vkDestroyDescriptorSetLayout(ctx.device, ctx.descriptorSetLayout, nullptr);
    vkDestroyDescriptorPool(ctx.device, ctx.descriptorPool, nullptr);
    for (auto framebuffer : ctx.framebuffers) vkDestroyFramebuffer(ctx.device, framebuffer, nullptr);
    for (auto imageView : ctx.swapchainImageViews) vkDestroyImageView(ctx.device, imageView, nullptr);
    vkDestroySwapchainKHR(ctx.device, ctx.swapchain, nullptr);
    vkDestroyRenderPass(ctx.device, ctx.renderPass, nullptr);
    vkDestroySurfaceKHR(ctx.instance, ctx.surface, nullptr);
    vkDestroyDevice(ctx.device, nullptr);
    vkDestroyInstance(ctx.instance, nullptr);
}