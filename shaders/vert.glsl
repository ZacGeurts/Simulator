#version 450
layout(location = 0) in vec3 pos;
layout(location = 1) in float color_intensity;
layout(location = 0) out vec3 fragColor;
layout(binding = 0) uniform UBO {
    mat4 view_projection;
    float vis_scale, vis_color_intensity;
    int display_mode;
} ubo;
void main() {
    gl_Position = ubo.view_projection * vec4(pos * ubo.vis_scale, 1.0);
    fragColor = vec3(color_intensity * ubo.vis_color_intensity, 0.5 * ubo.vis_color_intensity, (1.0 - color_intensity) * ubo.vis_color_intensity);
}
