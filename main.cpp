#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <cstdio>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include <cassert>
#include <array>
//#include "UI.h"
#include "Geometry2D.h"

static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

class Window {
 public:
  Window(size_t width, size_t height) {
    glfwSetErrorCallback(glfw_error_callback);
    assert(glfwInit());
    window_ = glfwCreateWindow(width, height, "Dear ImGui GLFW+OpenGL2 example", nullptr, nullptr);
    assert(window_ != nullptr);
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // Enable vsync

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    io_ = &ImGui::GetIO();
    (void) *io_;
    io_->ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;       // Enable Keyboard Controls
    //io_->ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    io_->ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // Enable Docking
    io_->ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;         // Enable Multi-Viewport / Platform Windows
    //io_->ConfigViewportsNoAutoMerge = true;
    //io_->ConfigViewportsNoTaskBarIcon = true;

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    ImGuiStyle& style = ImGui::GetStyle();
    if (io_->ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
      style.WindowRounding = 0.0f;
      style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL2_Init();
    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);
  }

  void RenderWindow() {
    while (!glfwWindowShouldClose(window_)) {
      // Poll and handle events (inputs, window resize, etc.)
      // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
      // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
      // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
      // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
      glfwPollEvents();

      // Start the Dear ImGui frame
      ImGui_ImplOpenGL2_NewFrame();
      ImGui_ImplGlfw_NewFrame();
      ImGui::NewFrame();

      /// 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
      if (show_demo_window) {
        ImGui::ShowDemoWindow(&show_demo_window);
      }

      /// 2. Show a simple window that we create ourselves. We use a Begin/End pair to created a named window.
      {
        static float f = 0.0f;
        static int counter = 0;
        static char* previous_state = nullptr;

        ImGui::Begin("Hello, world!");

        if (ImGui::CollapsingHeader("Inspect")) {
          static const size_t number_tools = 3;
          static std::array<const char*, number_tools> inspection_names = {"move", "rotate", "delete"};
          static std::array<char, number_tools> inspection_tools = {false, false, false};
          for (size_t x = 0; x < inspection_tools.size(); ++x) { /// TODO variable layout
            if (x > 0) { ImGui::SameLine(); }
            ImGui::PushID(x);
            if (ImGui::Selectable(inspection_names[x], inspection_tools[x] != 0, 0, ImVec2(50, 50))) {
              if (previous_state) {
                *previous_state ^= 1;
              }
              previous_state = &inspection_tools[x];
              *previous_state ^= 1;
            }
            ImGui::PopID();
          }
        }

        if (ImGui::CollapsingHeader("Construct")) {
          static const size_t number_tools = 4;
          static std::array<const char*, number_tools>
              inspection_names = {"midpoint", "bisector", "tangents", "projection"};
          static std::array<char, number_tools> inspection_tools = {false, false, false, false};
          for (size_t x = 0; x < inspection_tools.size(); ++x) { /// TODO variable layout
            if (x > 0) { ImGui::SameLine(); }
            ImGui::PushID(x);
            if (ImGui::Selectable(inspection_names[x], inspection_tools[x] != 0, 0, ImVec2(50, 50))) {
              if (previous_state) {
                *previous_state ^= 1;
              }
              previous_state = &inspection_tools[x];
              *previous_state ^= 1;
            }
            ImGui::PopID();
          }
        }

        if (ImGui::CollapsingHeader("Create")) {
          static const size_t number_tools = 3;
          static std::array<const char*, number_tools>
              inspection_names = {"point", "line", "segment"};
          static std::array<char, number_tools> inspection_tools = {false, false, false};
          for (size_t x = 0; x < inspection_tools.size(); ++x) { /// TODO variable layout
            if (x > 0) { ImGui::SameLine(); }
            ImGui::PushID(x);
            if (ImGui::Selectable(inspection_names[x], inspection_tools[x] != 0, 0, ImVec2(50, 50))) {
              if (previous_state) {
                *previous_state ^= 1;
              }
              previous_state = &inspection_tools[x];
              *previous_state ^= 1;
            }
            ImGui::PopID();
          }
        }

        ImGui::Text("This is some useful text.");
        ImGui::Checkbox("Demo Window", &show_demo_window);
        ImGui::Checkbox("Another Window", &show_another_window);

        ImGui::SliderFloat("float", &f, 0.0f, 1.0f);
        ImGui::ColorEdit3("clear color", (float*) &clear_color);

        if (ImGui::Button("Button")) {
          counter++;
        }
        ImGui::SameLine();
        ImGui::Text("counter = %d", counter);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                    1000.0f / ImGui::GetIO().Framerate,
                    ImGui::GetIO().Framerate);
        ImGui::End();
      }

      /// 3. Show another simple window.
      if (show_another_window) {
        ImGui::Begin("Another Window", &show_another_window);
        ImGui::Text("Hello from another window!");
        if (ImGui::Button("Close Me")) {
          show_another_window = false;
        }
        ImGui::End();
      }

      /// Rendering

      ImGui::Render();
      int display_w, display_h;
      glfwGetFramebufferSize(window_, &display_w, &display_h);
      glViewport(0, 0, display_w, display_h);
      glClearColor(clear_color.x * clear_color.w,
                   clear_color.y * clear_color.w,
                   clear_color.z * clear_color.w,
                   clear_color.w);
      glClear(GL_COLOR_BUFFER_BIT);

      ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

      glBegin(GL_QUADS); // Start drawing a quad primitive
      glColor3f(1.0, 0.0, 0.0);   // red
      glVertex2f(0.0f, 1.0f); // The bottom left corner
      glColor3f(0.f, 0.0f, 1.0f);
      glVertex2f(1.0f, 1.0f); // The top left corner
      glColor3f(1.0, 0.0, 0.0);   // red
      glVertex2f(1.0f, 0); // The top right corner
      glColor3f(0.0, 1.0, 0.0);   // red
      glVertex2f(0, 0); // The bottom right corner

      glEnd();
      // MyApp.Render(); /// TODO

      // Update and Render additional Platform Windows
      // (Platform functions may change the current OpenGL context, so we save/restore it to make it easier to paste this code elsewhere.
      //  For this specific demo app we could also call glfwMakeContextCurrent(window) directly)
      if (io_->ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        GLFWwindow* backup_current_context = glfwGetCurrentContext();
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
        glfwMakeContextCurrent(backup_current_context);
      }

      glfwSwapBuffers(window_);
    }
  }

  ~Window() {
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window_);
    glfwTerminate();
  }

 private:
  GLFWwindow* window_ = nullptr;
  ImGuiIO* io_ = nullptr;

  bool show_demo_window = true;
  bool show_another_window = true;
  ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
};

int main(int, char**) {
  Window window(2000, 1000);
  window.RenderWindow();
}
