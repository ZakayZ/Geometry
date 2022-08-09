//
// Created by Artem Novikov on 04.08.2022.
//

#ifndef GEOMETRY__UI_H_
#define GEOMETRY__UI_H_

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <cstdio>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include "Tools.h"
#include "UserInputs.h"

static float kToolSize = 50.f;

template <size_t ToolCount>
class Section {
 public:
  Section(const char* section_name,
          const std::array<const char*, ToolCount>& tool_names,
          const std::array<char, ToolCount>& tool_states,
          const std::array<ToolType, ToolCount>& tool_types)
      : section_name_(section_name), tool_names_(tool_names), tool_states_(tool_states), tool_types_(tool_types) {}

  std::array<char, ToolCount>& GetStates() { return tool_states_; }

  ToolType GenerateSection(char*& active_state) {
    ToolType got_switch = ToolType::None;
    if (ImGui::CollapsingHeader(section_name_)) {
      ImVec2 window_min = ImGui::GetWindowContentRegionMin();
      ImVec2 window_max = ImGui::GetWindowContentRegionMax();
      float width = window_max.x - window_min.x;
      size_t columns_size = std::max(1, std::min(int(width / kToolSize), int(ToolCount)));
      for (size_t tool = 0; tool < ToolCount; ++tool) {
        if (tool % columns_size > 0) { ImGui::SameLine(); }
        ImGui::PushID(tool);
        if (ImGui::Selectable(tool_names_[tool], tool_states_[tool] != 0, 0, ImVec2(kToolSize, kToolSize))) {
          if (*active_state) {
            *active_state = false;
          }
          active_state = &tool_states_[tool];
          *active_state = true;
          got_switch = tool_types_[tool];
        }
        ImGui::PopID();
      }
    }
    return got_switch;
  }

  constexpr static size_t GetToolCount() {
    return ToolCount;
  }

 private:
  const char* section_name_;
  std::array<const char*, ToolCount> tool_names_;
  std::array<char, ToolCount> tool_states_;
  std::array<ToolType, ToolCount> tool_types_;
};

struct ToolSections {
  Section<3> inspection = Section<3>("Inspection",
                                     {"move", "navigate", "delete"},
                                     {false, false, false},
                                     {ToolType::Inspect_Move, ToolType::Inspect_Navigate, ToolType::Inspect_Delete});
  Section<2> construction = Section<2>("Construction",
                                       {"midpoint", "project"},
                                       {false, false},
                                       {ToolType::Construct_Midpoint, ToolType::Construct_Project});
  Section<3> creation = Section<3>("Creation",
                                   {"point", "line", "segment"},
                                   {false, false, false},
                                   {ToolType::Create_Point, ToolType::Create_Line, ToolType::Create_Segment});
};

class UI {
 public:
  UI(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    io_ = &ImGui::GetIO();
    SetupIO();

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    if (io_->ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
      style.WindowRounding = 0.0f;
      style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    active_state_ = &sections_.inspection.GetStates()[0];
    *active_state_ = true;
  }

  ~UI() {
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
  }

  void SetUserInputs(UserInputs& user_inputs) {
    user_inputs_ = &user_inputs;
  }

  void StartFrame() {
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    MakeUI();
  }

  void Render() {
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
  }

  void EndFrame() {
    if (io_->ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
      GLFWwindow* backup_current_context = glfwGetCurrentContext();
      ImGui::UpdatePlatformWindows();
      ImGui::RenderPlatformWindowsDefault();
      glfwMakeContextCurrent(backup_current_context);
    }
  }

  ImGuiIO& GetIO() { return *io_; }

 private:
  void SetupIO() {
    io_->ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;       // Enable Keyboard Controls
    //io_->ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    io_->ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // Enable Docking
    io_->ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;         // Enable Multi-Viewport / Platform Windows
    //io_->ConfigViewportsNoAutoMerge = true;
    //io_->ConfigViewportsNoTaskBarIcon = true;
  }

  void MakeUI() {
    MakeDemo();
    ImGui::Begin("Geometry Tools");

    ToolType new_type = ToolType::None;
    ToolType type;

    type = MakeInspectSection();
    new_type = type != ToolType::None ? type : new_type;

    type = MakeConstructSection();
    new_type = type != ToolType::None ? type : new_type;

    type = MakeCreateSection();
    new_type = type != ToolType::None ? type : new_type;

    user_inputs_->SetToolType(new_type);

    ImGui::End();
  }

  void MakeDemo() {
    if (show_demo_window) {
      ImGui::ShowDemoWindow(&show_demo_window);
    }
  }

  ToolType MakeInspectSection() { return sections_.inspection.GenerateSection(active_state_); }

  ToolType MakeConstructSection() { return sections_.construction.GenerateSection(active_state_); }

  ToolType MakeCreateSection() { return sections_.creation.GenerateSection(active_state_); }

  bool show_demo_window = true;
  ImGuiIO* io_ = nullptr;
  ToolSections sections_;
  char* active_state_ = nullptr;
  UserInputs* user_inputs_ = nullptr;
};

#endif //GEOMETRY__UI_H_
