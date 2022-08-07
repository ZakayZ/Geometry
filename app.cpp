//
// Created by Artem Novikov on 08.08.2022.
//

#include "GeometryApp.h"

static Vector2i window_size = {1000, 1000};

int main() {
  AppManager app_manager(window_size[0], window_size[1]);
  auto& app = app_manager.GetApp();
  app.RunApp();
  app_manager.DestroyApp();
}

