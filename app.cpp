//
// Created by Artem Novikov on 08.08.2022.
//

#include "GeometryApp.h"

static Vector2i window_size = {1000, 1000};

int main() {
  GeometryApp app(window_size[0], window_size[1]);
  app.RunApp();
}

