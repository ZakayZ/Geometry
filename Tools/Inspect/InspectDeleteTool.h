//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_
#define GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_

#include "Tool.h"

class InspectDeleteTool : public Tool { /// TODO
 public:
  InspectDeleteTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}
 private:
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_
