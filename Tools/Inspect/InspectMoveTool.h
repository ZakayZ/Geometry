//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_
#define GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_

#include "Tool.h"

class InspectMoveTool : public Tool {
 public:
  InspectMoveTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  virtual ~InspectMoveTool() = default;

  virtual void ProcessPressed(const Point2f& clicked_pos) {
    auto selected_objects = output_geometry_.Selected(clicked_pos);
    if (selected_objects.empty()) { return; }
  }

  virtual Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) { return cursor_pos; }

  virtual void ProcessDown(const Point2f& cursor_pos) {

  }

  virtual void ProcessReleased(const Point2f& cursor_pos) {
    ProcessDown(cursor_pos);
    Refresh();
  }

  virtual void Refresh() {
    object_to_move_ = nullptr;
    connected_objects_.clear();
  }

 private:
  std::shared_ptr<Void2f> object_to_move_ = nullptr;
  std::list<std::shared_ptr<Void2f>> connected_objects_;
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_
