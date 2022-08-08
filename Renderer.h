//
// Created by Artem Novikov on 06.08.2022.
//

#ifndef GEOMETRY__RENDERER_H_
#define GEOMETRY__RENDERER_H_

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <OpenGL/gl.h>

#include "GeometryEntities/Transform.h"
#include "Geometry2D.h"

struct Style {
  Vector3f entity_color = {1.f, 1.f, 1.f};
  bool hidden = false;
};

struct Vertex2f {
  Point2f position;
  Vector3f color;
};

static const float kUnitLength = 50.f; /// size of a unit in pixels

class Renderer {
 public:
  Renderer(int window_width, int window_height) /// in pixels
      : window_box_(-Vector2f(window_width * 0.5f, window_height * 0.5f) / kUnitLength,
                    Vector2f(window_width * 0.5f, window_height * 0.5f) / kUnitLength),
        fit_to_screen_({{2.f / window_box_.GetSide(0), 0.f}, {0.f, 2.f / window_box_.GetSide(1)}}, {0.f, 0.f}) {}

  ~Renderer() = default;

  void RegisterWindowShift(const Vector2f& shift) {
    window_box_.ApplyTransform(Translate(shift.Scaled(scale_) / kUnitLength));
    RecalculateTransform();
  }

  void RegisterWindowScale(const Vector2f& new_window_size) {
    auto pivot = window_box_.GetLeft();
    pivot[1] = window_box_.GetRight()[1];
    window_box_.GetRight() = pivot + Vector2f(new_window_size[0] * scale_[0], 0) / kUnitLength;
    window_box_.GetLeft() = pivot + Vector2f(0, -new_window_size[1] * scale_[1]) / kUnitLength;
    RecalculateTransform();
  }

  void RegisterCoordinateSystemShift(const Vector2f& shift) {
    window_box_.ApplyTransform(Translate(shift.Scaled(scale_) / kUnitLength)); /// TODO fix if doesnt work
    RecalculateTransform();
  }

  void RegisterCoordinateSystemScale(const Vector2f& scale) {
    auto center = window_box_.GetCenter();
    window_box_.ApplyTransform(Scale(scale, center));
    scale_ *= scale;
    RecalculateTransform();
  }

  void Render(const std::shared_ptr<Void2f>& object, const Style& style) const {
    switch (object->GetType()) {
      case Entity::Point: {
        Render(static_cast<const Point2f&>(*object), style);
        break;
      }
      case Entity::Line: {
        Render(static_cast<const Line2f&>(*object), style);
        break;
      }
      case Entity::Segment: {
        Render(static_cast<const Segment2f&>(*object), style);
        break;
      }
      case Entity::Bezier: {
        /// TODO
        break;
      }
    }
  }

  void RenderCircle(
      const Point2f& center, float radius, const Vector3f& color = {1.f, 1.f, 1.f}, size_t vertex_count = 30) const {
    glBegin(GL_POLYGON);
    float delta = 2 * M_PI / vertex_count;
    for (size_t i = 0; i < vertex_count; ++i) {
      auto position = center + Vector2f(std::cosf(i * delta), std::sinf(i * delta)) * radius;
      RenderVertex(position, color);
    }
    glEnd();
  }

  void RenderCoordinateSystem() {
    RenderCircle({0, 0}, 1.f);
  }

  void Render(const Point2f& point, const Style& style) const {
    if (!style.hidden) {
      RenderCircle(point, 0.1f * GetScale(), style.entity_color); /// TODO scales the dimensions of the point
    }
  }

  void Render(const Line2f& line, const Style& style) const {
    auto intersections = line.Intersection(window_box_);
    if (!style.hidden && intersections.has_value()) {
      glBegin(GL_LINES);
      RenderVertex(intersections.value().first, style.entity_color);
      RenderVertex(intersections.value().second, style.entity_color);
      glEnd();
    }
  }

  void Render(const Segment2f& segment, const Style& style) const {
    if (!style.hidden
        && (window_box_.Contains(segment.GetLeft()) || window_box_.Contains(segment.GetRight()))) {
      glBegin(GL_LINES);
      RenderVertex(segment.GetLeft(), style.entity_color);
      RenderVertex(segment.GetRight(), style.entity_color);
      glEnd();
    }
  }

  Point2f MapCursorToGeometry(const Point2f& click_position) const {
    return Inverted(fit_to_screen_)(click_position);
  }

  Point2f MapCursorToWindow(const Point2f& geometry_position) const {
    return fit_to_screen_(geometry_position);
  }

  float GetScale() const { return scale_[0]; } /// TODO Dimensions

 private:
  void RenderVertex(const Vertex2f& vertex) const {
    auto screen_position = fit_to_screen_(vertex.position);
    glVertex2f(screen_position[0], screen_position[1]);
    glColor3f(vertex.color[0], vertex.color[1], vertex.color[2]);
  }

  void RenderVertex(const Point2f& position, const Vector3f& color) const {
    auto screen_position = fit_to_screen_(position);
    glVertex2f(screen_position[0], screen_position[1]);
    glColor3f(color[0], color[1], color[2]);
  }

  void RecalculateTransform() {
    fit_to_screen_.GetMatrix()[0][0] = 2.f / window_box_.GetSide(0);
    fit_to_screen_.GetMatrix()[1][1] = 2.f / window_box_.GetSide(1);
    fit_to_screen_.GetShift() = -(fit_to_screen_.GetMatrix() * window_box_.GetCenter());
  }

  Vector2f scale_ = {1.f, 1.f}; /// scale of a coordinate system
  BoundaryBox2f window_box_; /// window in a coordinate system
  Transform2f fit_to_screen_; /// transforms coordinate system to that of OpenGL
};

#endif //GEOMETRY__RENDERER_H_
