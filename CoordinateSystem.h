//
// Created by Artem Novikov on 08.08.2022.
//

#ifndef GEOMETRY__COORDINATESYSTEM_H_
#define GEOMETRY__COORDINATESYSTEM_H_

#include "GeometryEntities/Point.h"

class CoordinateSystem {
 public:
  CoordinateSystem() = default;

  std::vector<Point2f> GenerateMajorLines(const Vector2f& left_corner, const Vector2f& right_corner) const {
    return GenerateLines(left_corner, right_corner);
  }

  std::vector<Point2f> GenerateMinorLines(const Vector2f& left_corner, const Vector2f& right_corner) const {
    return GenerateLines(left_corner, right_corner, 0.2f);
  }

  const Vector3f& GetMajorColor() const {
    return major_color_;
  }

  const Vector3f& GetMinorColor() const {
    return minor_color_;
  }

 private:
  static float ProperTrunc(float input) {
    return input >= 0.f ? std::truncf(input) : std::truncf(input) - 1.f;
  }

  static float CalculateStep(float size) {
    return std::powf(10.f, ProperTrunc(std::log10(size) - 1.f));
  }

  static Vector2f GetModifier(const Vector2i& number_lines) {
    auto one_dimension_modifier = [](int number_lines) {
      if (number_lines > 50) { return 5.f; }
      if (number_lines > 20) { return 2.f; }
      return 1.f;
    };
    return Vector2f(one_dimension_modifier(number_lines[0]), one_dimension_modifier(number_lines[1]));
  }

  static Vector2f GetStep(const Vector2f& visible_size) {
    return Vector2f(CalculateStep(visible_size[0]), CalculateStep(visible_size[1]));
  }

  static std::vector<Point2f> GenerateLines(
      const Vector2f& left_corner, const Vector2f& right_corner, float modifier = 1.f) {
    std::vector<Point2f> grid_lines_points;
    Vector2f grid_step = GetStep(right_corner - left_corner);
    Vector2i number_visible_lines = (right_corner - left_corner).InvertedScale(grid_step);
    grid_step *= GetModifier(number_visible_lines) * modifier;
    auto left_rounded_corner = GetRoundCorner(left_corner, grid_step);
    number_visible_lines = (right_corner - left_rounded_corner).InvertedScale(grid_step);
    grid_lines_points.reserve(2 * (number_visible_lines[0] + number_visible_lines[1] + 2));

    for (int x_num = 0; x_num < number_visible_lines[0] + 1; ++x_num) {
      grid_lines_points.emplace_back(left_rounded_corner[0] + x_num * grid_step[0], left_corner[1]);
      grid_lines_points.emplace_back(left_rounded_corner[0] + x_num * grid_step[0], right_corner[1]);
    }

    for (int y_num = 0; y_num < number_visible_lines[1] + 1; ++y_num) {
      grid_lines_points.emplace_back(left_corner[0], left_rounded_corner[1] + y_num * grid_step[1]);
      grid_lines_points.emplace_back(right_corner[0], left_rounded_corner[1] + y_num * grid_step[1]);
    }

    return grid_lines_points;
  }

  static Vector2f GetRoundCorner(const Vector2f& left_corner, const Vector2f& delta) {
    Vector2f rounded_corner = left_corner;
    rounded_corner.InvertedScale(delta);
    rounded_corner = Vector2f(ProperTrunc(rounded_corner[0]), ProperTrunc(rounded_corner[1]));
    rounded_corner.Scale(delta);
    return rounded_corner;
  }

  Vector3f major_color_ = {0.1f, 0.1f, 0.1f};
  Vector3f minor_color_ = {0.4f, 0.4f, 0.4f};
};

#endif //GEOMETRY__COORDINATESYSTEM_H_
