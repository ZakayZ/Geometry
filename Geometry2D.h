//
// Created by Artem Novikov on 04.08.2022.
//

#ifndef GEOMETRY__GEOMETRY2D_H_
#define GEOMETRY__GEOMETRY2D_H_

#include <list>
#include <memory>

#include "GeometryEntities/Void.h"
#include "GeometryEntities/Point.h"
#include "GeometryEntities/Segment.h"
#include "GeometryEntities/Line.h"

template <typename T>
class Geometry2D {
 public:
  using iterator = typename std::list<std::shared_ptr<Void2<T>>>::iterator;
  using const_iterator = typename std::list<std::shared_ptr<Void2<T>>>::const_iterator;

  Geometry2D() = default;
  ~Geometry2D() = default;

  iterator begin() { return geometric_objects_.begin(); }
  const_iterator begin() const { return geometric_objects_.cbegin(); }
  const_iterator cbegin() const { return geometric_objects_.cbegin(); }
  iterator end() { return geometric_objects_.end(); }
  const_iterator end() const { return geometric_objects_.cend(); }
  const_iterator cend() const { return geometric_objects_.cend(); }

  auto GetData() const { return geometric_objects_; }

  template <typename U>
  void Push(const U& object) {
    geometric_objects_.push_back(std::make_shared<U>(object));
  }

  void Erase(const std::shared_ptr<Void2<T>> object_ptr) {
    for (auto it = geometric_objects_.begin(); it != geometric_objects_.end(); ++it) {
      if (object_ptr == *it) {
        geometric_objects_.erase(it);
        break;
      }
    }
  }

  std::list<std::shared_ptr<Void2<T>>> Selected(const Point2<T>& selection_position,
                                                const T& vicinity_distance = 0) const {
    std::list<std::shared_ptr<Void2<T>>> found;
    std::function<bool(const std::shared_ptr<Void2<T>>&)> is_close;
    if (vicinity_distance == 0) {
      is_close = [&selection_position](const std::shared_ptr<Void2<
          T>>& object_ptr) {
        return object_ptr->Contains(selection_position);
      };
    } else {
      is_close = [&vicinity_distance, &selection_position](const std::shared_ptr<
          Void2<T>>& object_ptr) {
        return object_ptr->SquaredDistance(selection_position) <= vicinity_distance;
      };
    }

    for (const auto& object_ptr : geometric_objects_) {
      if (is_close(object_ptr)) {
        found.push_back(object_ptr);
      }
    }
    return found;
  }

 private:
  std::list<std::shared_ptr<Void2<T>>> geometric_objects_;
};

#endif //GEOMETRY__GEOMETRY2D_H_
