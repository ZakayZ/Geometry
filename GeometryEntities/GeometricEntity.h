//
// Created by Artem Novikov on 05.07.2022.
//

#include <cstddef>
#include <iostream>
#include "Entities.h"

#ifndef GEOMERTY_GEOMETRYENTITIES_GEOMETRICENTITY_H_
#define GEOMERTY_GEOMETRYENTITIES_GEOMETRICENTITY_H_

class GeometryEntity {
 public:
  GeometryEntity() : data_(nullptr) {}

  template <typename T>
  explicit GeometryEntity(const T& data) : data_(new Derived<T>(data)) {}

  template <typename T>
  GeometryEntity(T&& data) : data_(new Derived<T>(data)) {}

  GeometryEntity(const GeometryEntity& other) = delete;

  GeometryEntity(GeometryEntity&& other) noexcept: data_(other.data_) { other.data_ = nullptr; }

  GeometryEntity& operator=(const GeometryEntity& other) = delete;

  GeometryEntity& operator=(GeometryEntity&& other) noexcept {
    if (this != &other) {
      delete data_;
      data_ = other.data_;
      other.data_ = nullptr;
    }
    return *this;
  }

  ~GeometryEntity() { delete data_; }

  template <typename Type>
  Type& GetValue() { return dynamic_cast<Derived<Type>*>(data_)->GetValue(); }
  template <typename Type>
  const Type& GetValue() const { return dynamic_cast<Derived<Type>*>(data_)->GetValue(); }
  Entity GetType() const { return data_->GetType(); }
  size_t GetDimension() const { return data_->GetDimension(); }

 private:
  class Base {
   public:
    Base() = default;
    virtual Entity GetType() const = 0;
    virtual size_t GetDimension() const = 0;
    virtual ~Base() = default;
  };

  template <typename Type>
  class Derived : public Base {
   public:
    template <typename... Args>
    Derived(Args&& ... args) : value_(std::forward<Args>(args)...) {}
    Type& GetValue() { return value_; }
    const Type& GetValue() const { return value_; }
    Entity GetType() const override { return value_.GetType(); }
    size_t GetDimension() const override { return value_.GetDimension(); }
    ~Derived() override = default;
   private:
    Type value_;
  };

  template <typename T, typename... Args>
  friend GeometryEntity MakeGeometryEntity(Args&& ... args);

  Base* data_;
};

template <typename T, typename... Args>
GeometryEntity MakeGeometryEntity(Args&& ... args) {
  GeometryEntity geometry_entity;
  geometry_entity.data_ = new GeometryEntity::Derived<T>(std::forward<Args>(args)...);
  return geometry_entity;
}

#endif //GEOMERTY_GEOMETRYENTITIES_GEOMETRICENTITY_H_
