//
// Created by Artem Novikov on 02.07.2022.
//

#include <type_traits>

#pragma once

template <class Container>
struct is_iterable_helper {
 public:
  template <class Cont>
  static auto f(int) -> decltype(std::declval<typename Cont::iterator>(), std::true_type());

  template <class...>
  static auto f(...) -> std::false_type;
 public:
  static const bool value = decltype(f<Container>(0))::value;
};

template <class Container>
static const bool is_iterable_v = is_iterable_helper<Container>::value;
