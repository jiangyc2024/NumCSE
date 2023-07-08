///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2017 SAM, D-MATH
/// Author(s): Till Ehrengruber <tille@student.ethz.ch>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////
#ifndef UTILS_VOID_CB
#define UTILS_VOID_CB

#include <cassert>

/**
 * Functor behaving like an unset function pointer (nullptr) with
 * arbitrary number of arguments.
 *
 * Usage:
 * ```
 * #include "Utils/void_cb.hpp"
 *
 * template <typename CB=void_cb>
 * void f(CB callback=nullptr) {
 *   // ...
 *   if (callback != nullptr) {
 *     callback("whatever arguments your callback takes");
 *   }
 *   // ...
 * }
 * ```
 */
struct void_cb {
  explicit constexpr void_cb(void* dummy = nullptr) {}

  template <typename... T>
  void operator()(T... /*unused*/) {
    assert(false);
  }

  explicit constexpr operator void*() const { return nullptr; }
};
#endif
