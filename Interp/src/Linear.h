#ifndef INTERP_LINEAR_GUARD_H
#define INTERP_LINEAR_GUARD_H

#include <algorithm>
#include <concepts>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interp {

template <typename xIter, typename yIter>
requires requires(xIter x, yIter y) {
  requires RealFloatingPointIterator<xIter>;
  requires FloatingPointIterator<yIter>;
  { (*x) + (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
  { (*x) * (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
}
class Linear {
 public:
  // Define some class member types
  using x_value_type = std::iter_value_t<xIter>;
  using y_value_type = std::iter_value_t<yIter>;

  // Declare the constructor.
  Linear(xIter, xIter, yIter);

  // Declare the application operator.
  y_value_type operator()(x_value_type) const;

 private:
  // Iterators to the function data.
  xIter xStart;
  xIter xFinish;
  yIter yStart;
};

template <typename xIter, typename yIter>
Linear<xIter, yIter>::Linear(xIter xStart, xIter xFinish, yIter yStart)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {}

template <typename xIter, typename yIter>
Linear<xIter, yIter>::y_value_type Linear<xIter, yIter>::operator()(
    const x_value_type x) const {
  auto iter = std::upper_bound(xStart, xFinish, x);
  if (iter == xStart) {
    return *yStart;
  }
  auto i2 = std::distance(xStart, iter);
  auto i1 = i2 - 1;
  if (iter == xFinish) {
    return *std::next(yStart, i1);
  } else {
    auto x1 = *std::prev(iter);
    auto x2 = *iter;
    auto y1 = *std::next(yStart, i1);
    auto y2 = *std::next(yStart, i2);
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
  }
}

}  // namespace Interp

#endif  //  INTERP_LINEAR_GUARD_H
