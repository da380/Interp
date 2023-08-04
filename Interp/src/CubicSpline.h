#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <vector>

namespace Interp {

template <typename xIter, typename yIter>
class CubicSpline {
 public:
  // Define some class member types
  using x_value_type = xIter::value_type;
  using y_value_type = yIter::value_type;

  // Declare the constructor.
  CubicSpline(xIter, xIter, yIter);

  // Declare the application operator.
  y_value_type operator()(x_value_type);

 private:
  using Vector = Eigen::Matrix<y_value_type, Eigen::Dynamic, 1>;

  // Iterators to the function data.
  xIter xStart;
  xIter xFinish;
  yIter yStart;

  // Cubic spline coefficients
  Vector ypp;
};

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {}

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::y_value_type CubicSpline<xIter, yIter>::operator()(
    x_value_type x) {
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

#endif  //  INTERP_CUBIC_SPLINE_GUARD_H
