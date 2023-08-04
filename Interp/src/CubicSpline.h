#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <Eigen/Core>
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
  using xVector = Eigen::Matrix<x_value_type, Eigen::Dynamic, 1>;
  using yVector = Eigen::Matrix<y_value_type, Eigen::Dynamic, 1>;

  // Store iterators to the function data.
  xIter xStart, xFinish;
  yIter yStart;
};

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter XFinish,
                                       yIter yStart)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
  std::cout << "Hello\n";
}

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::y_value_type CubicSpline<xIter, yIter>::operator()(
    x_value_type x) {
  return 0;
}

}  // namespace Interp

#endif  //  INTERP_CUBIC_SPLINE_GUARD_H
