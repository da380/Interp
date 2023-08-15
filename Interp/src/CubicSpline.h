#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <algorithm>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interp {

// Enum class for boundary condition types.
enum class CubicSplineBC { Free, Clamped };

template <typename xIter, typename yIter>
requires InterpolationIteratorPair<xIter, yIter>
class CubicSpline {
    public:
     // Define some class member types
     using x_value_t = std::iter_value_t<xIter>;
     using y_value_t = std::iter_value_t<yIter>;

     // General constructor.
     CubicSpline(xIter, xIter, yIter, CubicSplineBC, y_value_t, CubicSplineBC,
                 y_value_t);

     // Constructor for natural spline.
     CubicSpline(xIter, xIter, yIter);

     // Constructor when boundary conditions are the same.
     CubicSpline(xIter, xIter, yIter, CubicSplineBC, y_value_t, y_value_t);

     // Evaluate interpolating function.
     y_value_t operator()(x_value_t) const;

    private:
     using Vector = Eigen::Matrix<y_value_t, Eigen::Dynamic, 1>;
     using Matrix = Eigen::SparseMatrix<y_value_t>;
     // Iterators to the function data.
     xIter xStart;
     xIter xFinish;
     yIter yStart;

     // Cubic spline coefficients
     Vector ypp;
};

// Definition of the main constructor.
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, CubicSplineBC left,
                                       y_value_t ypl, CubicSplineBC right,
                                       y_value_t ypr)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
     // Dimension of the linear system.
     const auto n = std::distance(xStart, xFinish);
     assert(n > 1);
     assert(std::is_sorted(xStart, xFinish));

     // Set up the sparse matrix.
     Matrix A(n, n);
     A.reserve(Eigen::VectorXi::Constant(n, 3));

     // set some constants
     constexpr auto oneThird =
         static_cast<x_value_t>(1) / static_cast<x_value_t>(3);
     constexpr auto oneSixth =
         static_cast<x_value_t>(1) / static_cast<x_value_t>(6);

     // Add in the upper diagonal.
     if (left == CubicSplineBC::Clamped) {
          A.insert(0, 1) = oneSixth * (xStart[1] - xStart[0]);
     }
     for (int i = 1; i < n - 1; i++) {
          A.insert(i, i + 1) = oneSixth * (xStart[i] - xStart[i - 1]);
     }

     // Add in the lower diagonal.
     for (int i = 0; i < n - 2; i++) {
          A.insert(i + 1, i) = oneSixth * (xStart[i + 1] - xStart[i]);
     }
     if (right == CubicSplineBC::Clamped) {
          A.insert(n - 1, n - 2) = oneSixth * (xStart[n - 1] - xStart[n - 2]);
     }

     // Add in the diagonal.
     if (left == CubicSplineBC::Free) {
          A.insert(0, 0) = static_cast<y_value_t>(1.0);
     } else {
          A.insert(0, 0) = oneThird * (xStart[1] - xStart[0]);
     }
     for (int i = 1; i < n - 1; i++) {
          A.insert(i, i) = oneThird * (xStart[i + 1] - xStart[i - 1]);
     }
     if (right == CubicSplineBC::Free) {
          A.insert(n - 1, n - 1) = static_cast<y_value_t>(1.0);
     } else {
          A.insert(n - 1, n - 1) = oneThird * (xStart[n - 2] - xStart[n - 1]);
     }

     // Finalise the matrix construction.
     A.makeCompressed();

     // Set the right hand side.
     ypp = Vector(n);
     if (left == CubicSplineBC::Free) {
          ypp(0) = static_cast<y_value_t>(0);
     } else {
          ypp(0) = (yStart[1] - yStart[0]) / (xStart[1] - xStart[0]) - ypl;
     }
     for (int i = 1; i < n - 1; i++) {
          ypp(i) = (yStart[i + 1] - yStart[i]) / (xStart[i + 1] - xStart[i]) -
                   (yStart[i] - yStart[i - 1]) / (xStart[i] - xStart[i - 1]);
     }
     if (right == CubicSplineBC::Free) {
          ypp(n - 1) = static_cast<y_value_t>(0);
     } else {
          ypp(n - 1) = (yStart[n - 1] - yStart[n - 2]) /
                           (xStart[n - 1] - xStart[n - 2]) -
                       ypr;
     }

     // Solve the linear system.
     Eigen::SparseLU<Matrix> solver;
     solver.compute(A);
     ypp = solver.solve(ypp);
     assert(solver.info() == Eigen::Success);
}

// Definition of the constructor for natural splines.
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart)
    : CubicSpline(xStart, xFinish, yStart, CubicSplineBC::Free, 0,
                  CubicSplineBC::Free, 0) {}

// Definition of the constructor when boundary conditions are the same.
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, CubicSplineBC both,
                                       y_value_t ypl, y_value_t ypr)
    : CubicSpline(xStart, xFinish, yStart, both, ypl, both, ypr) {}

// Evaluation of the interpolating function.
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::y_value_t CubicSpline<xIter, yIter>::operator()(
    x_value_t x) const {
  // Find the first element larger than x.
  auto iter = std::upper_bound(xStart, xFinish, x);
  // Adjust the iterator if out of range.
  if (iter == xStart) ++iter;
  if (iter == xFinish) --iter;
  // Perform the interpolation.
  constexpr auto oneSixth =
      static_cast<x_value_t>(1) / static_cast<x_value_t>(6);
  auto i2 = std::distance(xStart, iter);
  auto i1 = i2 - 1;
  auto x1 = xStart[i1];
  auto x2 = xStart[i2];
  auto h = x2 - x1;
  auto a = (x2 - x) / h;
  auto b = (x - x1) / h;
  return a * yStart[i1] + b * yStart[i2] +
         ((a * a * a - a) * ypp(i1) + (b * b * b - b) * ypp(i2)) *
             (h * h * oneSixth);
};

}   // namespace Interp

#endif   //  INTERP_CUBIC_SPLINE_GUARD_H
