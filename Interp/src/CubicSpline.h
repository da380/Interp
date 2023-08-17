#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <algorithm>
#include <Eigen/Core>
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
    // y_value_t deriv(x_value_t) const;

   private:
    using Vector = Eigen::Matrix<y_value_t, Eigen::Dynamic, 1>;
    using Matrix = Eigen::SparseMatrix<y_value_t>;

    xIter xS;   // Iterator to the start of the x-values.
    xIter xF;   // Iterator to the end of hte x-values.
    yIter yS;   // Iterator to the start of the y-values.

    Vector ypp;   // Cubic spline coefficients
};

// Definition of the main constructor.
template <typename xIter, typename yIter>
  requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS,
                                       CubicSplineBC left, y_value_t ypl,
                                       CubicSplineBC right, y_value_t ypr)
    : xS{xS}, xF{xF}, yS{yS} {
    // Dimension of the linear system.
    const auto n = std::distance(xS, xF);
    assert(n > 1);
    assert(std::is_sorted(xS, xF));

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
        A.insert(0, 1) = oneSixth * (xS[1] - xS[0]);
    }
    for (int i = 1; i < n - 1; ++i) {
        A.insert(i, i + 1) = oneSixth * (xS[i] - xS[i - 1]);
    }

    // Add in the lower diagonal.
    for (int i = 0; i < n - 2; ++i) {
        A.insert(i + 1, i) = oneSixth * (xS[i + 1] - xS[i]);
    }
    if (right == CubicSplineBC::Clamped) {
        A.insert(n - 1, n - 2) = oneSixth * (xS[n - 1] - xS[n - 2]);
    }

    // Add in the diagonal.
    if (left == CubicSplineBC::Free) {
        A.insert(0, 0) = static_cast<y_value_t>(1.0);
    } else {
        A.insert(0, 0) = oneThird * (xS[1] - xS[0]);
    }
    for (int i = 1; i < n - 1; ++i) {
        A.insert(i, i) = oneThird * (xS[i + 1] - xS[i - 1]);
    }
    if (right == CubicSplineBC::Free) {
        A.insert(n - 1, n - 1) = static_cast<y_value_t>(1.0);
    } else {
        A.insert(n - 1, n - 1) = oneThird * (xS[n - 2] - xS[n - 1]);
    }

    // Finalise the matrix construction.
    A.makeCompressed();

    // Set the right hand side.
    ypp = Vector(n);
    if (left == CubicSplineBC::Free) {
        ypp(0) = static_cast<y_value_t>(0);
    } else {
        ypp(0) = (yS[1] - yS[0]) / (xS[1] - xS[0]) - ypl;
    }
    for (int i = 1; i < n - 1; i++) {
        ypp(i) = (yS[i + 1] - yS[i]) / (xS[i + 1] - xS[i]) -
                 (yS[i] - yS[i - 1]) / (xS[i] - xS[i - 1]);
    }
    if (right == CubicSplineBC::Free) {
        ypp(n - 1) = static_cast<y_value_t>(0);
    } else {
        ypp(n - 1) = (yS[n - 1] - yS[n - 2]) / (xS[n - 1] - xS[n - 2]) - ypr;
    }

    // Solve the linear system.
    Eigen::BiCGSTAB<Matrix> solver;
    solver.compute(A);
    ypp = solver.solve(ypp);
    assert(solver.info() == Eigen::Success);
}

// Definition of the constructor for natural splines.
template <typename xIter, typename yIter>
  requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS)
    : CubicSpline(xS, xF, yS, CubicSplineBC::Free, 0, CubicSplineBC::Free, 0) {}

// Definition of the constructor when boundary conditions are the same.
template <typename xIter, typename yIter>
  requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xS, xIter xF, yIter yS,
                                       CubicSplineBC both, y_value_t ypl,
                                       y_value_t ypr)
    : CubicSpline(xS, xF, yS, both, ypl, both, ypr) {}

// Evaluation of the interpolating function.
template <typename xIter, typename yIter>
  requires InterpolationIteratorPair<xIter, yIter>
CubicSpline<xIter, yIter>::y_value_t
CubicSpline<xIter, yIter>::operator()(x_value_t x) const {
    // Find the first element larger than x.
    auto iter = std::upper_bound(xS, xF, x);
    // Adjust the iterator if out of range.
    if (iter == xS)
        ++iter;
    if (iter == xF)
        --iter;
    // Perform the interpolation.
    constexpr auto oneSixth =
        static_cast<x_value_t>(1) / static_cast<x_value_t>(6);
    auto i2 = std::distance(xS, iter);
    auto i1 = i2 - 1;
    auto x1 = xS[i1];
    auto x2 = xS[i2];
    auto h = x2 - x1;
    auto a = (x2 - x) / h;
    auto b = (x - x1) / h;
    return a * yS[i1] + b * yS[i2] +
           ((a * a * a - a) * ypp(i1) + (b * b * b - b) * ypp(i2)) * h * h *
               oneSixth;
};

// // Evaluation of the derivative.
// template <typename xIter, typename yIter>
// CubicSpline<xIter, yIter>::y_value_t
// CubicSpline<xIter, yIter>::deriv(x_value_t x) const {
//     // Find the first element larger than x.
//     auto iter = std::upper_bound(xS, xF, x);
//     // Adjust the iterator if out of range.
//     if (iter == xS)
//         ++iter;
//     if (iter == xF)
//         --iter;
//     // Perform the interpolation.
//     constexpr auto oneSixth =
//         static_cast<x_value_t>(1) / static_cast<x_value_t>(6);
//     auto i2 = std::distance(xS, iter);
//     auto i1 = i2 - 1;
//     auto x1 = xS[i1];
//     auto x2 = xS[i2];
//     auto h = x2 - x1;
//     auto a = (x2 - x) / h;
//     auto b = (x - x1) / h;
//     return a * yS[i1] + b * yS[i2] +
//            ((a * a * a - a) * ypp(i1) + (b * b * b - b) * ypp(i2)) * h * h *
//                oneSixth;
// };

}   // namespace Interp

#endif   //  INTERP_CUBIC_SPLINE_GUARD_H
