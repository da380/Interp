#pragma once

#include "../Concepts.h"
#include "Function.h"
#include "Utility.h"

#include <algorithm>
#include <cassert>
#include <variant>

#include <Eigen/Core>
#include <Eigen/SparseCholesky>

namespace Interpolation {

namespace OneDimensional {

// Define classes for boundary conditions.
template <RealOrComplexFloatingPoint Scalar = double> class CubicSplineBC {
  public:
    constexpr CubicSplineBC() = default;

    constexpr CubicSplineBC(Scalar value) : _value{value} {}

    constexpr auto Clamped() const {
        return std::holds_alternative<Scalar>(_value);
    }

    constexpr auto operator()() const { return std::get<Scalar>(_value); }

  private:
    std::variant<std::monostate, Scalar> _value;
};

// Forward declare class.
template <RealView XView, RealOrComplexView YView>
    requires DataViews1D<XView, YView>
class CubicSpline;

// Set traits.
namespace Internal {

template <RealView XView, RealOrComplexView YView>
struct Traits<CubicSpline<XView, YView>> {
    using Real = std::ranges::range_value_t<XView>;
    using Scalar = std::ranges::range_value_t<YView>;
};
}   // namespace Internal

// Define the class.
template <RealView XView, RealOrComplexView YView>
    requires DataViews1D<XView, YView>
class CubicSpline : public Function<CubicSpline<XView, YView>> {

  public:
    using Real = typename Function<CubicSpline<XView, YView>>::Real;
    using Scalar = typename Function<CubicSpline<XView, YView>>::Scalar;

    CubicSpline() = delete;

    CubicSpline(XView x, YView y,
                CubicSplineBC<Scalar> &&left = CubicSplineBC<Scalar>(),
                CubicSplineBC<Scalar> &&right = CubicSplineBC<Scalar>())
        : _x{x}, _y{y} {
        assert(_x.size() == _y.size());
        assert(std::ranges::is_sorted(_x));

        // Set up the sparse matrix.
        const auto n = _x.size();
        Matrix A(n, n);
        A.reserve(Eigen::VectorXi::Constant(n, 3));

        // set some constants
        constexpr auto oneThird = static_cast<Real>(1) / static_cast<Real>(3);
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);

        // Add in the upper diagonal.
        if (left.Clamped()) {
            A.insert(0, 1) = oneSixth * (_x[1] - _x[0]);
        }
        for (int i = 1; i < n - 1; ++i) {
            A.insert(i, i + 1) = oneSixth * (_x[i] - _x[i - 1]);
        }

        // Add in the lower diagonal.
        for (int i = 0; i < n - 2; ++i) {
            A.insert(i + 1, i) = oneSixth * (_x[i + 1] - _x[i]);
        }
        if (right.Clamped()) {
            A.insert(n - 1, n - 2) = oneSixth * (_x[n - 1] - _x[n - 2]);
        }

        // Add in the diagonal.
        if (!left.Clamped()) {
            A.insert(0, 0) = 1;
        } else {
            A.insert(0, 0) = oneThird * (_x[1] - _x[0]);
        }
        for (int i = 1; i < n - 1; ++i) {
            A.insert(i, i) = oneThird * (_x[i + 1] - _x[i - 1]);
        }
        if (!right.Clamped()) {
            A.insert(n - 1, n - 1) = 1;
        } else {
            A.insert(n - 1, n - 1) = oneThird * (_x[n - 1] - _x[n - 2]);
        }

        // Finalise the matrix construction.
        A.makeCompressed();

        // Set the right hand side.
        Vector rhs(n);
        if (!left.Clamped()) {
            rhs(0) = 0;
        } else {
            rhs(0) = (_y[1] - _y[0]) / (_x[1] - _x[0]) - left();
        }
        for (int i = 1; i < n - 1; i++) {
            rhs(i) = (_y[i + 1] - _y[i]) / (_x[i + 1] - _x[i]) -
                     (_y[i] - _y[i - 1]) / (_x[i] - _x[i - 1]);
        }
        if (!right.Clamped()) {
            rhs(n - 1) = 0;
        } else {
            rhs(n - 1) =
                right() - (_y[n - 1] - _y[n - 2]) / (_x[n - 1] - _x[n - 2]);
        }

        // Solve the linear system.
        Eigen::SimplicialLDLT<Matrix> solver;
        solver.compute(A);
        _ypp = solver.solve(rhs);
        assert(solver.info() == Eigen::Success);
    }

    template <int N>
    auto Evaluate(Real x) const
        requires(N <= 2)
    {
        constexpr auto oneSixth = static_cast<Real>(1) / static_cast<Real>(6);
        auto i2 = std::distance(_x.begin(), Find(_x, x));
        auto i1 = i2 - 1;
        auto x1 = _x[i1];
        auto x2 = _x[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        if constexpr (N == 0) {
            return a * _y[i1] + b * _y[i2] +
                   ((a * a * a - a) * _ypp(i1) + (b * b * b - b) * _ypp(i2)) *
                       h * h * oneSixth;
        } else if constexpr (N == 1) {
            return (_y[i2] - _y[i1]) / h + oneSixth * h *
                                               ((-3 * a * a + 1) * _ypp(i1) +
                                                (3 * b * b - 1) * _ypp(i2));
        } else {
            return a * _ypp(i1) + b * _ypp(i2);
        }
    }

  private:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix = Eigen::SparseMatrix<Scalar>;
    XView _x;
    YView _y;

    Vector _ypp;
};

// Deduction guides for construction from ranges.
template <std::ranges::viewable_range XView, std::ranges::viewable_range YView>
CubicSpline(XView &&,
            YView &&) -> CubicSpline<std::ranges::views::all_t<XView>,
                                     std::ranges::views::all_t<YView>>;

}   // namespace OneDimensional

}   // namespace Interpolation
