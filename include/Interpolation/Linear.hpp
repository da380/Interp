#ifndef INTERPOLATION_LINEAR_GUARD_H
#define INTERPOLATION_LINEAR_GUARD_H

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <vector>

#include "NumericConcepts/Iterators.hpp"

namespace Interpolation {

template <NumericConcepts::RealIterator xIter,
          NumericConcepts::RealOrComplexIterator yIter>
    requires NumericConcepts::SameIteratorPrecision<xIter, yIter>
class Linear {
  public:
    // Define some class member types
    using x_value_type = std::iter_value_t<xIter>;
    using y_value_type = std::iter_value_t<yIter>;

    // Default constructor.
    Linear() = default;

    // Declare the constructor.
    Linear(xIter xS, xIter xF, yIter yS) : _xS{xS}, _xF{xF}, _yS{yS} {}

    // Declare the application operator.
    y_value_type operator()(const x_value_type x) const {
        // Find first element larger than x.
        auto iter = std::upper_bound(_xS, _xF, x);
        // Adjust the iterator if out of range.
        if (iter == _xS)
            ++iter;
        if (iter == _xF)
            --iter;
        // Perform the interpolation.
        auto i2 = std::distance(_xS, iter);
        auto i1 = i2 - 1;
        auto x1 = _xS[i1];
        auto x2 = _xS[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        return a * _yS[i1] + b * _yS[i2];
    }

    // Declare the derivative operator.
    y_value_type Derivative(const x_value_type x) const {
        // Find first element larger than x.
        auto iter = std::upper_bound(_xS, _xF, x);
        // Adjust the iterator if out of range.
        if (iter == _xS)
            ++iter;
        if (iter == _xF)
            --iter;
        // Perform the interpolation.
        auto i2 = std::distance(_xS, iter);
        auto i1 = i2 - 1;
        auto x1 = _xS[i1];
        auto x2 = _xS[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        return (_yS[i2] - _yS[i1]) / h;
    }

  private:
    xIter _xS;   // Iterator to start of x values.
    xIter _xF;   // Iterator to end of x values.
    yIter _yS;   // Iterator to start of y values.
};

namespace Ranges {

namespace NC = NumericConcepts;
template <NC::RealView X, NC::RealOrComplexView Y>
    requires NC::SameRangePrecision<X, Y>
class Linear {
  public:
    // Define some class member types
    using Real = std::ranges::range_value_t<X>;
    using Scalar = std::ranges::range_value_t<Y>;

    // Default constructor.
    Linear() = default;

    // General constructor.
    Linear(X x, Y y) : _x{x}, _y{y} {}

    // Evaluate interpolating function.
    auto operator()(Real x) const {
        // Find the first element larger than x.
        auto it = std::ranges::upper_bound(_x, x);
        // Adjust the iterator if out of range.
        if (it == _x.begin())
            ++it;
        if (it == _x.end())
            --it;
        // Perform the interpolation.
        auto i2 = std::distance(_x.begin(), it);
        auto i1 = i2 - 1;
        auto x1 = _x[i1];
        auto x2 = _x[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        return a * _y[i1] + b * _y[i2];
    }

    auto Derivative(Real x) const {
        // Find the first element larger than x.
        auto it = std::ranges::upper_bound(_x, x);
        // Adjust the iterator if out of range.
        if (it == _x.begin())
            ++it;
        if (it == _x.end())
            --it;
        // Perform the interpolation.
        auto i2 = std::distance(_x.begin(), it);
        auto i1 = i2 - 1;
        auto x1 = _x[i1];
        auto x2 = _x[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        return (_y[i2] - _y[i1]) / h;
    }

  private:
    X _x;   // View to the ordinates
    Y _y;   // View to the values.
};

// Deductions guides for construction from ranges.
template <NC::RealRange X, NC::RealOrComplexRange Y>
Linear(X &&, Y &&)
    -> Linear<std::ranges::views::all_t<X>, std::ranges::views::all_t<Y>>;
}   // namespace Ranges

}   // namespace Interpolation

#endif   //  INTERPOLATION_LINEAR_GUARD_H
