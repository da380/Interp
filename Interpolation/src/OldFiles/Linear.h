#ifndef INTERPOLATION_LINEAR_GUARD_H
#define INTERPOLATION_LINEAR_GUARD_H

#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interpolation {

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class Linear {
  public:
    // Define some class member types
    using x_value_type = std::iter_value_t<xIter>;
    using y_value_type = std::iter_value_t<yIter>;

    // Declare the constructor.
    Linear(xIter, xIter, yIter);

    // Declare the application operator.
    y_value_type operator()(x_value_type) const;

    // Declare the derivative operator.
    y_value_type Derivative(x_value_type) const;

  private:
    xIter _xS;   // Iterator to start of x values.
    xIter _xF;   // Iterator to end of x values.
    yIter _yS;   // Iterator to start of y values.
};

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::Linear(xIter xS, xIter xF, yIter yS)
    : _xS{xS}, _xF{xF}, _yS{yS} {}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::y_value_type
Linear<xIter, yIter>::operator()(const x_value_type x) const {
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

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Linear<xIter, yIter>::y_value_type
Linear<xIter, yIter>::Derivative(const x_value_type x) const {
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

}   // namespace Interpolation

#endif   //  INTERPOLATION_LINEAR_GUARD_H
