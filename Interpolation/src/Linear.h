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
  using x_value_t = std::iter_value_t<xIter>;
  using y_value_t = std::iter_value_t<yIter>;

  // Declare the constructor.
  Linear(xIter, xIter, yIter);

  // Declare the application operator.
  y_value_t operator()(x_value_t) const;

  // Declare the derivative operator.
  y_value_t Derivative(x_value_t) const;

 private:
  // Iterators to the function data.
  xIter xS;
  xIter xF;
  yIter yS;
};

template <typename xIter, typename yIter>
requires InterpolationIteratorPair<xIter, yIter> Linear<xIter, yIter>::Linear(
    xIter xS, xIter xF, yIter yS)
    : xS{xS}, xF{xF}, yS{yS} {}

template <typename xIter, typename yIter>
requires InterpolationIteratorPair<xIter, yIter> Linear<xIter, yIter>::y_value_t
Linear<xIter, yIter>::operator()(const x_value_t x) const {
  // Find first element larger than x.
  auto iter = std::upper_bound(xS, xF, x);
  // Adjust the iterator if out of range.
  if (iter == xS) ++iter;
  if (iter == xF) --iter;
  // Perform the interpolation.
  auto i2 = std::distance(xS, iter);
  auto i1 = i2 - 1;
  auto x1 = xS[i1];
  auto x2 = xS[i2];
  auto h = x2 - x1;
  auto a = (x2 - x) / h;
  auto b = (x - x1) / h;
  return a * yS[i1] + b * yS[i2];
}

template <typename xIter, typename yIter>
requires InterpolationIteratorPair<xIter, yIter> Linear<xIter, yIter>::y_value_t
Linear<xIter, yIter>::Derivative(const x_value_t x)
const {
  // Find first element larger than x.
  auto iter = std::upper_bound(xS, xF, x);
  // Adjust the iterator if out of range.
  if (iter == xS) ++iter;
  if (iter == xF) --iter;
  // Perform the interpolation.
  auto i2 = std::distance(xS, iter);
  auto i1 = i2 - 1;
  auto x1 = xS[i1];
  auto x2 = xS[i2];
  auto h = x2 - x1;
  auto a = (x2 - x) / h;
  auto b = (x - x1) / h;
  return (yS[i2] - yS[i1]) / h;
}

}  // namespace Interpolation

#endif  //  INTERPOLATION_LINEAR_GUARD_H
