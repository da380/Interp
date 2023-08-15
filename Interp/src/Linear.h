#ifndef INTERP_LINEAR_GUARD_H
#define INTERP_LINEAR_GUARD_H

#include <algorithm>
#include <concepts>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interp {

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
Linear<xIter, yIter>::y_value_t
Linear<xIter, yIter>::operator()(const x_value_t x) const {
     // Find first element larger than x.
     auto iter = std::upper_bound(xStart, xFinish, x);
     // Adjust the iterator if out of range.
     if (iter == xStart)
          ++iter;
     if (iter == xFinish)
          --iter;
     // Perform the interpolation.
     auto i2 = std::distance(xStart, iter);
     auto i1 = i2 - 1;
     auto x1 = xStart[i1];
     auto x2 = xStart[i2];
     auto y1 = yStart[i1];
     auto y2 = yStart[i2];
     return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

}   // namespace Interp

#endif  //  INTERP_LINEAR_GUARD_H
