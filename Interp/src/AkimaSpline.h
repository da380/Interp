#ifndef INTERP_AKIMA_SPLINE_GUARD_H
#define INTERP_AKIMA_SPLINE_GUARD_H

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interp {
template <typename xIter, typename yIter>
     requires InterpolationIteratorPair<xIter, yIter>
class Akima {
    public:
     // Define some class member types
     using x_value_t = std::iter_value_t<xIter>;
     using y_value_t = std::iter_value_t<yIter>;

     // Declare the constructor.
     Akima(xIter, xIter, yIter);

     // Declare the application operator.
     y_value_t operator()(x_value_t) const;

    private:
     // Iterators to the function data.
     xIter xS;
     xIter xF;
     yIter yS;

     // m and s values
     std::vector<y_value_t> slopes;
     std::vector<y_value_t> splslopes;
};

template <typename xIter, typename yIter>
     requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::Akima(xIter xS, xIter xF, yIter yS)
    : xS{xS}, xF{xF}, yS{yS} {
     // Dimension of the linear system.
     const auto n = std::distance(xS, xF);
     assert(n > 1);
     assert(std::is_sorted(xS, xF));

     // fill out slopes
     for (int i = 0; i < n - 1; ++i) {
          slopes.push_back((yS[i + 1] - yS[i]) / (xS[i + 1] - xS[i]));
     }

     // fill out splslopes
     double onehalf = static_cast<double>(2);
     splslopes.push_back(slopes[0]);
     splslopes.push_back((slopes[0] + slopes[1]) * onehalf);
     for (int i = 2; i < n - 2; ++i) {
          y_value_t a = std::sqrt(std::norm(slopes[i + 1] - slopes[i]));

          y_value_t b = std::sqrt(std::norm(slopes[i - 1] - slopes[i - 2]));
          // if (a == 0 && b != 0) {
          //    splslopes.push_back(slopes[i]);
          // } else if (a != 0 && b == 0) {
          //     splslopes.push_back(slopes[i - 1]);
          //	} else if (a == 0 && b == 0) {
          //     splslopes.push_back((slopes[i - 1] + slopes[i]) * onehalf);
          // } else {
          //	  splslopes.push_back((a * slopes[i - 1] + b * slopes[i]) / (a +
          // b));
          //    }
          splslopes.push_back((a * slopes[i - 1] + b * slopes[i]) / (a + b));
     }
     splslopes.push_back((slopes[n - 3] + slopes[n - 2]) * onehalf);
     splslopes.push_back(slopes[n - 2]);
}

template <typename xIter, typename yIter>
     requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::y_value_t
Akima<xIter, yIter>::operator()(x_value_t x) const {
     // Find the first element larger than x.
     auto iter = std::upper_bound(xS, xF, x);
     const auto m = std::distance(xS, iter);

     // Adjust the iterator if out of range.
     if (iter == xS)
          ++iter;
     if (iter == xF)
          --iter;
     // Perform the interpolation.
     auto a = iter[-1];
     auto b = splslopes[m - 2];
     auto c =
         (static_cast<y_value_t>(3.0) * slopes[m - 2] -
          static_cast<y_value_t>(2.0) * splslopes[m - 2] - splslopes[m - 1]) /
         (iter[0] - iter[-1]);
     auto d = (splslopes[m - 2] + splslopes[m - 1] -
               static_cast<y_value_t>(2.0) * slopes[m - 2]) /
              ((iter[0] - iter[-1]) * (iter[0] - iter[-1]));
     return a +
            (x - iter[-1]) * (b + (x - iter[-1]) * (c + d * (x - iter[-1])));
}

}   // namespace Interp

#endif
