// #ifndef INTERP_AKIMA_SPLINE_GUARD_H
// #define INTERP_AKIMA_SPLINE_GUARD_H

// #include <Eigen/IterativeLinearSolvers>
// #include <Eigen/SparseCore>
// #include <algorithm>
// #include <eigen3/Eigen/Core>
// #include <iostream>
// #include <iterator>
// #include <vector>

// #include "Concepts.h"

// namespace Interp {
// template <typename xIter, typename yIter>
//     requires InterpolationIteratorPair<xIter, yIter>
// class Akima {
//    public:
//     // Define some class member types
//     using x_value_t = std::iter_value_t<xIter>;
//     using y_value_t = std::iter_value_t<yIter>;

//     // Declare the constructor.
//     Akima(xIter, xIter, yIter);

//     // Declare the application operator.
//     y_value_t operator()(x_value_t) const;

//    private:
//     // Iterators to the function data.
//     xIter xS;
//     xIter xF;
//     yIter yS;

//     // m and s values
//     std::vector<y_value_t> slopes;
//     std::vector<y_value_t> splslopes;
// };

// template <typename xIter, typename yIter>
// Akima<xIter, yIter>::Akima(xIter xS, xIter xF, yIter yS)
//     : xS{xS}, xF{xF}, yS{yS} {
//     // Dimension of the linear system.
//     const auto n = std::distance(xS, xF);
//     assert(n > 1);
//     assert(std::is_sorted(xS, xF));

//     // fill out slopes
//     for (int i = 0; i < n - 1; ++i) {
//         slopes.push_back((yS[i + 1] - yS[i]) / (x[i + 1] - x[i]));
//     }

//     // fill out splslopes
//     double onehalf = static_cast<double>(2);
//     splslopes.push_back(slopes[0]);
//     splslopes.push_back((slopes[0] + slopes[1]) * onehalf);
//     for (int i = 2; i < n - 2; ++i) {
//         y_value_t a = slopes[i + 1] - slopes[i];

//         y_value_t b = slopes[i - 1] - slopes[i - 2];
//         if (a < 0) {
//             a = -a;
//         }
//         if (b < 0) {
//             b = -b;
//         }
//         if (a == 0 && b != 0) {
//             splslopes.push_back(slopes[i]);
//         } else if (a != 0 && b == 0) {
//             splslopes.push_back(slopes[i - 1]);
//         } else if (a == 0 && b == 0) {
//             splslopes.push_back((slopes[i - 1] + slopes[i]) * onehalf);
//         } else {
//             splslopes.push_back((a * slopes[i - 1] + b * slopes[i]) / (a +
//             b))
//         }
//     }
//     splslopes.push_back((slopes[n - 3] + slopes[n - 2]) * onehalf);
//     splslopes.push_back(slopes[n - 2]);
// }

// }   // namespace Interp

// #endif