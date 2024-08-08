#ifndef INTERPOLATION_AKIMA_SPLINE_GUARD_H
#define INTERPOLATION_AKIMA_SPLINE_GUARD_H

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include "Concepts.h"

namespace Interpolation {
template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
class Akima {
  public:
    // Define some class member types
    using x_value_type = std::iter_value_t<xIter>;
    using y_value_type = std::iter_value_t<yIter>;

    // Declare the constructor.
    Akima(xIter, xIter, yIter);

    // Declare the application operator.
    y_value_type operator()(x_value_type) const;
    y_value_type deriv(x_value_type) const;

  private:
    // Iterators to the function data.
    xIter _xS;
    xIter _xF;
    yIter _yS;

    // m and s values
    std::vector<y_value_type> m;
    std::vector<y_value_type> s;
};

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::Akima(xIter xS, xIter xF, yIter yS)
    : _xS{xS}, _xF{xF}, _yS{yS} {
    // Dimension of the linear system.
    const auto n = std::distance(_xS, _xF);
    assert(n > 1);
    assert(std::is_sorted(_xS, _xF));

    // fill out m
    for (int i = 0; i < n - 1; ++i) {
        m.push_back((_yS[i + 1] - yS[i]) / (_xS[i + 1] - _xS[i]));
    }

    // fill out s
    double onehalf = static_cast<double>(1) / static_cast<double>(2);
    s.push_back(m[0]);
    s.push_back((m[0] + m[1]) * onehalf);
    for (int i = 2; i < n - 2; ++i) {
        // y_value_type a = std::sqrt(std::norm(m[i + 1] - m[i]));
        y_value_type a = std::abs(m[i + 1] - m[i]);

        //     y_value_type b = std::sqrt(std::norm(m[i - 1] - m[i -
        //     2]));
        y_value_type b = std::abs(m[i - 1] - m[i - 2]);
        // std::cout << m[i - 1] - m[i - 2] << "   " << b << std::endl;
        // if (a == 0 && b != 0) {
        //    s.push_back(m[i]);
        // } else if (a != 0 && b == 0) {
        //     s.push_back(m[i - 1]);
        //	} else if (a == 0 && b == 0) {
        //     s.push_back((m[i - 1] + m[i]) * onehalf);
        // } else {
        //	  s.push_back((a * m[i - 1] + b * m[i]) / (a +
        // b));
        //    }
        if ((a + b) == 0) {
            s.push_back((m[i - 1] + m[i]) / 2);
        } else {
            s.push_back((a * m[i - 1] + b * m[i]) / (a + b));
        }
    }
    s.push_back((m[n - 3] + m[n - 2]) * onehalf);
    s.push_back(m[n - 2]);
    // std::cout << "Distance: " << std::distance(xS, xF) << std::endl;
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::y_value_type
Akima<xIter, yIter>::operator()(x_value_type x) const {
    // Find the first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    const auto n = std::distance(_xS, iter);
    // std::cout << "n: " << n << std::endl;
    // std::cout << "iter[-1]: " << iter[-1] << std::endl;
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    auto a = _yS[n - 1];
    auto b = s[n - 1];
    auto c = (static_cast<y_value_type>(3.0) * m[n - 1] -
              static_cast<y_value_type>(2.0) * s[n - 1] - s[n]) /
             (iter[0] - iter[-1]);
    auto d = (s[n - 1] + s[n] - static_cast<y_value_type>(2.0) * m[n - 1]) /
             ((iter[0] - iter[-1]) * (iter[0] - iter[-1]));
    return a + (x - iter[-1]) * (b + (x - iter[-1]) * (c + d * (x - iter[-1])));
}

template <typename xIter, typename yIter>
    requires InterpolationIteratorPair<xIter, yIter>
Akima<xIter, yIter>::y_value_type
Akima<xIter, yIter>::deriv(x_value_type x) const {
    // Find the first element larger than x.
    auto iter = std::upper_bound(_xS, _xF, x);
    const auto n = std::distance(_xS, iter);
    // std::cout << "n: " << n << std::endl;
    // std::cout << "iter[-1]: " << iter[-1] << std::endl;
    // Adjust the iterator if out of range.
    if (iter == _xS)
        ++iter;
    if (iter == _xF)
        --iter;
    // Perform the interpolation.
    // auto a = yS[n - 1];
    auto b = s[n - 1];
    auto c = (static_cast<y_value_type>(3.0) * m[n - 1] -
              static_cast<y_value_type>(2.0) * s[n - 1] - s[n]) /
             (iter[0] - iter[-1]);
    auto d = (s[n - 1] + s[n] - static_cast<y_value_type>(2.0) * m[n - 1]) /
             ((iter[0] - iter[-1]) * (iter[0] - iter[-1]));
    return (b + (x - iter[-1]) * (2.0 * c + 3.0 * d * (x - iter[-1])));
}

}   // namespace Interpolation

#endif   // INTERPOLATION_AKIMA_SPLINE_GUARD_H
