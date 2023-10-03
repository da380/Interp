#ifndef INTERPOLATION_LAGRANGE_GUARD_H
#define INTERPOLATION_LAGRANGE_GUARD_H

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <concepts>
#include <iterator>
#include <numeric>
#include <ranges>
#include <vector>

#include "Concepts.h"

namespace Interpolation {

template <RealFloatingPointIterator I>
class LagrangePolynomial {
 public:
  using value_t = std::iter_value_t<I>;

  // Constructors.
  LagrangePolynomial() = delete;
  LagrangePolynomial(I start, I finish)
      : n{std::distance(start, finish)}, X{start} {
    assert(n > 0);
    assert(std::is_sorted(start, finish));
  }

  // Evaluation of ith function at x.
  value_t operator()(int i, value_t x) const {
    auto prod1 = static_cast<value_t>(1);
    auto prod2 = static_cast<value_t>(1);
    for (int j = 0; j < n; j++) {
      if (i != j) {
        prod1 *= x - X[j];
        prod2 *= X[i] - X[j];
      }
    }
    return prod1 / prod2;
  }

  // Derivative function.
  value_t Derivative(int i, value_t x) const {
    auto hp = static_cast<value_t>(0);
    auto prod2 = static_cast<value_t>(1);
    for (int j = 0; j < n; j++) {
      if (i != j) {
        auto prod1 = static_cast<value_t>(1);
        for (int k = 0; k < n; k++) {
          if (k != i && k != j) prod1 *= x - X[k];
        }
        hp += prod1;
        prod2 *= X[i] - X[j];
      }
    }
    return hp / prod2;
  }

 private:
  std::ptrdiff_t n;  // Number of nodes.
  I X;               // Iterator to the start of the nodes.
};

template <typename xIter, typename yIter>
requires InterpolationIteratorPair<xIter, yIter>
class Lagrange {
 public:
  // Define some class member types
  using x_value_t = std::iter_value_t<xIter>;
  using y_value_t = std::iter_value_t<yIter>;

  // Constructors.
  Lagrange(xIter xS, xIter xF, yIter yS)
      : xS{xS}, xF{xF}, yS{yS}, h{LagrangePolynomial(xS, xF)} {}

  // Return number of points
  auto size() const { return std::distance(xS, xF); }

  // Evaluation functions
  y_value_t operator()(x_value_t x) const {
    auto y = static_cast<y_value_t>(0);
    for (int i = 0; i < size(); i++) {
      y += h(i, x) * yS[i];
    }
    return y;
  }
  y_value_t Derivative(x_value_t x) const {
    auto yp = static_cast<y_value_t>(0);
    for (int i = 0; i < size(); i++) {
      yp += h.Derivative(i, x) * yS[i];
    }
    return yp;
  }

 private:
  xIter xS;  // Iterator to the start of the x-values.
  xIter xF;  // Iterator to the end of hte x-values.
  yIter yS;  // Iterator to the start of the y-values.

  LagrangePolynomial<xIter> h;  // Lagrange Polynomial for interpolation.
};

}  // namespace Interpolation

#endif
