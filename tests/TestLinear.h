#ifndef INTERP_TEST_LINEAR_GUARD_H
#define INTERP_TEST_LINEAR_GUARD_H

#include <Interp/Linear>
#include <complex>
#include <concepts>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "MakeData1D.h"

template <std::floating_point Float>
int LinearCheckNodalValues() {
  using namespace Interp;

  auto func = [](Float x) {
    return x * std::sin(2 * x);
  };

  std::vector<Float> x, y;
  MakeData1D(x, y, func);

  Linear f(x.begin(), x.end(), y.begin());

  constexpr auto eps = 100 * std::numeric_limits<Float>::epsilon();
  for (int i = 0; i < x.size(); i++) {
    auto xx = x[i];
    auto f1 = func(xx);
    auto f2 = f(xx);
    auto diff = std::abs(f2 - f2) / std::abs(f1);
    if (diff > eps) return 1;
  }

  return 0;
}

#endif  // INTERP_TEST_LINEAR_GUARD_H
