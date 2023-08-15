#ifndef INTERP_MAKE_DATA_1D_GUARD_H
#define INTERP_MAKE_DATA_1D_GUARD_H

template <typename xVec, typename yVec, typename Function>
void MakeData1D(xVec x, yVec y, Function func) {
  using x_value_t = xVec::value_type;

  // generate a random size for the data
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> d(10, 1000);
  int n = d(gen);

  // Set the x values
  const x_value_t x1 = 0.0;
  const x_value_t x2 = std::numbers::pi_v<x_value_t>;
  std::uniform_real_distribution<x_value_t> dr{x1, x2};
  std::generate_n(std::back_inserter(x), n, [&]() { return dr(gen); });
  std::sort(x.begin(), x.end());

  // Set the y values
  std::transform(x.begin(), x.end(), std::back_inserter(y), func);
}

#endif  // INTERP_MAKE_DATA_1D_GUARD_H
