
#include <Interp/All>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

int
main() {
     using namespace Interp;
     using Float = double;
     using x_value_t = Float;
     // using y_value_t = std::complex<Float>;
     using y_value_t = Float;

     int n = 10;
     std::vector<x_value_t> x(n);
     std::vector<y_value_t> y(n);

     auto func = [](x_value_t x) {
          if constexpr (ComplexFloatingPoint<y_value_t>) {
               return std::complex<x_value_t>(std::sin(x), std::cos(x));
          }
          if constexpr (RealFloatingPoint<y_value_t>) {
               return std::sin(x);
          }
     };

     for (int i = 0; i < n; i++) {
          x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
          y[i] = func(x[i]);
     }

     auto f = CubicSpline(x.begin(), x.end(), y.begin());

     int m = 50;
     std::ofstream file("CubicSpline.out");
     for (int i = 0; i < m; i++) {
          auto x = static_cast<double>(i) / static_cast<double>(m - 1);
          file << x << " " << func(x) << " " << f(x) << std::endl;
     }
}
