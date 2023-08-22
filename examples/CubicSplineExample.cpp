
#include <Interp/All>
#include <algorithm>
#include <chrono>
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

     int n = 11;
     std::vector<x_value_t> x(n);
     std::vector<y_value_t> y(n);
     const x_value_t pi = 3.1415926535;

     auto func = [](x_value_t x) {
          if constexpr (ComplexFloatingPoint<y_value_t>) {
               return std::complex<x_value_t>(std::sin(x), std::cos(x));
          }
          if constexpr (RealFloatingPoint<y_value_t>) {
               return std::sin(x);
          }
     };

     for (int i = 0; i < n; i++) {
          x[i] = pi * static_cast<double>(i) / static_cast<double>(n - 1);
          y[i] = func(x[i]);
     }

     using namespace std::chrono;
     auto start = high_resolution_clock::now();
     auto f = CubicSpline(x.begin(), x.end(), y.begin());
     auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);
     std::cout << duration.count() / static_cast<double>(1000000) << std::endl;
     start = high_resolution_clock::now();
     auto g = Akima(x.begin(), x.end(), y.begin());
     stop = high_resolution_clock::now();
     duration = duration_cast<microseconds>(stop - start);
     std::cout << duration.count() / static_cast<double>(1000000) << std::endl;
     //     const char *path =
     //     "/home/adcm2/raidam/Interp/output_files/cubic_spline.in";
     // std::ofstream file(path);
     // for (int i = 0; i < n; i++) {
     //      // auto x = static_cast<double>(i) / static_cast<double>(m -
     //      1); file << x[i] << ";" << y[i] << std::endl;
     // }
     // file.close();

     const char *path2 = "/home/adcm2/raidam/Interp/output_files/C+Akima.out";
     std::ofstream file2(path2);
     int m = 51;
     for (int i = 0; i < m; i++) {
          auto xx = pi * static_cast<double>(i) / static_cast<double>(m - 1);
          file2 << xx << ";" << func(xx) << ";" << f(xx) << ";" << g(xx)
                << std::endl;
          // std::cout << f(xx) << "     " << g(xx) << std::endl;
     }
     file2.close();

     // g(0.09);
     // g(0.1);
}
