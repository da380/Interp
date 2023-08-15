
#include <Interp/All>
#include <algorithm>
#include <cmath>
#include <fstream>
// #include <ios>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

int
main() {
     using namespace Interp;

     int n = 11, n2 = 101;
     std::vector<double> x(n), y(n), z(n), w(n), xx(n2);

     double xtest;

     // x is (0,...,1)
     // y = x
     // z = x^2
     for (int i = 0; i < n2; i++) {
          xx[i] = static_cast<double>(i) / static_cast<double>(n2 - 1);
     }
     for (int i = 0; i < n; i++) {
          x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
          y[i] = x[i];
          z[i] = x[i] * x[i];
          w[i] = x[i] * z[i] - z[i] + y[i];
          // std::cout << x[i] << std::endl;
     }

     // constructing CubicSpline object with x, y, use default natural BCs
     CubicSpline f(x.begin(), x.end(), y.begin());
     xtest = 0.1;
     // std::cout << "The true value of x is: " << xtest << std::endl;
     // std::cout << f(xtest) << std::endl;

     // constructing CubicSpline object with x, z, demonstrating all
     // combinations
     CubicSpline g(x.begin(), x.end(), z.begin(), "natural", 2, 2);
     CubicSpline h(x.begin(), x.end(), z.begin(), "clamped", 0, 2);
     CubicSpline m(x.begin(), x.end(), z.begin(), "natural", 2, "clamped", 2);
     // std::cout << "The true value of x^2 is: " << xtest * xtest << std::endl;
     // std::cout << "Natural and natural: " << g(xtest) << std::endl;
     // std::cout << "Clamped and clamped: " << h(xtest) << std::endl;
     // std::cout << "Natural and clamped: " << m(xtest) << std::endl;

     // constructing CubicSpline object with x, w, demonstrating all
     // combinations
     CubicSpline o(x.begin(), x.end(), w.begin(), "natural", 0, 0);
     CubicSpline p(x.begin(), x.end(), w.begin(), "clamped", 1, 2);
     CubicSpline q(x.begin(), x.end(), w.begin(), "natural", -2, "clamped", 2);
     // std::cout << "The true value of 1.1*x^3 - 2*x^2 + x is: "
     //           << xtest * (xtest * (1.1 * xtest - 2) + 1) << std::endl;
     // std::cout << "Natural and natural: " << o(xtest) << std::endl;
     // std::cout << "Clamped and clamped: " << p(xtest) << std::endl;
     // std::cout << "Natural and clamped: " << q(xtest) << std::endl;

     // std::cout << a << std::endl;

     // Write output
     const char *path =
         "/home/alex/Documents/c++/Interp/output_files/cubic_spline.out";
     std::ofstream MyFile(path, std::ios::trunc);
     // MyFile.open();

     // check
     for (int i = 0; i < n2; i++) {
          // MyFile << std::fixed;
          MyFile << xx[i] << ";";
          // MyFile << w[i] << ":";
          MyFile << p(xx[i]) << std::endl;

          // MyFile << xx[i];
          // MyFile << f(xx[i]) << std::endl;
          // std::cout << f(xx[i]) << std::endl;
     }
     MyFile.close();

     const char *path2 =
         "/home/alex/Documents/c++/Interp/output_files/cubic_spline.in";
     std::ofstream MyFile2(path2, std::ios::trunc);
     // MyFile.open();

     for (int i = 0; i < n; i++) {
          // MyFile << std::fixed;
          MyFile2 << x[i] << ";";
          MyFile2 << w[i] << std::endl;

          // MyFile << xx[i];
          // MyFile << f(xx[i]) << std::endl;
          // std::cout << f(xx[i]) << std::endl;
     }
     MyFile2.close();
}
