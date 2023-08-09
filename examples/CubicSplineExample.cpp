
#include <Interp/All>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

int
main() {
     using namespace Interp;

     int n = 11;
     std::vector<double> x(n), y(n);
     double xtest;
     std::string bcarg;

     // auto func = [](double x) { return std::sin(x); };

     for (int i = 0; i < n; i++) {
          x[i] = static_cast<double>(i);
          // std::cout << x[i] << std::endl;
          y[i] = x[i] * x[i];
          // y[i] = func(x[i]);
     }
     //  std::cout << x[9] << std::endl;

     // std::cout << x[0] << std::endl;
     // std::cout << x[10] << std::endl;
     // for (int i = 0; i < n; i++) {
     //      y[i] = x[i];
     // }
     // std::cout << x[0] << std::endl;
     // std::cout << x[10] << std::endl;
     bcarg = "natural";
     CubicSpline f(x.begin(), x.end(), y.begin(), "clamped", 0, "clamped", 20);
     xtest = 0.65;
     // std::cout << f(xtest) << std::endl;
     // std::cout << xtest * xtest << std::endl;
     // std::cout << a << std::endl;
}
