
#include <Interp/All>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

int main() {
  using namespace Interp;

  int n = 10;
  std::vector<double> x(n), y(n);

  auto func = [](double x) { return std::sin(x); };

  for (int i = 0; i < n; i++) {
    x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
    y[i] = func(x[i]);
  }

  CubicSpline f(x.begin(), x.end(), y.begin());
  std::cout << a << std::endl;


