
#include "Interpolation/CubicSpline.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <numeric>
#include <random>
#include <vector>

int
main() {

    using namespace Interpolation;

    int n = 10;
    std::vector<double> x(n);
    std::vector<double> y(n);

    auto func = [](double x) { return std::sin(x); };

    for (int i = 0; i < n; i++) {
        x[i] = static_cast<double>(i) / static_cast<double>(n - 1);
        y[i] = func(x[i]);
    }

    // auto f = CubicSpline(x.begin(), x.end(), y.begin());
    auto f = Ranges::CubicSpline(x, y);

    int m = 50;
    std::ofstream file("Cubic.out");
    for (int i = 0; i < m; i++) {
        auto x = static_cast<double>(i) / static_cast<double>(m - 1);
        file << x << " " << func(x) << " " << f(x) << std::endl;
    }
}
