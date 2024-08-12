#include <Interpolation/All>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <span>
#include <vector>

int
main() {

    using Real = double;

    auto a = std::vector<double>{1, 1};

    std::cout << Interpolation::OneDimensional::PolynomialEvaluate(a, 1)
              << std::endl;
}