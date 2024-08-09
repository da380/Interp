#include <Interpolation/All>
#include <algorithm>
#include <iostream>
#include <ranges>
#include <span>
#include <vector>

int
main() {

    using Real = double;

    auto n = 10;

    auto x = std::vector<Real>();
    auto y = std::vector<Real>();

    std::ranges::copy(std::ranges::views::iota(0, n), std::back_inserter(x));

    std::ranges::transform(x, std::back_inserter(y), [](auto x) { return x; });

    auto f = Interpolation::OneDimensional::Linear(x, y);

    auto g = Interpolation::OneDimensional::CubicSpline(x, y);

    auto h = Interpolation::OneDimensional::Derivative(g);

    for (auto p : x)
        std::cout << g(p) - p << " " << h(p) - 1 << std::endl;
}