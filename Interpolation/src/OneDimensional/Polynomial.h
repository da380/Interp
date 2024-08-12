#pragma once

#include "../Concepts.h"
#include "Function.h"

#include <algorithm>
#include <ranges>

namespace Interpolation {

namespace OneDimensional {

template <RealOrComplexFloatingPoint _Scalar> class Polynomial;

namespace Internal {

template <RealOrComplexFloatingPoint _Scalar>
struct Traits<Polynomial<_Scalar>> {
    using Real = _Scalar;
    using Scalar = _Scalar;
};
}   // namespace Internal

template <RealOrComplexFloatingPoint _Scalar>
class Polynomial : public Function<Polynomial<_Scalar>> {
  public:
    using Real = typename Function<Polynomial<_Scalar>>::Real;
    using Scalar = typename Function<Polynomial<_Scalar>>::Scalar;

    template <RealOrComplexView View> Polynomial(View a) {
        std::ranges::copy(a, std::back_inserter(_a));
    }

    auto Order() const { return _a.size() - 1; }

  private:
    std::vector<Scalar> _a;
};

// Deduction guide to build from ranges.
template <RealOrComplexView View>
Polynomial(View) -> Polynomial<std::ranges::range_value_t<View>>;

template <RealOrComplexView View>

auto
PolynomialEvaluate(View a, std::ranges::range_value_t<View> x) {
    using Scalar = std::ranges::range_value_t<View>;
    return std::ranges::fold_right_last(
               a, [x](auto b, auto c) { return b + x * c; })
        .value_or(Scalar{0});
}

template <std::ranges::viewable_range Range>
auto
PolynomialEvaluate(Range &a, std::ranges::range_value_t<Range> x) {
    return PolynomialEvaluate(std::ranges::views::all(a), x);
}

}   // namespace OneDimensional

}   // namespace Interpolation
