#pragma once

#include "../Concepts.h"
#include "Function.h"
#include "Utility.h"

#include <algorithm>
#include <cassert>

#include <ranges>

namespace Interpolation {

namespace OneDimensional {

template <RealView XView, RealOrComplexView YView>
    requires DataViews1D<XView, YView>
class Linear;

namespace Internal {

template <RealView XView, RealOrComplexView YView>
struct Traits<Linear<XView, YView>> {
    using Real = std::ranges::range_value_t<XView>;
    using Scalar = std::ranges::range_value_t<YView>;
};

}   // namespace Internal

template <RealView XView, RealOrComplexView YView>
    requires DataViews1D<XView, YView>
class Linear : public Function<Linear<XView, YView>> {

  public:
    using Real = typename Function<Linear<XView, YView>>::Real;
    using Scalar = typename Function<Linear<XView, YView>>::Scalar;

    Linear() = delete;

    Linear(XView x, YView y) : _x{x}, _y{y} {
        assert(_x.size() == _y.size());
        assert(std::ranges::is_sorted(_x));
    }

    template <int N>
    auto Evaluate(Real x) const
        requires(N == 0)
    {
        auto i2 = std::distance(_x.begin(), Find(_x, x));
        auto i1 = i2 - 1;
        auto x1 = _x[i1];
        auto x2 = _x[i2];
        auto h = x2 - x1;
        auto a = (x2 - x) / h;
        auto b = (x - x1) / h;
        return a * _y[i1] + b * _y[i2];
    }

  private:
    XView _x;
    YView _y;
};

// Deduction guides for construction from ranges.
template <std::ranges::viewable_range XView, std::ranges::viewable_range YView>
Linear(XView &&, YView &&) -> Linear<std::ranges::views::all_t<XView>,
                                     std::ranges::views::all_t<YView>>;

}   // namespace OneDimensional
}   // namespace Interpolation
