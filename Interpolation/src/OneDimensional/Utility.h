#pragma once

#include "../Concepts.h"
#include <algorithm>
#include <ranges>

namespace Interpolation {

namespace OneDimensional {

// Locate point within sorted 1D array.
template <RealView Abscissa,
          RealFloatingPoint Real = std::ranges::range_value_t<Abscissa>>
auto
Find(Abscissa x, Real xp) {
    auto iter = std::ranges::upper_bound(x, xp);
    if (iter == x.begin())
        ++iter;
    if (iter == x.end())
        --iter;
    return iter;
}
}   // namespace OneDimensional
}   // namespace Interpolation
