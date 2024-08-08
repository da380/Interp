#ifndef INTERPOLATION_LINEAR_GUARD_H
#define INTERPOLATION_LINEAR_GUARD_H

#include "Concepts.h"
#include "Function1D.h"
#include <type_traits>

// #include <algorithm>
// #include <cassert>
// #include <concepts>
// #include <iterator>
// #include <vector>

namespace Interpolation {

template <RealFloatingPoint Real, RealOrComplexFloatingPoint Scalar>
class Linear;

namespace Internal {

template <RealFloatingPoint _Real, RealOrComplexFloatingPoint _Scalar>
struct Traits<Linear<_Real, _Scalar>> {
    using Real = _Real;
    using Scalar = _Scalar;
};

}   // namespace Internal

template <RealFloatingPoint _Real, RealOrComplexFloatingPoint _Scalar>
class Linear : public Function1D<Linear<_Real, _Scalar>> {
  public:
    using Real = _Real;
    using Scalar = _Scalar;

    // template <std::integral N> auto Evaluate<N>(Real x) const { return 0; }

    Linear() = default;
};

}   // namespace Interpolation

#endif   // INTERPOLATION_LINEAR_GUARD_H`