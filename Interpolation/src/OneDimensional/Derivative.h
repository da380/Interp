#pragma once

#include "../Concepts.h"
#include "Function.h"

namespace Interpolation {

namespace OneDimensional {

template <typename Derived>
    requires std::derived_from<Derived, Function<Derived>>
class Derivative;

namespace Internal {

template <typename Derived> struct Traits<Derivative<Derived>> {
    using Real = typename Derived::Real;
    using Scalar = typename Derived::Scalar;
};

}   // namespace Internal

template <typename Derived>
    requires std::derived_from<Derived, Function<Derived>>
class Derivative : public Function<Derivative<Derived>> {
  public:
    using Real = typename Function<Derivative<Derived>>::Real;
    using Scalar = typename Function<Derivative<Derived>>::Scalar;

    Derivative(Derived &f) : _f{f} {}

    template <int N> auto Evaluate(Real x) const {
        return _f.template Evaluate<N + 1>(x);
    }

  private:
    Derived &_f;
};

}   // namespace OneDimensional

}   // namespace Interpolation
