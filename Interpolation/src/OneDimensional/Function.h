#pragma once

#include "Traits.h"
#include <cassert>

namespace Interpolation {

namespace OneDimensional {

template <typename Derived> class Function {
  public:
    using Real = typename Internal::Traits<Derived>::Real;
    using Scalar = typename Internal::Traits<Derived>::Scalar;

    template <int N> auto Evaluate(Real x) const {
        return GetDerived().template Evaluate<N>(x);
    }

    auto operator()(Real x) const { return Evaluate<0>(x); }

  private:
    auto &GetDerived() const { return static_cast<const Derived &>(*this); }
    auto &GetDerived() { return static_cast<Derived &>(*this); }
};

}   // namespace OneDimensional

}   // namespace Interpolation
