#ifndef INTERPOLATION_FUNCTION_1D_GUARD_H
#define INTERPOLATION_FUNCTION_1D_GUARD_H

#include "Traits.h"
#include <cassert>

namespace Interpolation {

template <typename Derived> class Function1D {
  public:
    using Real = typename Internal::Traits<Derived>::Real;
    using Scalar = typename Internal::Traits<Derived>::Scalar;

    /*
    template <std::integral N> auto Evaluate(Real x) const {
        return GetDerived().Evaluate<N>(x);
    }
*/

    auto operator()(Real x) const { return 0; }

  private:
    auto &GetDerived() const { return static_cast<const Derived &>(*this); }
    auto &GetDerived() { return static_cast<Derived &>(*this); }
};

}   // namespace Interpolation

#endif   // INTERPOLATION_BASE_1D_GUARD_H