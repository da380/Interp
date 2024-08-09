#pragma once

#include "../Concepts.h"
#include "Function.h"

#include <type_traits>

namespace Interpolation {

namespace OneDimensional {

template <typename Derived0, typename Derived1>
    requires requires() {
        requires std::derived_from<Derived0, Function<Derived0>>;
        requires std::derived_from<Derived1, Function<Derived1>>;
        requires std::same_as<typename Derived0::Real, typename Derived1::Real>;
        requires std::convertible_to<typename Derived0::Scalar,
                                     typename Derived1::Scalar> ||
                     std::convertible_to<typename Derived1::Scalar,
                                         typename Derived0::Scalar>;
    }
class Sum;

namespace Internal {

template <typename Derived0, typename Derived1>
struct Traits<Sum<Derived0, Derived1>> {
    using Real = typename Derived0::Real;
    using Scalar = std::invoke_result<std::plus<>, typename Derived0::Scalar,
                                      typename Derived1::Scalar>;
};

}   // namespace Internal

template <typename Derived0, typename Derived1>
    requires requires() {
        requires std::derived_from<Derived0, Function<Derived0>>;
        requires std::derived_from<Derived1, Function<Derived1>>;
        requires std::same_as<typename Derived0::Real, typename Derived1::Real>;
        requires std::convertible_to<typename Derived0::Scalar,
                                     typename Derived1::Scalar> ||
                     std::convertible_to<typename Derived1::Scalar,
                                         typename Derived0::Scalar>;
    }
class Sum : public Function<Sum<Derived0, Derived1>> {
  public:
    using Real = typename Function<Sum<Derived0, Derived1>>::Real;
    using Scalar = typename Function<Sum<Derived0, Derived1>>::Scalar;

    Sum(Derived0 &f0, Derived1 &f1) : _f0{f0}, _f1{f1} {}

    template <int N> auto Evaluate(Real x) const {
        return _f0.template Evaluate<N>(x) + _f1.template Evaluate<N>(x);
    }

  private:
    Derived0 &_f0;
    Derived1 &_f1;
};

}   // namespace OneDimensional

}   // namespace Interpolation
