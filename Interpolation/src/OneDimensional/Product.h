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
class Product;

namespace Internal {

template <typename Derived0, typename Derived1>
struct Traits<Product<Derived0, Derived1>> {
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
class Product : public Function<Product<Derived0, Derived1>> {
  public:
    using Real = typename Function<Product<Derived0, Derived1>>::Real;
    using Scalar = typename Function<Product<Derived0, Derived1>>::Scalar;

    Product(Derived0 &f0, Derived1 &f1) : _f0{f0}, _f1{f1} {}

    template <int N>
    auto Evaluate(Real x) const
        requires(N >= 0 && N <= 2)
    {
        if constexpr (N == 0) {
            return _f0.template Evaluate<0>(x) * _f1.template Evaluate<0>(x);
        } else if constexpr (N == 1) {
            return _f0.template Evaluate<1>(x) * _f1.template Evaluate<0>(x) +
                   _f0.template Evaluate<0>(x) * _f1.template Evaluate<1>(x);
        } else if constexpr (N == 2) {
            return _f0.template Evaluate<2>(x) * _f1.template Evaluate<0>(x) +
                   Scalar{2} * _f0.template Evaluate<1>(x) *
                       _f1.template Evaluate<1>(x) +
                   _f0.template Evaluate<0>(x) * _f1.template Evaluate<2>(x);
        }
    }

  private:
    Derived0 &_f0;
    Derived1 &_f1;
};

}   // namespace OneDimensional

}   // namespace Interpolation
