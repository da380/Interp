#pragma once

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>

namespace Interpolation {

template <typename T> struct RemoveComplexHelper {
    using value_type = T;
};

template <typename T> struct RemoveComplexHelper<std::complex<T>> {
    using value_type = T;
};

template <typename T>
using RemoveComplex = typename RemoveComplexHelper<T>::value_type;
template <typename T> struct IsComplexFloatingPoint : public std::false_type {};

template <typename T>
struct IsComplexFloatingPoint<std::complex<T>>
    : public std::bool_constant<std::is_floating_point_v<T>> {};

template <typename T>
concept RealFloatingPoint = std::floating_point<T>;

template <typename T>
concept ComplexFloatingPoint =
    IsComplexFloatingPoint<std::remove_const_t<T>>::value;

template <typename T>
concept RealOrComplexFloatingPoint =
    RealFloatingPoint<T> or ComplexFloatingPoint<T>;

template <typename T>
concept RealView = requires() {
    requires std::ranges::view<T>;
    requires std::ranges::random_access_range<T>;
    requires RealFloatingPoint<std::ranges::range_value_t<T>>;
};

template <typename T>
concept ComplexView = requires() {
    requires std::ranges::view<T>;
    requires std::ranges::random_access_range<T>;
    requires ComplexFloatingPoint<std::ranges::range_value_t<T>>;
};

template <typename T>
concept RealOrComplexView = RealView<T> || ComplexView<T>;

template <typename Abscissa, typename Ordinate>
concept DataViews1D = requires() {
    requires RealView<Abscissa>;
    requires RealOrComplexView<Ordinate>;
    requires std::same_as<std::ranges::range_value_t<Abscissa>,
                          RemoveComplex<std::ranges::range_value_t<Ordinate>>>;
};
}   // namespace Interpolation
