#ifndef INTERPOLATION_CONCEPTS_GUARD_H
#define INTERPOLATION_CONCEPTS_GUARD_H

#include <complex>
#include <concepts>
#include <iterator>
#include <type_traits>

namespace Interpolation {

// Concepts for real or complex floating point types.
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

// Concepts for iterators with real or complex floating point values.
template <typename T>
concept RealFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires RealFloatingPoint<std::iter_value_t<T>>;
};

template <typename T>
concept ComplexFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires ComplexFloatingPoint<std::iter_value_t<T>>;
};

template <typename T>
concept RealOrComplexFloatingPointIterator = requires() {
    requires std::random_access_iterator<T>;
    requires RealOrComplexFloatingPoint<std::iter_value_t<T>>;
};

template <typename xIter, typename yIter>
concept InterpolationIteratorPair = requires(xIter x, yIter y) {
    requires RealFloatingPointIterator<xIter>;
    requires RealOrComplexFloatingPointIterator<yIter>;
    requires std::convertible_to<std::iter_value_t<xIter>,
                                 std::iter_value_t<yIter>>;
    { (*x) + (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*x) * (*y) } -> std::convertible_to<std::iter_value_t<yIter>>;
    { (*y) / (*x) } -> std::convertible_to<std::iter_value_t<yIter>>;
};

}   // namespace Interpolation

#endif   //  INTERPOLATION_CONCEPTS_GUARD_H
