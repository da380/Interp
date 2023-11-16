#ifndef INTERPOLATION_POLYNOMIAL_GUARD_H
#define INTERPOLATION_POLYNOMIAL_GUARD_H

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <random>
#include <vector>

#include "Concepts.h"

namespace Interpolation {

template <typename T>
requires RealOrComplexFloatingPoint<T>
class Polynomial1D {
 public:
  // Type alias for the scalar.
  using value_type = T;

  // Construct from std::vector.
  Polynomial1D(std::vector<T> a) : _a{a} {}

  // Construct from std::initializer list.
  Polynomial1D(std::initializer_list<T> list) : _a{std::vector<T>{list}} {}

  // Return a real random polynomial of given degree.
  static Polynomial1D Random(int n)
  requires RealFloatingPoint<T>
  {
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<T> d{};
    std::vector<T> a;
    std::generate_n(std::back_inserter(a), n + 1, [&]() {
      if constexpr (RealFloatingPoint<T>) {
        return d(gen);
      }
    });
    return Polynomial1D(a);
  }

  // Return a complex random polynomial of given degree.
  static Polynomial1D Random(int n)
  requires ComplexFloatingPoint<T>
  {
    using S = typename T::value_type;
    std::random_device rd{};
    std::mt19937_64 gen{rd()};
    std::normal_distribution<S> d{};
    std::vector<T> a;
    std::generate_n(std::back_inserter(a), n + 1, [&]() {
      return T{d(gen), d(gen)};
    });
    return Polynomial1D(a);
  }

  // Returns the degree.
  auto Degree() const { return _a.size() - 1; }

  // Evaluates the function.
  T operator()(T x) const {
    return std::accumulate(_a.rbegin(), _a.rend(), static_cast<T>(0),
                           [x](auto p, auto a) { return p * x + a; });
  }

  // Evaluates the derivative.
  T Derivative(T x) const {
    return std::accumulate(_a.rbegin(), std::prev(_a.rend()), static_cast<T>(0),
                           [x, m = Degree()](auto p, auto a) mutable {
                             return p * x + static_cast<T>(m--) * a;
                           });
  }

  // Evaluates the primative.
  T Primative(T x) const {
    return std::accumulate(_a.rbegin(), _a.rend(), static_cast<T>(0),
                           [x, m = Degree() + 1](auto p, auto a) mutable {
                             return p * x + a * x / static_cast<T>(m--);
                           });
  }

  // Returns integral over [a,b].
  T Integrate(T a, T b) const { return Primative(b) - Primative(a); }

//output 
std::vector<T> polycoeff() {
  return this->_a;
};


//addition operator
  template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> operator+(FLOAT b) {
    std::vector<T> myval = _a;
    myval[0] += b;
    Polynomial1D<T> res{myval};
    return res;
  };

//subtraction operator
  template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> operator-(FLOAT b) {
    std::vector<T> myval = _a;
    myval[0] -= b;
    Polynomial1D<T> res{myval};
    return res;
  };

  // multiplication operator
  template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> operator*(FLOAT b) {
    std::vector<T> myval;
    typename std::vector<T>::iterator iter = _a.begin();
    for (iter; iter < _a.end(); ++iter){
      myval.push_back(b * iter[0]);
    }
    Polynomial1D<T> res{myval};
    return res;
  };

// division operator
template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> operator/(FLOAT b) {
    std::vector<T> myval;
    typename std::vector<T>::iterator iter = _a.begin();
    for (iter; iter < _a.end(); ++iter){
      myval.push_back(iter[0]/b);
    }
    Polynomial1D<T> res{myval};
    return res;
  };

// template<typename FLOAT> 
//   requires RealOrComplexFloatingPoint<FLOAT>
//   Polynomial1D<T> operator+(Polynomial1D<FLOAT>  b) {
//     auto maxnum = std::max(this->Degree(), b.Degree());
//     auto minnum = std::min(this->Degree(), b.Degree());
//     std::vector<T> myval;
//     for (int idx = 0; idx < maxnum+1; ++idx){
//       if (idx < this->Degree() && idx < b.Degree()){
//         myval.push_back(_a[idx]+b.polycoeff()[idx]);
//       } else if (idx > this->Degree() && idx < b.Degree()){
//         myval.push_back(b.polycoeff()[idx]);
//       } else if (idx < this->Degree() && idx > b.Degree()){
//         myval.push_back(_a[idx]);
//     }}
//     // myval[0] += b;
//     Polynomial1D<T> res{myval};
//     return res;
//   };


  

 private:
  std::vector<T> _a;  // Vector of polynomial coefficients.
};

}  // namespace Interpolation

#endif  // INTERPOLATION_POLYNOMIAL_GUARD_H
