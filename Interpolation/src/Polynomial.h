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

//Constructor default
Polynomial1D() = default;

  // Construct from std::vector.
  Polynomial1D(std::vector<T> a) : _a{a} {}

  // Construct from std::initializer list.
  Polynomial1D(std::initializer_list<T> list) : _a{std::vector<T>{list}} {}

//
   Polynomial1D(const Polynomial1D&) = default;
   Polynomial1D(Polynomial1D&&) = default;

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

//output coefficient vector
std::vector<T> polycoeff() const {
  return this->_a;
};
T polycoeff(int i) const{
  if (i > this->Degree()){
    return 0.0;
  } else{
return this->_a[i];
  }
  
};


//addition operator
template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> &operator+=(FLOAT  b) {
    _a[0]+=b;
    return *this;
    };

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
  Polynomial1D<T> &operator-=(FLOAT  b) {
    _a[0]-=b;
    return *this;
    };
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
  Polynomial1D<T> &operator*=(FLOAT  b) {
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] *= b;
    };
    return *this;
    };
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
  Polynomial1D<T> &operator/=(FLOAT  b) {
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] /= b;
    };
    return *this;
    };
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

template<typename FLOAT> 
  requires RealOrComplexFloatingPoint<FLOAT>
  Polynomial1D<T> &operator+=(const Polynomial1D<FLOAT>&  b) {

    bool sdeg = (this->Degree() < b.Degree());
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] += b.polycoeff(idx);
    };
    if (sdeg){
    for (int idx = this->Degree() +1; idx < b.Degree()+1; ++idx){
        _a.push_back(b.polycoeff(idx));
    }; };
    return *this;
    };


//ostream
friend std::ostream& operator<<(std::ostream& os, const Polynomial1D<T>& obj)
{
  // Write obj to stream
  for (int idx = 0; idx < obj.Degree() + 1; ++idx){
    os << obj.polycoeff(idx) << " ";
  }
  return os;
};

  

 private:
  std::vector<T> _a;  // Vector of polynomial coefficients.
};



}  // namespace Interpolation

#endif  // INTERPOLATION_POLYNOMIAL_GUARD_H
