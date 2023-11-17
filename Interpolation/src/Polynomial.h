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
  Polynomial1D(std::vector<T>& a) : _a{a} {}
   Polynomial1D(std::vector<T>&& a) : _a{std::move(a)} {}

  // Construct from std::initializer list.
  Polynomial1D(std::initializer_list<T> list) : _a{std::vector<T>{list}} {}

//
  
   Polynomial1D(const Polynomial1D&) = default;
   Polynomial1D(Polynomial1D&&) = default;


  //defining equals
  template<typename FLOAT>
  requires std::is_convertible_v<FLOAT,T>
  Polynomial1D(Polynomial1D<FLOAT> polinit) {
    // this->_a = polinit.polycoeff();
    std::vector<FLOAT> mtmp = polinit.polycoeff();
    std::transform(mtmp.begin(), mtmp.end(), std::back_inserter(this->_a), [](double val) -> float { return (float)val; });
    // _a(mtmp.begin(), mtmp.end());
  }
  template<typename FLOAT>
  requires std::is_convertible_v<FLOAT,T>
  Polynomial1D<T>& operator=(Polynomial1D<FLOAT>& polinit){
    this->_a = polinit.polycoeff();
    return *this;
  }
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
  } else if (i < 0) {
    return 0.0;
    } else{
return this->_a[i];
  }
  
};

// self operators, ie ()= operators
//addition 
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator+=(FLOAT  b) {
    _a[0]+=b;
    return *this;
    };

//subtraction 
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator-=(FLOAT  b) {
    _a[0]-=b;
    return *this;
    };
  

  // multiplication 
  template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator*=(FLOAT  b) {
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] *= b;
    };
    return *this;
    };
  

// division 
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator/=(FLOAT  b) {
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] /= b;
    };
    return *this;
    };


//+= operator with polynomial
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
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

//-= operator with polynomial
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator-=(const Polynomial1D<FLOAT>&  b) {

    bool sdeg = (this->Degree() < b.Degree());
    for (int idx = 0; idx < this->Degree() +1; ++idx){
      _a[idx] -= b.polycoeff(idx);
    };
    if (sdeg){
    for (int idx = this->Degree() +1; idx < b.Degree()+1; ++idx){
        _a.push_back(-b.polycoeff(idx));
    }; };
    return *this;
    };

//*= operator with polynomial
template<typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Polynomial1D<T> &operator*=(const Polynomial1D<FLOAT>&  b) {
    int maxc = this->Degree() + b.Degree();

     std::vector<T> newcoeff(maxc + 1,0.0);
    for (int idx = 0; idx < maxc + 1; ++idx){
      for (int idx2 = 0; idx2 < idx+1; ++idx2){
        newcoeff[idx] += this->polycoeff(idx2) * b.polycoeff(idx - idx2);
      }
    };
    this->_a = std::move(newcoeff);
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

//non-member operators
template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Interpolation::Polynomial1D<T> operator+(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval += b;
    return myval;
  };

  template<typename T,typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Interpolation::Polynomial1D<T> operator-(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval -= b;
    return myval;
  };

  template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
   Interpolation::Polynomial1D<T> operator*(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval *= b;
    return myval;
  };

template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
   Interpolation::Polynomial1D<T> operator/(Interpolation::Polynomial1D<T> a, FLOAT b) {
    Interpolation::Polynomial1D<T> myval = a;
    myval /= b;
    return myval;
  };

  template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Interpolation::Polynomial1D<T> operator+(const Interpolation::Polynomial1D<T>&  a, const Interpolation::Polynomial1D<FLOAT>&  b) {
    int maxval = std::max(a.Degree(), b.Degree());
    std::vector<T> myvec(maxval + 1, 0.0);
    for (int idx = 0; idx < maxval + 1; ++idx){
      myvec[idx] = a.polycoeff(idx) + b.polycoeff(idx);
    };
    Interpolation::Polynomial1D<T> myval{myvec};
    return myval;
    };

    template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Interpolation::Polynomial1D<T> operator-(const Interpolation::Polynomial1D<T>&  a, const Interpolation::Polynomial1D<FLOAT>&  b) {
    int maxval = std::max(a.Degree(), b.Degree());
    std::vector<T> myvec(maxval + 1, 0.0);
    for (int idx = 0; idx < maxval + 1; ++idx){
      myvec[idx] = a.polycoeff(idx) - b.polycoeff(idx);
    };
    Interpolation::Polynomial1D<T> myval{myvec};
    return myval;
    };

    template<typename T, typename FLOAT> 
  requires std::is_convertible_v<FLOAT, T>
  Interpolation::Polynomial1D<T> operator*(const Interpolation::Polynomial1D<T>&  a, const Interpolation::Polynomial1D<FLOAT>&  b){
    int maxc = a.Degree() + b.Degree();
     std::vector<T> newcoeff(maxc + 1,0.0);
    for (int idx = 0; idx < maxc + 1; ++idx){
      for (int idx2 = 0; idx2 < idx+1; ++idx2){
        newcoeff[idx] += a.polycoeff(idx2) * b.polycoeff(idx - idx2);
      }
    };
    Interpolation::Polynomial1D<T> myval{newcoeff};
    return myval;
    };

#endif  // INTERPOLATION_POLYNOMIAL_GUARD_H
