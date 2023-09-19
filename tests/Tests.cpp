#include <gtest/gtest.h>

#include "TestCubicSpline.h"
#include "TestLinear.h"

// Tests for linear interpolation

TEST(Linear, CheckRealSingle) {
  int i = LinearCheck<float, float>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckRealDouble) {
  int i = LinearCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckRealLongDouble) {
  int i = LinearCheck<long double, long double>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexSingle) {
  int i = LinearCheck<float, std::complex<float>>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexDouble) {
  int i = LinearCheck<double, std::complex<double>>();
  EXPECT_EQ(0, i);
}

TEST(Linear, CheckComplexLongDouble) {
  int i = LinearCheck<long double, std::complex<long double>>();
  EXPECT_EQ(0, i);
}

// Tests for cubic spline interpolation

TEST(CubicSpline, CheckRealSingle) {
  int i = CubicSplineCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckRealDouble) {
  int i = CubicSplineCheck<double, double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckRealLongDouble) {
  int i = CubicSplineCheck<long double, long double>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexSingle) {
  int i = CubicSplineCheck<float, std::complex<float>>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexDouble) {
  int i = CubicSplineCheck<double, std::complex<double>>();
  EXPECT_EQ(0, i);
}

TEST(CubicSpline, CheckComplexLongDouble) {
  int i = CubicSplineCheck<long double, std::complex<long double>>();
  EXPECT_EQ(0, i);
}
