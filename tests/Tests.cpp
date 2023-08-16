#include <gtest/gtest.h>

#include "TestLinear.h"

// Tests for linear interpolation
TEST(Linear, CheckNodalValues) {
  int i = LinearCheckNodalValues<double>();
  EXPECT_EQ(0, i);
}
