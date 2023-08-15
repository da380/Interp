#include <gtest/gtest.h>

#include "TestLinear.h"

// Tests for linear interpolation
TEST(Linear, CheckNodelValues) {
  int i = LinearCheckNodalValues<double>();
  EXPECT_EQ(0, i);
}
