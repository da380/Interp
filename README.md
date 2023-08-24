# Interp


c++ header-only template library for iterpolation. Linear algebra done using Eigen3. 

Cubic splines interpolation done using natural boundary conditions, ie that the second derivative is zero at the start and end points. 
Checks performed in CubicSplineExample.cpp given in ./examples/:
1. Linear function y = x, all results were as expected.
2. y = x^2. Tests all three methods.
3. y = x^3 - x^2 + x. Tests all three methods.

Implemented so that clamped and natural boundary conditions can be used. Together if so desired, ie lower BC clamped, upper natural, or vice versa. This is done using constructor overloading. If one uses just three inputs it does fully natural cubic spline interpolation. 




Hi

