
add_executable(LinearExample LinearExample.cpp)
target_link_libraries(LinearExample Interpolation)

add_executable(CubicSplineExample CubicSplineExample.cpp)
target_link_libraries(CubicSplineExample Interpolation)

add_executable(LagrangeExample LagrangeExample.cpp)
target_link_libraries(LagrangeExample Interpolation)









