#include <Interpolation/All>
#include <iostream>

using namespace Interpolation;

int
main() {

    using Real = double;

    auto f = Linear<Real, Real>();

    std::cout << f(2) << std::endl;
}