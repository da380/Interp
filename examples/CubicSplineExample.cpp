
#include <Eigen/Dense>
#include <iostream>

#include "Interp/Core"



int main() {

  using complex = std::complex<double>;
  using Matrix = Eigen::Matrix<complex, Eigen::Dynamic, Eigen::Dynamic>;

  auto a = Matrix(5,5);

  std::cout << a << std::endl;
  

}
