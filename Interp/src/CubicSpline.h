#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <algorithm>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <string>
#include <vector>

namespace Interp {

template <typename xIter, typename yIter>
class CubicSpline {
    public:
     // Define some class member types
     using x_value_type = xIter::value_type;
     using y_value_type = yIter::value_type;

     // Declare the full constructor.
     CubicSpline(xIter, xIter, yIter, std::string, y_value_type, std::string,
                 y_value_type);

     // Natural BC constructor
     CubicSpline(xIter, xIter, yIter);

     // Both BC same constructor
     CubicSpline(xIter, xIter, yIter, std::string, y_value_type, y_value_type);

     // Declare the application operator.
     y_value_type operator()(x_value_type);

    private:
     using VectorX = Eigen::Matrix<x_value_type, Eigen::Dynamic, 1>;
     using VectorY = Eigen::Matrix<y_value_type, Eigen::Dynamic, 1>;
     // Iterators to the function data.
     xIter xStart;
     xIter xFinish;
     yIter yStart;

     // Cubic spline coefficients
     VectorY ypp;
};

// Natural BC constructor
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart)
    : CubicSpline(xStart, xFinish, yStart, "natural", 0, "natural", 0){};
template <typename xIter, typename yIter>

// Both BC same constructor
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, std::string bcarg,
                                       y_value_type lowbc, y_value_type upbc)
    : CubicSpline(xStart, xFinish, yStart, bcarg, lowbc, bcarg, upbc){};

// Full constructor
template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, std::string lbcarg,
                                       y_value_type lowbc, std::string ubcarg,
                                       y_value_type upbc)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
     /////////////////////////////////////////////////////////////////////

     // Declarations
     int idxint;             // integer index
     const int mylen =
         xFinish - xStart;   // difference between start and finish

     // Using Eigen's inbuilt matrix type
     // Matrix M
     Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> mymatM =
         Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic>::Zero(
             mylen - 2, mylen - 2);
     // RHS
     Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> myRHS =
         Eigen::Matrix<x_value_type, Eigen::Dynamic, 1>::Zero(mylen - 2, 1);

     /////////////////////////////////////////////////////////////////////
     // filling out mymatM
     // main diagonal
     for (idxint = 0; idxint < mylen - 2; idxint++) {
          mymatM(idxint, idxint) = 2;
     };

     // right diagonal
     for (idxint = 0; idxint < mylen - 3; idxint++) {
          mymatM(idxint, idxint + 1) =
              (*(xStart + idxint + 2) - *(xStart + idxint + 1)) /
              (*(xStart + idxint + 2) - *(xStart + idxint));
     }

     // left diagonal
     for (idxint = 0; idxint < mylen - 4; idxint++) {
          mymatM(idxint + 1, idxint) = 1 - mymatM(idxint + 1, idxint + 2);
     }

     // final element of left diagonal
     mymatM(mylen - 3, mylen - 4) =
         (*(xFinish - 2) - *(xFinish - 3)) / (*(xFinish - 1) - *(xFinish - 3));

     /////////////////////////////////////////////////////////////////////
     // filling out RHS
     x_value_type mu;
     for (idxint = 0; idxint < mylen - 2; idxint++) {
          mu = (*(xStart + idxint + 1) - *(xStart + idxint)) /
               (*(xStart + idxint + 2) - *(xStart + idxint));

          myRHS(idxint) = 6.0 *
                          (*(yStart + idxint) + (mu * (*(yStart + idxint + 2)) -
                                                 *(yStart + idxint + 1)) /
                                                    (1.0 - mu)) *
                          mu /
                          ((*(xStart + idxint + 1) - *(xStart + idxint)) *
                           (*(xStart + idxint + 1) - *(xStart + idxint)));
     };

     /////////////////////////////////////////////////////////////////////
     // making small changes depending upon boundary conditions
     //  size ypp
     ypp.resize(mylen);

     // lower BC
     if (!lbcarg.compare("natural")) {
          // M doesn't change, only RHS does
          // set ypp_0 to be the inputted second derivative
          ypp(0) = lowbc;

          // correct RHS
          myRHS(0) = myRHS(0) - (1.0 - mymatM(0, 1)) * lowbc;

     } else if (!lbcarg.compare("clamped")) {
          // note: can't find ypp(0) here, as depends on ypp(1)
          // correct M
          mymatM(0, 0) = 1.5 + 0.5 * mymatM(0, 1);

          // correct RHS(0)
          myRHS(0) = myRHS(0) + 3 * (1 - mymatM(0, 1)) /
                                    (*(xStart + 1) - *xStart) *
                                    (lowbc - (*(yStart + 1) - *yStart) /
                                                 (*(xStart + 1) - *xStart));

     } else {
          std::cout << "Please check lower BC! Cubic spline calculated with "
                       "default natural BC \n";
     }

     // upper BC
     if (!ubcarg.compare("natural")) {
          // M doesn't change, only RHS does
          // set ypp(end) to be the inputted second derivative
          ypp(mylen - 1) = upbc;

          // correct RHS
          myRHS(mylen - 3) =
              myRHS(mylen - 3) - (1.0 - mymatM(mylen - 3, mylen - 4)) * upbc;

     } else if (!ubcarg.compare("clamped")) {
          // note: can't find ypp(end) here
          // correct M
          mymatM(mylen - 3, mylen - 3) =
              1.5 + 0.5 * mymatM(mylen - 3, mylen - 4);
          // correct RHS(end)
          myRHS(mylen - 3) =
              myRHS(mylen - 3) +
              3 * (1 - mymatM(mylen - 3, mylen - 4)) /
                  (*(xFinish - 1) - *(xFinish - 2)) *
                  ((*(yStart + mylen - 1) - *(yStart + mylen - 2)) /
                       (*(xFinish - 1) - *(xFinish - 2)) -
                   upbc);
     } else {
          std::cout << "Please check BC! Cubic spline calculated with "
                       "natural BC \n";
     }

     /////////////////////////////////////////////////////////////////////
     // LU decomposition
     // prepare for decomposition by declaring class dec. NOTE: sizing is done
     // using Eigen::Dynamic, as not known at compilation
     Eigen::FullPivLU<
         Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> >
         dec(mymatM);

     // actually doing LU decomposition
     Eigen::Matrix<y_value_type, Eigen::Dynamic, 1> ypps = dec.solve(myRHS);

     /////////////////////////////////////////////////////////////////////
     // finalising the second derivative
     // assigning middle segment
     ypp.segment(1, mylen - 2) = ypps;

     // correcting for clamped BCs
     if (!lbcarg.compare("clamped")) {
          ypp(0) = -0.5 * ypp(1) +
                   3.0 / (*(xStart + 1) - *xStart) *
                       ((*(yStart + 1) - *yStart) / (*(xStart + 1) - *xStart) -
                        lowbc);
     }
     if (!ubcarg.compare("clamped")) {
          ypp(mylen - 1) =
              -0.5 * ypp(mylen - 2) +
              3.0 / (*(xFinish - 1) - *(xFinish - 2)) *
                  (upbc - (*(yStart + mylen - 1) - *(yStart + mylen - 2)) /
                              (*(xFinish - 1) - *(xFinish - 2)));
     }

     // test, commented out
     // std::cout << ypp << std::endl;
};

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::y_value_type
CubicSpline<xIter, yIter>::operator()(x_value_type x) {
     // finding first element larger than x
     auto iter = std::upper_bound(xStart, xFinish, x);

     // error messages
     if (iter == xStart) {
          std::cout << "The x value lies below the domain!" << std::endl;
          return *yStart;
     }
     if (x > *(xFinish - 1)) {
          std::cout << "The x value lies above the domain!" << std::endl;
          auto j = xFinish - xStart;
          return *(yStart + j - 1);
     }

     // finding expansion coefficients
     auto i = iter - xStart - 1;
     auto a = *(yStart + i);
     auto b = (*(yStart + i + 1) - *(yStart + i)) / (*iter - *(iter - 1)) -
              (CubicSpline::ypp(i + 1) + 2.0 * ypp(i)) *
                  (*(iter) - *(iter - 1)) / 6.0;
     auto c = ypp(i) / 2.0;
     auto d = (ypp(i + 1) - ypp(i)) / (6.0 * (*iter - *(iter - 1)));

     // (x - x_{i - 1})
     auto xdiff = x - *(iter - 1);

     // finding value
     auto yval = a + xdiff * (b + xdiff * (c + d * xdiff));
     return yval;
};

}   // namespace Interp

#endif   //  INTERP_CUBIC_SPLINE_GUARD_H
