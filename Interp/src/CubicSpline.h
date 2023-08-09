#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

#include <algorithm>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <optional>
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

     // Declare the constructor with same at both ends
     CubicSpline(xIter, xIter, yIter, std::string, y_value_type, y_value_type);

     // Declare the constructor without final three arguments
     CubicSpline(xIter, xIter, yIter);

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

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
     CubicSpline(xStart, xFinish, yStart, "natural", 0.0, "natural", 0.0);
};

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, std::string bcarg,
                                       y_value_type lowbc, y_value_type upbc)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
     CubicSpline(xStart, xFinish, yStart, bcarg, lowbc, bcarg, upbc);
};

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::CubicSpline(xIter xStart, xIter xFinish,
                                       yIter yStart, std::string lbcarg,
                                       y_value_type lowbc, std::string ubcarg,
                                       y_value_type upbc)
    : xStart{xStart}, xFinish{xFinish}, yStart{yStart} {
     xIter idx;
     int i, j, mylen, ilen;
     const int mylen2 = xFinish - xStart;
     // int matlength = std::distance(xStart, xFinish) - 2;

     // for (idx = xStart; idx < xFinish; idx++) {
     // deltah.push_back(*(idx + 1) - *idx);
     // }
     i = 0;
     j = 0;
     // mylen = 0;
     ilen = xFinish - xStart;
     for (idx = xStart + 1; idx < xFinish - 1; idx++) {
          j += 1;
     };
     mylen = j;
     j = 0;
     Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> mymatM =
         Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic>::Zero(
             mylen, mylen);
     Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> myRHS =
         Eigen::Matrix<x_value_type, Eigen::Dynamic, 1>::Zero(mylen, 1);
     //  std::cout << mylen << std::endl;
     //  std::cout << mymatM << std::endl;
     //  std::cout << mymatM(0, 0) << std::endl;
     //  mymatM(mylen - 1, mylen - 1) = 5;
     //  std::cout << mymatM << std::endl;

     for (i = 0; i < mylen; i++) {
          mymatM(i, i) = 2;
     };
     i = 0;
     //  std::cout << "Hello, world 1!" << std::endl;
     for (idx = xStart + 1; idx < xFinish - 2; idx++) {
          mymatM(i, i + 1) = (*(idx + 1) - *(idx)) / (*(idx + 1) - *(idx - 1));
          i += 1;
     }
     //  std::cout << "Hello, world 2!" << std::endl;
     for (i = 0; i < mylen - 2; i++) {
          mymatM(i + 1, i) = 1 - mymatM(i + 1, i + 2);
     }
     //  std::cout << "Hello, world 3!" << std::endl;
     mymatM(mylen - 1, mylen - 2) =
         (*(xFinish - 2) - *(xFinish - 3)) / (*(xFinish - 1) - *(xFinish - 3));
     //  std::cout << mymatM << std::endl;

     // finding the rhs
     i = 0;
     x_value_type mu;
     for (idx = xStart; idx < xFinish - 2; idx++) {
          // mu = (*(idx + 1) - *idx) / (*(idx + 2) - *idx);
          // // std::cout << mu << std::endl;
          // // std::cout << *(idx + 2) << std::endl;
          // myRHS(i) = 6.0 *
          //            (*(idx) + (mu * (*(idx + 2)) - *(idx + 1)) /
          //            (1.0 - mu))
          //            * mu / ((*(idx + 1) - *idx) * (*(idx + 1) -
          //            *idx));
          // std::cout << (*(idx) + (mu * (*(idx + 2)) - *(idx + 1)) /
          // (1 - mu))
          // << std::endl;
          // std::cout << *idx + (mu * (*(idx + 2)) - *(idx + 1)) /
          // (1.0 - mu)
          //           << std::endl;

          mu = (*(xStart + i + 1) - *(xStart + i)) /
               (*(xStart + i + 2) - *(xStart + i));
          // std::cout << mu << std::endl;
          myRHS(i) =
              6.0 *
              (*(yStart + i) +
               (mu * (*(yStart + i + 2)) - *(yStart + i + 1)) / (1.0 - mu)) *
              mu /
              ((*(xStart + i + 1) - *(xStart + i)) *
               (*(xStart + i + 1) - *(xStart + i)));

          i += 1;
     };

     // print check
     //  std::cout << mymatM << std::endl;
     //  std::cout << myRHS << std::endl;

     // if the requirement is for clamped bcs:

     // finding ypp
     ypp.resize(mylen + 2);
     // const int mylenlu = ypp.size();
     // std::cout << ypp.size() << std::endl;
     // std::cout << myRHS.size() << std::endl;
     // std::cout << mymatM.size() << std::endl;

     // std::cout << xFinish - xStart << std::endl;
     //  ypp = mymatM.partialPivLu().solve(myRHS);

     // checking what type of cubic interpolant it wants:
     if (!lbcarg.compare("natural")) {
          // std::cout << "Cubic spline calculated with natural BC \n";
          ypp(0) = lowbc;
          myRHS(0) = myRHS(0) - (1.0 - mymatM(0, 1)) * lowbc;

     } else if (!lbcarg.compare("clamped")) {
          // std::cout << "Cubic spline calculated with clamped BC \n";
          mymatM(0, 0) = 1.5 + 0.5 * mymatM(0, 1);

          myRHS(0) = myRHS(0) + 3 * (1 - mymatM(0, 1)) /
                                    (*(xStart + 1) - *xStart) *
                                    (lowbc - (*(yStart + 1) - *yStart) /
                                                 (*(xStart + 1) - *xStart));

     } else {
          std::cout << "Please check lower BC! Cubic spline calculated with "
                       "natural BC \n";
     }
     if (!ubcarg.compare("natural")) {
          // std::cout << "Cubic spline calculated with natural BC \n";

          ypp(mylen + 1) = upbc;
          myRHS(mylen - 1) =
              myRHS(mylen - 1) - (1.0 - mymatM(mylen - 1, mylen - 2)) * upbc;

     } else if (!ubcarg.compare("clamped")) {
          // std::cout << "Cubic spline calculated with clamped BC \n";
          // mymatM(0, 0) = 1.5 + 0.5 * mymatM(0, 1);
          mymatM(mylen - 1, mylen - 1) =
              1.5 + 0.5 * mymatM(mylen - 1, mylen - 2);
          // myRHS(0) = myRHS(0) + 3 * (1 - mymatM(0, 1)) /
          //                           (*(xStart + 1) - *xStart) *
          //                           (lowbc - (*(yStart + 1) - *yStart) /
          //                                        (*(xStart + 1) - *xStart));
          myRHS(mylen - 1) =
              myRHS(mylen - 1) +
              3 * (1 - mymatM(mylen - 1, mylen - 2)) /
                  (*(xFinish - 1) - *(xFinish - 2)) *
                  ((*(yStart + ilen - 1) - *(yStart + ilen - 2)) /
                       (*(xFinish - 1) - *(xFinish - 2)) -
                   upbc);
     } else {
          std::cout << "Please check BC! Cubic spline calculated with "
                       "natural BC \n";
     }

     Eigen::FullPivLU<
         Eigen::Matrix<x_value_type, Eigen::Dynamic, Eigen::Dynamic> >
         dec(mymatM);
     // std::cout << "Hello, testing" << std::endl;
     Eigen::Matrix<y_value_type, Eigen::Dynamic, 1> ypps = dec.solve(myRHS);

     // for (i = 0; i < mylen; i++) {
     //      ypp(i + 1) = ypps(i);
     // }
     ypp.segment(1, mylen) = ypps;
     // std::cout << "Hello, testing 2" << std::endl;
     // std::cout << mymatM << std::endl;
     // std::cout << myRHS << std::endl;

     if (!lbcarg.compare("clamped")) {
          ypp(0) = -0.5 * ypp(1) +
                   3.0 / (*(xStart + 1) - *xStart) *
                       ((*(yStart + 1) - *yStart) / (*(xStart + 1) - *xStart) -
                        lowbc);
     }
     if (!ubcarg.compare("clamped")) {
          ypp(mylen + 1) =
              -0.5 * ypp(mylen) +
              3.0 / (*(xFinish - 1) - *(xFinish - 2)) *
                  (upbc - (*(yStart + ilen - 1) - *(yStart + ilen - 2)) /
                              (*(xFinish - 1) - *(xFinish - 2)));
     }
     // std::cout << ypp << std::endl;
};

// template <typename xIter, typename yIter>
// CubicSpline<xIter, yIter>::y_value_type CubicSpline<xIter,
// yIter>::operator()(
//     x_value_type x) {
//   auto iter = std::upper_bound(xStart, xFinish, x);
//   if (iter == xStart) {
//     return *yStart;
//   }
//   auto i2 = std::distance(xStart, iter);
//   auto i1 = i2 - 1;
//   if (iter == xFinish) {
//     return *std::next(yStart, i1);
//   } else {
//     auto x1 = *std::prev(iter);
//     auto x2 = *iter;
//     auto y1 = *std::next(yStart, i1);
//     auto y2 = *std::next(yStart, i2);
//     return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
//   }
// }

template <typename xIter, typename yIter>
CubicSpline<xIter, yIter>::y_value_type
CubicSpline<xIter, yIter>::operator()(x_value_type x) {
     auto iter = std::upper_bound(xStart, xFinish, x);
     // std::cout << iter - xStart << std::endl;
     if (iter == xStart) {
          std::cout << "The x value lies below the domain!" << std::endl;
          return *yStart;
     }
     if (x > *(xFinish - 1)) {
          std::cout << "The x value lies above the domain!" << std::endl;
          auto j = xFinish - xStart;
          return *(yStart + j - 1);
     }
     auto i = iter - xStart - 1;
     auto a = *(yStart + i);
     // std::cout << "Hello, testing 3" << std::endl;
     auto b = (*(yStart + i + 1) - *(yStart + i)) / (*iter - *(iter - 1)) -
              (ypp(i + 1) + 2.0 * ypp(i)) * (*(iter) - *(iter - 1)) / 6.0;
     // std::cout << "Hello, testing 4" << std::endl;
     auto c = ypp(i) / 2.0;
     auto d = (ypp(i + 1) - ypp(i)) / (6.0 * (*iter - *(iter - 1)));
     auto xdiff = x - *(iter - 1);
     auto yval = a + xdiff * (b + xdiff * (c + d * xdiff));
     return yval;
};

}   // namespace Interp

#endif   //  INTERP_CUBIC_SPLINE_GUARD_H
