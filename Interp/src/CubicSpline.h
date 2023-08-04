#ifndef INTERP_CUBIC_SPLINE_GUARD_H
#define INTERP_CUBIC_SPLINE_GUARD_H

namespace Interp {

template <typename OrdIter, typename CordIter>
class CubicSpline {
 public:
  using ordinate_value_type = OrdIter::value_type;
  using coordinate_value_type = CordIter::value_type;

 private:
};

}  // namespace Interp

#endif  //  INTERP_CUBIC_SPLINE_GUARD_H
