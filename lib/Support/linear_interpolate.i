// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%import "auto_derivative.i"
%{
#include "linear_interpolate.h"

// LinearInterpolate is a bit difficult to wrap directly for
// python. So add a adapter class.
namespace FullPhysics {
template<class TX, class TY> class LinearInterpolateWrap :
  public LinearInterpolate<TX, TY> {
public:
    typedef typename LinearInterpolate<TX, TY>::BehaviorOutOfRange BehaviorOutOfRange;
    LinearInterpolateWrap
      (const std::vector<TX>& X,
       const std::vector<TY>& Y
       )
      : LinearInterpolate<TX, TY>
      (X.begin(), X.end(), Y.begin()) {}
    LinearInterpolateWrap
      (const std::vector<TX>& X,
       const std::vector<TY>& Y,
       int Out_of_range)
      : LinearInterpolate<TX, TY>
	(X.begin(), X.end(), Y.begin(), (BehaviorOutOfRange) Out_of_range) {}
    TY __call__(const TX& x) const
      { return (*this)(x); }
};
}
%}

%fp_shared_ptr(FullPhysics::LinearInterpolateWrap<FullPhysics::AutoDerivative<double>, FullPhysics::AutoDerivative<double> >);

namespace FullPhysics {
template<class TX, class TY> class LinearInterpolateWrap {
public:
    enum {OUT_OF_RANGE_EXTRAPOLATE = 0, OUT_OF_RANGE_CLIP, OUT_OF_RANGE_ERROR};
    LinearInterpolateWrap
      (const std::vector<TX>& X,
       const std::vector<TY>& Y);
    LinearInterpolateWrap
      (const std::vector<TX>& X,
       const std::vector<TY>& Y,
       int Out_of_range);
    TY __call__(const TX& x) const;
};

%template(LinearInterpolateAutoDerivative) FullPhysics::LinearInterpolateWrap<FullPhysics::AutoDerivative<double>, FullPhysics::AutoDerivative<double> >;
}




