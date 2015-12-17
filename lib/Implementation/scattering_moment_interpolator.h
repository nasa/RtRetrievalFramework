#ifndef SCATTER_MOMENT_INTERPOLATOR_H
#define SCATTER_MOMENT_INTERPOLATOR_H
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include "auto_derivative.h"

namespace FullPhysics {
/****************************************************************//**
  The calculation of the phase function scattering matrix moments is
  fairly expensive. As an optimization, we interpolate between full
  calculated points to get the points in between. This class handles
  this interpolation.

  Normally, you don't directly use this class. Instead, Atmosphere
  uses this class to speed up its calculations.
*******************************************************************/

class ScatteringMomentInterpolate2Point {
public:
//-----------------------------------------------------------------------
/// Constructor, this takes the two wave numbers and scattering
/// matrices to interpolate.
//-----------------------------------------------------------------------
  ScatteringMomentInterpolate2Point(double Wn1, 
		      const blitz::Array<double, 2>& Pf_mom1,
		      double Wn2, 
		      const blitz::Array<double, 2>& Pf_mom2)
  : wn0(Wn1), pf0(Pf_mom1), delta_pf0((Pf_mom2 - Pf_mom1)/(Wn2 - Wn1))
  {
    range_max_check(Wn1, Wn2);
  }
//-----------------------------------------------------------------------
/// Interpolate the data to give the phase function scattering
/// moments for the given wave number. You can optionally specify the
/// number of moments and scattering matrix elements to return, the
/// default is to return all of them. This returns a matrix that is 
/// number_moment + 1 x number scattering elements in
/// size. 
//-----------------------------------------------------------------------

  blitz::Array<double, 2> 
    operator()(double wn, int nummom = -1, int numscat = -1) const
  {
    using namespace blitz;
    Range r1, r2;
    Range ra(Range::all());
    int s2;
    if(nummom == -1 || nummom - 1> pf0.rows())
      r1 = Range(0, pf0.rows() - 1);
    else
      r1 = Range(0, nummom);
    if(numscat == -1 || numscat > pf0.cols()) {
      r2 = ra;
      s2 = pf0.cols();
    } else {
      r2 = Range(0, numscat - 1);
      s2 = numscat;
    }
    int s1 = (nummom == -1 ? pf0.rows() : nummom + 1);
    Array<double, 2> res(s1, s2);
    if(nummom - 1 > pf0.rows())
      res = 0;
    res(r1, ra) = pf0(r1, r2) + delta_pf0(r1, r2) * (wn - wn0);
    return res;
  }
private:
  double wn0;
  blitz::Array<double, 2> pf0;
  blitz::Array<double, 2> delta_pf0;
};

class ScatteringMomentInterpolate {
public:
  template<class I1, class I2> ScatteringMomentInterpolate(I1 xstart, I1 xend,
							   I2 ystart)
  {
    while(xstart != xend) {
      double x0 = *xstart++;
      const blitz::Array<double, 2>& y0 = *ystart++;
      if(xstart == xend)
	return;
      double x1 = *xstart;
      const blitz::Array<double, 2>& y1 = *ystart;
      if(x1 <= x0)
	throw Exception("X needs to be sorted");
      inter[x1] = boost::shared_ptr<ScatteringMomentInterpolate2Point>
	(new ScatteringMomentInterpolate2Point(x0, y0, x1, y1));
    }
  }
  blitz::Array<double, 2> 
  operator()(double x, int nummom = -1, int nscatt = -1) const
  { 
    typedef
      std::map<double, boost::shared_ptr<ScatteringMomentInterpolate2Point> >::const_iterator Itype;
    Itype i = inter.lower_bound(x);
    if(i == inter.end())	// If we are past the upper end, 
				// then extrapolate
      return (*inter.rbegin()->second)(x, nummom, nscatt);
    else
      return (*i->second)(x, nummom, nscatt);
  }
private:
  std::map<double, boost::shared_ptr<ScatteringMomentInterpolate2Point> > inter;
};  

}
#endif
