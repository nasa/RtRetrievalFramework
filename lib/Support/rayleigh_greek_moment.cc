#include "rayleigh_greek_moment.h"
using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Return the Ralyeigh Greek Moment Array. This is in Fortran order,
/// so you can directly pass it to Fortran code.
//-----------------------------------------------------------------------

const blitz::Array<double, 2>& RayleighGreekMoment::array(double depolar_fact)
{
  static bool first_time = true;
  static Array<double, 2> res;

  if(first_time) {
    res.reference(Array<double, 2>(shape(3,6), ColumnMajorArray<2>()));
    res = 0;
    res(0,0) = 1;
    res(1,0) = 1e-11; // this is to ensure LIDORT doesn't freak out
    res(1,3) = 3*(1-2*depolar_fact)/(2+depolar_fact);
    res(2,0) = (1-depolar_fact)/(2+depolar_fact);
    res(2,4) = sqrt(6)*(1-depolar_fact)/(2+depolar_fact);
    res(2,1) = 6*(1-depolar_fact)/(2+depolar_fact);
    first_time = false;
  }
  return res;
}
