#ifndef RAYLEIGH_GREEK_MOMENT_H
#define RAYLEIGH_GREEK_MOMENT_H
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class provides the Rayleigh Greek Moments. This is really just
  a constant, but we put this in a central class to have the code for
  this in one place.
*******************************************************************/
class RayleighGreekMoment {
public:
  static const blitz::Array<double, 2>& array
    (double depolar_fact = 0.02790      // depolarization factor; for air is 0.02790 (young, 1980)
     );
};
}
#endif
