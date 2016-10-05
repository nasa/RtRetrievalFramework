#ifndef AEROSOL_PROPERTY_H
#define AEROSOL_PROPERTY_H
#include "printable.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This gives the Aerosol properties for an Aerosol.
*******************************************************************/
class AerosolProperty : public Printable<AerosolProperty> {
public:
  virtual ~AerosolProperty() {}

//-----------------------------------------------------------------------
/// Return extinction coefficient for the given wave number.
/// \param wn - Wavenumber
//-----------------------------------------------------------------------

  virtual double extinction_coefficient(double wn) const = 0;

//-----------------------------------------------------------------------
/// Return scattering coefficient for the given wave number.
/// \param wn - Wavenumber
//-----------------------------------------------------------------------

  virtual double scattering_coefficient(double wn) const = 0;

//-----------------------------------------------------------------------
/// Return phase function moments for the given wave number.
///
/// Note that we use the de Rooij convention for the scattering matrix
/// moments. 
///
/// \param wn Wavenumber
/// \param nmom Optional number of moments to return. Default is all
///     moments.
/// \param nscatt Optional number of scattering elements to
///     return. Default is all of them.
/// \return Phase function moment. This is nmom + 1 x number
/// scattering elements. 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> phase_function_moment(double wn, 
			int nmom = -1, int nscatt = -1) const = 0;

  virtual void print(std::ostream& Os) const { Os << "AerosolProperty";}
};
}
#endif
