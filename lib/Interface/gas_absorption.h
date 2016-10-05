#ifndef GAS_ABSORPTION_H
#define GAS_ABSORPTION_H
#include "printable.h"
#include "double_with_unit.h"
#include "auto_derivative_with_unit.h"
#include "array_ad.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class determine the gaseous absorption coefficient for a given
  wave number, temperature and pressure.
*******************************************************************/
class GasAbsorption: public Printable<GasAbsorption> {
public:
  virtual ~GasAbsorption() {}

//-----------------------------------------------------------------------
/// Return true if we have data for the given wave number. A
/// particular gas might not have absorption coefficients for all
/// spectral bands, e.g., ABSCO tables.
//-----------------------------------------------------------------------

  virtual bool have_data(double wn) const = 0;

//-----------------------------------------------------------------------
/// For some tables, we might have a broadener (e.g., "h2o"). This
/// returns the name of the broadener, if any.
//-----------------------------------------------------------------------
  
  virtual std::string broadener_name() const = 0;

//-----------------------------------------------------------------------
/// This interpolates the ABSCO data to give absorption cross section
/// for a given pressure, temperature, and broadener VMR.
///
/// \param Wn wave number
/// \param Press Pressure
/// \param Temp Temperature
/// \param Broadener_vmr Broadner VMR (e.g., H2O VMR). Not all
///   tables will make use of this information.
/// \return Absorption cross section in cm^2 / molecule
//-----------------------------------------------------------------------

  virtual DoubleWithUnit absorption_cross_section(double Wn, 
     const DoubleWithUnit& Press, 
     const DoubleWithUnit& Temp,
     const DoubleWithUnit& Broadener_vmr) const = 0;

//-----------------------------------------------------------------------
/// This interpolates the ABSCO data to give absorption cross section
/// for a given pressure, temperature, and broadener VMR.
///
/// \param Wn wave number
/// \param Press Pressure
/// \param Temp Temperature
/// \param Broadener_vmr Broadner VMR (e.g., H2O VMR). Not all
///   tables will make use of this information.
/// \return Absorption cross section in cm^2 / molecule
//-----------------------------------------------------------------------

  virtual AutoDerivativeWithUnit<double>
  absorption_cross_section(double Wn, 
    const DoubleWithUnit& Press, 
    const AutoDerivativeWithUnit<double>& Temp,
    const AutoDerivativeWithUnit<double>& Broadener_vmr) const = 0;

  virtual void print(std::ostream& Os) const {Os << "GasAbsorption";}
};
}
#endif
