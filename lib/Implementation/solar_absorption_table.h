#ifndef SOLAR_ABSORPTION_TABLE_H
#define SOLAR_ABSORPTION_TABLE_H
#include "solar_absorption_spectrum.h"
#include "hdf_file.h"
#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the solar absorption spectrum.

  This particular implementation reads a table of values, and
  interpolates to give a value.
*******************************************************************/

class SolarAbsorptionTable : public SolarAbsorptionSpectrum {  
public:
  SolarAbsorptionTable(const HdfFile& Hdf_static_input,
			 const std::string& Hdf_group);
  virtual ~SolarAbsorptionTable() {}
  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_absorption_spectrum(
     const SpectralDomain& spec_domain) const;
private:
  Unit domain_unit;
  LinearInterpolate<double, double> table;
  std::string hdf_file_name;
  std::string hdf_group;
};
}
#endif
