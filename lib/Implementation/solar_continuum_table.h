#ifndef SOLAR_CONTINUUM_TABLE_H
#define SOLAR_CONTINUUN_TABLE_H
#include "hdf_file.h"
#include "solar_continuum_spectrum.h"
#include "linear_interpolate.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This class calculates the solar continuum spectrum.

  This particular implementation uses a table to calculate the Solar
  Planck Function, doing a linear interpolation for points in between.
*******************************************************************/

class SolarContinuumTable : public SolarContinuumSpectrum {
public:
  SolarContinuumTable(const HdfFile& F,
		      const std::string& Hdf_group,
		      bool Convert_from_photon = true);
  virtual ~SolarContinuumTable() {}
  virtual void print(std::ostream& Os) const;
  virtual Spectrum solar_continuum_spectrum(
     const SpectralDomain& spec_domain) const;
private:
  Unit domain_unit, range_unit;
  LinearInterpolate<double, double> table;
  bool convert_from_photon;
  std::string hdf_file_name;
  std::string hdf_group;
};
}
#endif
