#ifndef AEROSOL_PROPERTY_HDF_H
#define AEROSOL_PROPERTY_HDF_H
#include "aerosol_property_imp_base.h"
#include "hdf_file.h"
#include "linear_interpolate.h"
#include "scattering_moment_interpolator.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  This gives the Aerosol properties for an Aerosol. This particular
  implementation reads the Aerosol properties from the HDF group in
  a HDF file. The fields "wave_number", "extinction_coefficient",
  "scattering_coefficient" and "phase_function_moment" are read.

  The HDF file supplies the particle properties for a few
  wavenumbers. We then linearly interpolate to get the aerosol
  properties for wavenumbers between these value.  If a wavenumber
  outside the range of the file is requested, then we extrapolate to
  get the value.
*******************************************************************/
class AerosolPropertyHdf : public AerosolPropertyImpBase {
public:
  AerosolPropertyHdf(const HdfFile& F, const std::string& Group_name);
  virtual ~AerosolPropertyHdf() {}
  virtual boost::shared_ptr<AerosolProperty> clone() const;
  virtual double extinction_coefficient(double wn) const 
  { return (*qext)(wn); }
  virtual double scattering_coefficient(double wn) const
  { return (*qscat)(wn); }
  virtual blitz::Array<double, 2> phase_function_moment(double wn, 
	int nmom = -1, int nscatt = -1) const
  { return (*pf)(wn, nmom, nscatt); }
  virtual void print(std::ostream& Os) const;
private:
  boost::shared_ptr<LinearInterpolate<double, double> > qext;
  boost::shared_ptr<LinearInterpolate<double, double> > qscat;
  boost::shared_ptr<ScatteringMomentInterpolate> pf;
  std::string hdf_file, hdf_group;
};
}
#endif
