#ifndef AEROSOL_PROPERTY_RH_HDF_H
#define AEROSOL_PROPERTY_RH_HDF_H
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

  This variation of AerosolProperty interpolates the aerosol
  properties by the relative humidity.
*******************************************************************/
class AerosolPropertyRhHdf : public AerosolPropertyImpBase {
public:
  AerosolPropertyRhHdf(const HdfFile& F, const std::string& Group_name,
		       const boost::shared_ptr<Pressure>& Press,
		       const boost::shared_ptr<RelativeHumidity>& Rh);
  virtual ~AerosolPropertyRhHdf() {}
  virtual boost::shared_ptr<AerosolProperty> clone() const;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) 
    const;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer_not_used(double wn) 
    const;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn)
    const;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer_not_used(double wn)
    const;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const;
  ArrayAd<double, 3> 
  phase_function_moment_each_layer_not_used(double wn, int nmom = -1, 
					    int nscatt = -1) const;
  virtual void print(std::ostream& Os) const;
private:
  boost::shared_ptr<RelativeHumidity> rh;
  std::vector<AutoDerivative<double> > rh_val;
  std::vector<double> rh_val_d; // Temporary
  std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > qext;
  std::vector<boost::shared_ptr<LinearInterpolate<double, double> > > qscat;
  std::vector<boost::shared_ptr<ScatteringMomentInterpolate> > pf;
  std::string hdf_file, hdf_group;
};
}
#endif
