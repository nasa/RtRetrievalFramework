// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "sub_state_vector_array.h"
#include "merra_aerosol.h"
#include "absorber.h"
#include "temperature.h"
#include "altitude.h"
%}
%base_import(generic_object)
%import "hdf_file.i"
%import "double_with_unit.i"
%import "pressure.i"
%import "aerosol.i"
%import "composite_initial_guess.i"
%fp_shared_ptr(FullPhysics::MerraAerosol);

namespace FullPhysics {
class MerraAerosol: public GenericObject {
public:
  MerraAerosol(const HdfFile& Merra_climatology,
	       const HdfFile& Aerosol_property,
	       DoubleWithUnit Latitude, DoubleWithUnit Longitude,
	       const boost::shared_ptr<Pressure> &Press,
	       const boost::shared_ptr<RelativeHumidity> &Rh,
	       const blitz::Array<double, 2>& Aerosol_cov,
	       double Max_aod = 0.2,
	       double Exp_aod = 0.8,
	       int Min_types = 2,
	       int Max_types = 4,
	       bool Linear_aod = false,
	       bool Relative_humidity_aerosol = false,
	       double Max_residual = 0.005,
	       double Reference_wn=1e4/0.755);
  %python_attribute(aerosol, boost::shared_ptr<Aerosol>); 
  %python_attribute(initial_guess, boost::shared_ptr<InitialGuessBuilder>);
  %python_attribute(number_merra_particle, int);
};
}
