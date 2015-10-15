// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "merra_aerosol.h"
%}
%base_import(generic_object)
%import "hdf_file.i"
%import "double_with_unit.i"
%import "pressure.i"
%import "aerosol_optical.i"
%import "composite_initial_guess.i"
%fp_shared_ptr(FullPhysics::MerraAerosol);

namespace FullPhysics {
class MerraAerosol: public GenericObject {
public:
  MerraAerosol(const HdfFile& Merra_climatology,
	       const HdfFile& Aerosol_property,
	       DoubleWithUnit Latitude, DoubleWithUnit Longitude,
	       const boost::shared_ptr< Pressure > &Press,
	       const blitz::Array<double, 2>& Aerosol_cov,
	       double Max_aod = 0.2,
	       double Exp_aod = 0.8,
	       int Min_types = 2,
	       int Max_types = 4,
	       double Max_residual = 0.005,
	       double Reference_wn=1e4/0.755);
  %python_attribute(aerosol, boost::shared_ptr<AerosolOptical>); 
  %python_attribute(initial_guess, boost::shared_ptr<InitialGuessBuilder>);
  %python_attribute(number_merra_particle, int);
};
}
