// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "aerosol_property_hdf.h"
%}
%base_import(aerosol_property)
%import "hdf_file.i"
%fp_shared_ptr(FullPhysics::AerosolPropertyHdf)
namespace FullPhysics {
class AerosolPropertyHdf : public AerosolProperty {
public:
  AerosolPropertyHdf(const HdfFile& F, const std::string& Group_name);
  virtual double extinction_coefficient(double wn) const;
  virtual double scattering_coefficient(double wn) const;
  virtual blitz::Array<double, 2> phase_function_moment(double wn, 
		int nmom = -1, int nscatt = -1) const;

};
}

