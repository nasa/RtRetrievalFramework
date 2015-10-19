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
  virtual boost::shared_ptr<AerosolProperty> clone() const;
  virtual boost::shared_ptr<AerosolProperty> 
  clone(const boost::shared_ptr<Pressure>& Press) const;
  virtual ArrayAd<double, 1> extinction_coefficient_each_layer(double wn) const;
  virtual ArrayAd<double, 1> scattering_coefficient_each_layer(double wn) const;
  virtual ArrayAd<double, 3> 
  phase_function_moment_each_layer(double wn, int nmom = -1, 
				   int nscatt = -1) const;

};
}

