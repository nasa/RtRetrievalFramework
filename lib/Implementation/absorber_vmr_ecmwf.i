// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "absorber_vmr_ecmwf.h"
%}
%base_import(absorber_vmr_scaled)
%import "ecmwf.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AbsorberVmrEcmwf)
namespace FullPhysics {
class AbsorberVmrEcmwf : public AbsorberVmrScaled {
public:
  AbsorberVmrEcmwf(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Scale,                         
		   bool Scale_flag,
		   const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;
  %python_attribute(scale_factor, double)
  %python_attribute(scale_uncertainty, double)
};
}

