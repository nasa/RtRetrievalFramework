// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "absorber_vmr_level_scaled.h"
%}
%base_import(absorber_vmr_scaled)
%import "pressure.i"
%fp_shared_ptr(FullPhysics::AbsorberVmrLevelScaled)
namespace FullPhysics {
class AbsorberVmrLevelScaled : public AbsorberVmrScaled {
public:
  AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
			 const blitz::Array<double, 1>& Vmr_profile,
			 double Scale,                         
			 bool Scale_flag,
			 const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;
  %python_attribute(scale_factor, double)
  %python_attribute(scale_uncertainty, double)
};
}

