// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%{
#include "absorber_vmr_scaled.h"
%}
%base_import(absorber_vmr_imp_base)

%fp_shared_ptr(FullPhysics::AbsorberVmrScaled)
namespace FullPhysics {
class AbsorberVmrScaled : public AbsorberVmrImpBase {
public:
  AbsorberVmrScaled(const boost::shared_ptr<Pressure>& Press,
		    double Scale,                         
		    bool Scale_flag,
		    const std::string& Gas_name);
  virtual boost::shared_ptr<AbsorberVmr> clone();
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(scale_factor, double)
  %python_attribute(scale_uncertainty, double)
protected:
  virtual void calc_vmr() const;
};
}

