// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "absorber_vmr_shape.h"
%}
%base_import(absorber_vmr_imp_base)
%import "pressure.i"

%fp_shared_ptr(FullPhysics::AbsorberVmrShape)
namespace FullPhysics {
class AbsorberVmrShape : public AbsorberVmrImpBase {
public:
    AbsorberVmrShape(const blitz::Array<double, 1> VMR_base,
                     const blitz::Array<double, 2> Shape_profile,
                     const blitz::Array<double, 1> Shape_scaling,
                     const boost::shared_ptr<Pressure>& Press,
                     const blitz::Array<bool, 1>& Shape_flag,
                     const std::string& Gas_name,
                     const bool Log_profiles = false);
    virtual boost::shared_ptr<AbsorberVmr> clone() const;
    virtual boost::shared_ptr<AbsorberVmr> clone(const boost::shared_ptr<Pressure>& Press) const;
    %python_attribute(sub_state_identifier, std::string);
    virtual std::string state_vector_name_i(int i) const;
    %python_attribute(vmr_base, blitz::Array<double, 1>);
    %python_attribute(shape_profile, blitz::Array<double, 2>);
    %python_attribute(shape_scaling, blitz::Array<double, 1>);
    %python_attribute(pressure_profile, blitz::Array<double, 1>);
};
}
