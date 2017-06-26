// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "noise_model.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::NoiseModel);

namespace FullPhysics {
class NoiseModel : public GenericObject {
public:
  virtual ~NoiseModel() {}
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const = 0;
  std::string print_to_string() const;
};
}
