// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "perturbation.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::Perturbation);

namespace FullPhysics {
class Perturbation : public GenericObject {
public:
  virtual ~Perturbation();
  std::string print_to_string() const;
  %python_attribute_abstract(perturbation, blitz::Array<double, 1>);
};
}
