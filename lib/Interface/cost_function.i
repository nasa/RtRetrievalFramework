// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "cost_function.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::CostFunction);

namespace FullPhysics {
class CostFunction : public GenericObject {
public:
  virtual ~CostFunction();
  std::string print_to_string() const;
  virtual void cost_function(const blitz::Array<double, 1>& X,
			blitz::Array<double, 1>& OUTPUT,
			blitz::Array<double, 1>& OUTPUT,
			blitz::Array<double, 2>& OUTPUT) const = 0;
};
}
