// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "initial_guess.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::InitialGuess);

namespace FullPhysics {
class InitialGuess : public GenericObject {
public:
  virtual ~InitialGuess();
  std::string print_to_string() const;
  %python_attribute_abstract(initial_guess, blitz::Array<double, 1>)
  %python_attribute(apriori, virtual blitz::Array<double, 1>)
  %python_attribute(apriori_covariance, virtual blitz::Array<double, 2>)
};
}
