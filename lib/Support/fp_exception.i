// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "fp_exception.h"
%}

%base_import(generic_object)

%init {
  FullPhysics::no_gsl_abort();
}

%fp_shared_ptr(FullPhysics::Exception)
namespace FullPhysics {
  class Exception : public GenericObject {
  public:
    Exception(const std::string& W);
    std::string print_to_string() const;
    const char* what();
  };

void no_gsl_abort();
}


