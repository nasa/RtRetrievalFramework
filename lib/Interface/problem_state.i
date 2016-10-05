// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "problem_state.h"
%}

%base_import(generic_object)
%fp_shared_ptr(FullPhysics::ProblemState);

namespace FullPhysics {
class ProblemState : public GenericObject {
public:
  ProblemState();
  ProblemState(const ProblemState& s);
  virtual ~ProblemState();
  virtual void set(const ProblemState& s);
  virtual void clear();

  virtual bool parameters_different(const blitz::Array<double, 1>& x) const;
  %python_attribute_with_set(parameters, blitz::Array<double, 1>);
  %python_attribute(parameter_size, int);
  %python_attribute_abstract(expected_parameter_size, int);
  virtual void assert_parameter_set_correctly() const;
  virtual void assert_parameter_correct(const blitz::Array<double, 1>& x) const;
  std::string print_to_string() const;
};
}
