// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include <std_vector.i>
%include "common.i"
%{
#include "register_output_base.h"
%}
%import "output.i"
%fp_shared_ptr(FullPhysics::RegisterOutputBase)

namespace FullPhysics {
// Allow these classes to be derived from in Python.
%feature("director") RegisterOutputBase;

class RegisterOutputBase {
public:
  virtual ~RegisterOutputBase();
  std::string print_to_string() const;
  virtual std::string desc() const;
  virtual void register_output(const boost::shared_ptr<Output>& out) const = 0;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}
%template(vector_register_output) std::vector<boost::shared_ptr<FullPhysics::RegisterOutputBase> >;

