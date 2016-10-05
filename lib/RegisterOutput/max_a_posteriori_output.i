// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "max_a_posteriori_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "max_a_posteriori.i"

%fp_shared_ptr(FullPhysics::MaxAPosterioriOutput);

namespace FullPhysics {
class MaxAPosterioriOutput : public RegisterOutputBase {
public:
  MaxAPosterioriOutput(const boost::shared_ptr<MaxAPosteriori>& MAP, 
		       bool Write_jacobian = false);

  virtual ~MaxAPosterioriOutput();
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
};
}


