// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "level_1b_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "level_1b.i"

%fp_shared_ptr(FullPhysics::Level1bOutput);

namespace FullPhysics {
class Level1bOutput : public RegisterOutputBase {
public:
  Level1bOutput(const boost::shared_ptr<Level1b>& F);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


