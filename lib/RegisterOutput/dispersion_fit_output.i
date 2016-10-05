// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "dispersion_fit_output.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "dispersion_fit.i"

%fp_shared_ptr(FullPhysics::DispersionFitOutput);

namespace FullPhysics {
class DispersionFitOutput : public RegisterOutputBase {
public:
  DispersionFitOutput(const boost::shared_ptr<DispersionFit>& D);
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}


