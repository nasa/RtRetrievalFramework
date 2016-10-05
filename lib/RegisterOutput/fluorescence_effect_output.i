// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "fluorescence_effect_output.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}
%base_import(register_output_base)
%import "output.i"
%import "fluorescence_effect.i"

%fp_shared_ptr(FullPhysics::FluorescenceEffectOutput);

namespace FullPhysics {
class FluorescenceEffectOutput : public RegisterOutputBase {
public:
  FluorescenceEffectOutput(const boost::shared_ptr<FluorescenceEffect>& Fluor);
  virtual ~FluorescenceEffectOutput();
  virtual void register_output(const boost::shared_ptr<Output>& out) const;
  virtual void register_output_apriori(const boost::shared_ptr<Output>& out) const;
};
}
