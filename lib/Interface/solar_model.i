// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include <std_vector.i>
%include "common.i"

%{
#include "solar_model.h"
#include "pressure.h"
#include "sub_state_vector_array.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}

%base_import(spectrum_effect)

%fp_shared_ptr(FullPhysics::SolarModel);

namespace FullPhysics {
class SolarModel : public SpectrumEffect {
public:
  virtual ~SolarModel();
  virtual Spectrum apply_solar_model(const Spectrum& Spec) const;
  virtual Spectrum solar_spectrum(const SpectralDomain& Spec_domain) const = 0;
  virtual void apply_effect(Spectrum& Spec, const ForwardModelSpectralGrid& Forward_model_grid) const;
};
}
%template(vector_solar_model) std::vector<boost::shared_ptr<FullPhysics::SolarModel> >;
