// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "fluorescence_effect.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
%}
%base_import(spectrum_effect_imp_base)
%import "rt_atmosphere.i"
%import "spectrum_sampling.i"
%import "double_with_unit.i"
%import "stokes_coefficient.i"
%fp_shared_ptr(FullPhysics::FluorescenceEffect);

namespace FullPhysics {
class FluorescenceEffect : public SpectrumEffectImpBase {
public:
  FluorescenceEffect(const blitz::Array<double, 1>& Coeff,
                     const blitz::Array<bool, 1>& Used_flag,
                     const boost::shared_ptr<RtAtmosphere>& Atm,
		     const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
                     const DoubleWithUnit& Sza, 
                     const int Spec_index,
                     const DoubleWithUnit& Reference,
                     const Unit& Retrieval_unit);
  virtual void apply_effect(Spectrum& Spec,
		    const ForwardModelSpectralGrid& Forward_model_grid) const;
  virtual boost::shared_ptr<SpectrumEffect> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  %python_attribute(name, std::string)
  %python_attribute(fluorescence_at_reference, double)
  %python_attribute(fluorescence_at_reference_uncertainty, double)
  %python_attribute(fluorescence_slope, double)
  %python_attribute(fluorescence_slope_uncertainty, double)
 };
}
