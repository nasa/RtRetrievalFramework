// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "forward_model_spectral_grid.h"
#include "instrument.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
%}

%base_import(generic_object)
%import "spectral_domain.i"
%import "spectrum.i"
%import "instrument.i"
%import "spectral_window.i"
%import "spectrum_sampling.i"

%fp_shared_ptr(FullPhysics::ForwardModelSpectralGrid);

namespace FullPhysics {
class ForwardModelSpectralGrid  : public GenericObject {
public:
  ForwardModelSpectralGrid(
   const boost::shared_ptr<Instrument>& Inst,
   const boost::shared_ptr<SpectralWindow>& Spectral_window,
   const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling);
  std::string print_to_string() const;
  %python_attribute(number_spectrometer, virtual int);
  const SpectralDomain low_resolution_grid(int Spec_index) const;
  const SpectralDomain high_resolution_grid(int Spec_index) const;
  const SpectralDomain high_resolution_interpolated_grid(int Spec_index) const;
  Spectrum interpolate_spectrum(const Spectrum& Spec_in, int Spec_index) const;
  const std::vector<int> pixel_list(int Spec_index) const;
};
}
