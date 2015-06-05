// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "oco_forward_model.h"
#include "ils_instrument.h"
#include "sub_state_vector_array.h"
#include "pressure.h"
%}
%base_import(forward_model)
%base_import(named_spectrum)
%import "instrument.i"
%import "spectral_window.i"
%import "level_1b.i"
%import "radiative_transfer.i"
%import "spectrum_sampling.i"
%import "state_vector.i"
%import "spectrum_effect.i"
%import "forward_model_spectral_grid.i"

%fp_shared_ptr(FullPhysics::OcoForwardModel);

namespace FullPhysics {
class OcoForwardModel : public ForwardModel,
   public Observable<boost::shared_ptr<FullPhysics::NamedSpectrum> > {
public:
  OcoForwardModel(
   const boost::shared_ptr<Instrument>& Inst,
   const boost::shared_ptr<SpectralWindow>& Spectral_window,
   const boost::shared_ptr<Level1b>& Level_1b,
   const boost::shared_ptr<RadiativeTransfer>& Rt,
   const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling,
   const boost::shared_ptr<StateVector>& Sv,
   const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& Spectrum_effect = 
		  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >());
  virtual ~OcoForwardModel();
  %python_attribute(state_vector, boost::shared_ptr<StateVector>)
  virtual Spectrum radiance(int Spec_index, bool Skip_jacobian = false) 
    const;
  virtual Spectrum measured_radiance(int Spec_index) const;
  %python_attribute_with_set(instrument, boost::shared_ptr<Instrument>)
  %python_attribute_with_set(spectral_window, boost::shared_ptr<SpectralWindow>)
  %python_attribute_with_set(level_1b, boost::shared_ptr<Level1b>)
  %python_attribute_with_set(radiative_transfer, boost::shared_ptr<RadiativeTransfer>)
  %python_attribute(spectrum_sampling, boost::shared_ptr<SpectrumSampling>)
  %python_attribute(spectral_grid, boost::shared_ptr<ForwardModelSpectralGrid>)
  Spectrum apply_spectrum_corrections(const Spectrum& highres_spec, int Spec_index) const;

  virtual void add_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs); 
  virtual void remove_observer(Observer<boost::shared_ptr<FullPhysics::NamedSpectrum> >& Obs);

  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int Spec_index) const;

  // vector of vector for SpectrumEffect is kind of a pain in
  // python. So just brute force a conversion.
  %extend {
    int _speceff_size() { return (int) $self->spectrum_effect().size(); }
    int speceff_size2(int i) { return (int) $self->spectrum_effect()[i].size();}
    boost::shared_ptr<SpectrumEffect> speceff_val(int i, int j)
    {return $self->spectrum_effect()[i][j];}
  }
  %pythoncode {
@property
def speceff_size(self):
  return self._speceff_size()

@property
def spectrum_effect(self):
  res = []
  for i in range(self.speceff_size):
     res2 = []
     for j in range(self.speceff_size2(i)):
        res2.append(self.speceff_val(i,j))
     res.append(res2)
  return res
  }
};
}
