#ifndef OCO_FORWARD_MODEL_H
#define OCO_FORWARD_MODEL_H
#include "forward_model.h"
#include "forward_model_spectral_grid.h"
#include "level_1b.h"
#include "radiative_transfer.h"
#include "spectrum_effect.h"
#include "named_spectrum.h"
#include "spectral_window.h"
#include "instrument.h"
#include "spectrum_sampling.h"
#include "observer.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This is the forward model used form GOSAT/OCO. This is fairly
  general, we may want to rename this at some point if we use it for
  other instruments.
*******************************************************************/

class OcoForwardModel : public ForwardModel,
  public Observable<boost::shared_ptr<NamedSpectrum> > {
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
  virtual ~OcoForwardModel() {}
  virtual void setup_grid()
  {
    g.reset(new ForwardModelSpectralGrid(inst, swin, spectrum_sampling_));
  }
  virtual boost::shared_ptr<StateVector> state_vector() const { return statev; }
  virtual int number_spectrometer() const {return swin->number_spectrometer();}
  virtual std::string hdf_band_name(int Spec_index) const 
  { return inst->hdf_band_name(Spec_index); }
  virtual SpectralDomain spectral_domain(int Spec_index) const
  { 
    if(!g)
      throw Exception ("setup_grid needs to be called before calling spectral_domain");
    return g->low_resolution_grid(Spec_index);
  }
  virtual SpectralDomain::TypePreference 
  spectral_domain_type_preference() const
  {
    return inst->pixel_spectral_domain(0).type_preference();
  }
  virtual Spectrum radiance(int Spec_index, bool Skip_jacobian = false) 
    const;
  virtual Spectrum measured_radiance(int Spec_index) const;
  virtual void print(std::ostream& Os) const;

  // These are useful for python, so make available here.
  const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& 
  spectrum_effect() const {return spec_effect; }
  const boost::shared_ptr<Instrument>& instrument() const
  { return inst; }
  void instrument(const boost::shared_ptr<Instrument>& V)
  { inst = V;}
  const boost::shared_ptr<SpectralWindow>& spectral_window() const
  { return swin; }
  void spectral_window(const boost::shared_ptr<SpectralWindow>& V)
  { swin = V; }
  const boost::shared_ptr<Level1b>& level_1b() const 
  { return l1b; }
  void level_1b(const boost::shared_ptr<Level1b>& V)
  { l1b = V; }
  const boost::shared_ptr<RadiativeTransfer>& radiative_transfer() const
  { return rt; }
  void radiative_transfer(const boost::shared_ptr<RadiativeTransfer>& V)
  { rt = V; }
  const boost::shared_ptr<SpectrumSampling>& spectrum_sampling()
    const { return spectrum_sampling_;}
  void spectrum_sampling(const boost::shared_ptr<SpectrumSampling>& V)
  { spectrum_sampling_ = V; }
  
  Spectrum apply_spectrum_corrections(const Spectrum& highres_spec, int Spec_index) const;
  const boost::shared_ptr<ForwardModelSpectralGrid>& spectral_grid() const
  { return g; }
  /// Required observable functions
  virtual void add_observer(Observer<boost::shared_ptr<NamedSpectrum> > & Obs) 
  { add_observer_do(Obs); }
  virtual void remove_observer(Observer<boost::shared_ptr<NamedSpectrum> >& Obs) 
  { remove_observer_do(Obs); }

  void notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int Spec_index) const;

private:
  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > spec_effect;
  boost::shared_ptr<Instrument> inst;
  boost::shared_ptr<SpectralWindow> swin;
  boost::shared_ptr<Level1b> l1b;
  boost::shared_ptr<RadiativeTransfer> rt;
  boost::shared_ptr<SpectrumSampling> spectrum_sampling_;
  boost::shared_ptr<StateVector> statev;
  boost::shared_ptr<ForwardModelSpectralGrid> g;
};
}
#endif
