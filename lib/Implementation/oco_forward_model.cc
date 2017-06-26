#include "oco_forward_model.h"
#include "ostream_pad.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(OcoForwardModel, ForwardModel)
.def(luabind::constructor<
       const boost::shared_ptr<Instrument>&,
       const boost::shared_ptr<SpectralWindow>&,
       const boost::shared_ptr<Level1b>&,
       const boost::shared_ptr<RadiativeTransfer>&,
       const boost::shared_ptr<SpectrumSampling>&,
       const boost::shared_ptr<StateVector>&,
       const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& >())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

OcoForwardModel::OcoForwardModel(
      const boost::shared_ptr<Instrument>& Inst,
      const boost::shared_ptr<SpectralWindow>& Spectral_window,
      const boost::shared_ptr<Level1b>& Level1b,
      const boost::shared_ptr<RadiativeTransfer>& Rt,
      const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling,
      const boost::shared_ptr<StateVector>& Sv,
      const std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > >& Spectrum_effect)
: spec_effect(Spectrum_effect), inst(Inst), swin(Spectral_window), 
  l1b(Level1b), rt(Rt), spectrum_sampling_(Spectrum_sampling), statev(Sv)
{
  if(spec_effect.size() == 0)
    spec_effect.resize(number_spectrometer());
  if(number_spectrometer() != (int) spec_effect.size())
    throw Exception("Spectrum effect needs to be the same size as the number of spectrometers");
}

// See bass class for description.
Spectrum OcoForwardModel::radiance
(int Spec_index, bool Skip_jacobian) const
{
  if(!g)
    throw Exception ("setup_grid needs to be called before calling radiance");
  range_check(Spec_index, 0, number_spectrometer());
  Spectrum highres_spec =
    rt->reflectance(g->high_resolution_grid(Spec_index), Spec_index, 
		    Skip_jacobian);
  notify_spectrum_update(highres_spec, "high_res_rt", Spec_index);

  return apply_spectrum_corrections(highres_spec, Spec_index);
}

//-----------------------------------------------------------------------
/// Applies corrections and modeling to modeled spectrum:
/// * Interpolation to uniform grid
/// * Application of spectrum effects
/// * Application of instrument model
//-----------------------------------------------------------------------

Spectrum OcoForwardModel::apply_spectrum_corrections(const Spectrum& highres_spec, int Spec_index) const 
{
  if(!g)
    throw Exception ("setup_grid needs to be called before calling apply_spectrum_corrections");
  Spectrum highres_spec_intepolated =
    g->interpolate_spectrum(highres_spec, Spec_index);
  notify_spectrum_update(highres_spec_intepolated, "high_res_interpolated", Spec_index);

  BOOST_FOREACH(const boost::shared_ptr<SpectrumEffect>& i, spec_effect[Spec_index]) {
    i->apply_effect(highres_spec_intepolated, *g);
    notify_spectrum_update(highres_spec_intepolated, "high_res_spec_effect_" + i->name(), Spec_index);
  }

  return inst->apply_instrument_model(highres_spec_intepolated, 
				      g->pixel_list(Spec_index), Spec_index);
}

// See bass class for description.
Spectrum OcoForwardModel::measured_radiance
(int Spec_index) const
{
  if(!g)
    throw Exception ("setup_grid needs to be called before calling measured_radiance");
  range_check(Spec_index, 0, number_spectrometer());
  Spectrum full(inst->pixel_spectral_domain(Spec_index),
		l1b->radiance(Spec_index));
  const std::vector<int>& plist = g->pixel_list(Spec_index);
  Array<double, 1> res_d((int) plist.size());
  ArrayAd<double, 1> res_r((int) plist.size(), 
			   full.spectral_range().data_ad().number_variable());
  Array<double, 1> uncer;
  if(full.spectral_range().uncertainty().rows() > 0)
    uncer.resize(res_r.rows());
  for(int i = 0; i < res_d.rows(); ++i) {
    res_d(i) = full.spectral_domain().data()(plist[i]);
    AutoDerivative<double> t = full.spectral_range().data_ad()(plist[i]);
    res_r(i) = t;
    if(uncer.rows() > 0)
      uncer(i) = full.spectral_range().uncertainty()(plist[i]);
  }
  return Spectrum(SpectralDomain(res_d, full.spectral_domain().units()),
		  SpectralRange(res_r, full.spectral_range().units(),
				uncer));
}

void OcoForwardModel::notify_spectrum_update(const Spectrum& updated_spec, const std::string& spec_name, int Spec_index) const {
  if (olist.size() > 0) {
    const_cast<OcoForwardModel *>(this)->notify_update_do(boost::shared_ptr<NamedSpectrum>(new NamedSpectrum(updated_spec, spec_name, Spec_index)));
  }  
}

void OcoForwardModel::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "OcoForwardModel:\n";
  Os << "  Input file list:\n";
  opad << input_file_description() << "\n";
  opad.strict_sync();
  Os << "  Spectral Window:\n";
  opad << *swin << "\n";
  opad.strict_sync();
  Os << "  Spectrum Sampling:\n";
  opad << *spectrum_sampling_ << "\n";
  opad.strict_sync();
  Os << "  Level 1b:\n";
  opad << *l1b << "\n";
  opad.strict_sync();
  Os << "  Instrument:\n";
  opad << *inst << "\n";
  opad.strict_sync();
  Os << "  Radiative Transfer:\n";
  opad << *rt << "\n";
  opad.strict_sync();
  for(int i = 0; i < number_spectrometer(); ++i) {
    Os << "  Spectrum Effect[" << i << "]:\n";
    BOOST_FOREACH(boost::shared_ptr<SpectrumEffect> se, spec_effect[i])
      opad << *se << "\n";
    opad.strict_sync();
  }
  if (statev->state().rows() > 0) {
    Os << "  State Vector:\n";
    opad << *statev << "\n";
    opad.strict_sync();
  }
}
