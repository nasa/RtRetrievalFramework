#include "ils_instrument.h"
#include "ostream_pad.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IlsInstrument, Instrument)
.def(luabind::constructor<const std::vector<boost::shared_ptr<Ils> >&>())
.def(luabind::constructor<const std::vector<boost::shared_ptr<Ils> >&,
			  const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >& >())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

IlsInstrument::IlsInstrument(const std::vector<boost::shared_ptr<Ils> >& Ils_list,
			     const std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > >&
			     Instrument_correction)
: ils_(Ils_list), inst_corr(Instrument_correction)
{
  if(inst_corr.size() == 0)
    inst_corr.resize(ils_.size());
  if(ils_.size() != inst_corr.size())
    throw Exception("Ils and Instrument_correction need to be the same size");
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    i->add_observer(*this);
  BOOST_FOREACH(std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
		inst_corr) {
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection>& j, i)
      j->add_observer(*this);
  }
}

boost::shared_ptr<Instrument> IlsInstrument::clone() const
{
  std::vector<boost::shared_ptr<Ils> > ils_vec;
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    ils_vec.push_back(i->clone());
  std::vector<std::vector<boost::shared_ptr<InstrumentCorrection> > > 
    inst_corr_vec;
  BOOST_FOREACH(const std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
		inst_corr) {
    std::vector<boost::shared_ptr<InstrumentCorrection> > t;
    BOOST_FOREACH(const boost::shared_ptr<InstrumentCorrection>& j, i)
      t.push_back(j->clone());
    inst_corr_vec.push_back(t);
  }
  
  return boost::shared_ptr<Instrument>(new IlsInstrument(ils_vec, 
							 inst_corr_vec));
}

// See base class for description.
void IlsInstrument::notify_add(StateVector& Sv) 
{ 
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    Sv.add_observer(*i); 
  BOOST_FOREACH(std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
		inst_corr) {
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection>& j, i)
      Sv.add_observer(*j);
  }
}
// See base class for description.

void IlsInstrument::notify_remove(StateVector& Sv) 
{ 
  BOOST_FOREACH(boost::shared_ptr<Ils> i, ils_)
    Sv.remove_observer(*i);
  BOOST_FOREACH(std::vector<boost::shared_ptr<InstrumentCorrection> >& i, 
		inst_corr) {
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection>& j, i)
      Sv.remove_observer(*j);
  }
}

Spectrum IlsInstrument::apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const 
{
  range_check(Spec_index, 0, number_spectrometer());

  SpectralDomain full = pixel_spectral_domain(Spec_index);
  blitz::Array<double, 1> res_sd((int) Pixel_list.size());
  for(int i = 0; i < res_sd.rows(); ++i)
    res_sd(i) = full.data()(Pixel_list[i]);
  SpectralDomain res_dom(res_sd, full.units());

  // Gain some speed advantages when running things without autoderivatives
  // if they are not needed.
  SpectralRange res_sr;
  if(High_resolution_spectrum.spectral_range().data_ad().
       number_variable() > 0) {
    ArrayAd<double, 1> rad_ad =
      ils_[Spec_index]->apply_ils
      (High_resolution_spectrum.spectral_domain().data(), 
       High_resolution_spectrum.spectral_range().data_ad(),
       Pixel_list);
    res_sr = SpectralRange(rad_ad, 
			   High_resolution_spectrum.spectral_range().units());
  } else {
    Array<double, 1> rad =
      ils_[Spec_index]->apply_ils
      (High_resolution_spectrum.spectral_domain().data(), 
       High_resolution_spectrum.spectral_range().data(),
       Pixel_list);
    res_sr = SpectralRange(rad, 
			   High_resolution_spectrum.spectral_range().units());
  }
  BOOST_FOREACH(const boost::shared_ptr<InstrumentCorrection>& i, 
		inst_corr[Spec_index])
    i->apply_correction(ils_[Spec_index]->pixel_grid(), Pixel_list, res_sr);

  return Spectrum(res_dom, res_sr);
}

void IlsInstrument::print(std::ostream& Os) const 
{ 
  Os << "IlsInstrument:\n";
  OstreamPad opad(Os, "    ");
  for(int i = 0; i < number_spectrometer(); ++i) {
    Os << "  " << band_name(i) << ":\n";
    opad << *ils_[i] << "\n";
    BOOST_FOREACH(boost::shared_ptr<InstrumentCorrection> ic, inst_corr[i])
      opad << *ic << "\n";
    opad.strict_sync();
  }
}
