#include "forward_model.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
typedef
void (ForwardModel::*ifd)(const std::string&);

#include "register_lua.h"
REGISTER_LUA_CLASS(ForwardModel)
.def("input_file_description", ((ifd) &ForwardModel::input_file_description))
.def("setup_grid", &ForwardModel::setup_grid)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// The range of radiance() that corresponds to a particular
/// band. Note a particular issue, the range may well be empty if a
/// band is not used at all. This is a useful edge case, but
/// unfortunately blitz::Range does not support empty ranges. As a
/// simple work around, we use the boost::optional class to return a
/// value only if the range is not empty.
//-----------------------------------------------------------------------

boost::optional<Range> ForwardModel::pixel_range(int Spec_index) const
{
  range_check(Spec_index, 0, number_spectrometer());
  int sind = 0;
  for(int i = 0; i < Spec_index; ++i)
    sind += spectral_domain(i).data().rows();
  int nrow = spectral_domain(Spec_index).data().rows();
  if(nrow > 0)
    return boost::optional<Range>(Range(sind, sind + nrow - 1));
  else
    return boost::optional<Range>();
}

//-----------------------------------------------------------------------
/// This is radiance, all stacked together as one long spectrum (so
/// band 0, followed by band 1, etc.).
//-----------------------------------------------------------------------

Spectrum ForwardModel::radiance_all(bool Skip_jacobian) const
{
  std::vector<Spectrum> sall;
  std::vector<Range> prall;
  for(int i = 0; i < number_spectrometer(); ++i) {
    boost::optional<Range> pr = pixel_range(i);
    if(pr) {
      sall.push_back(radiance(i, Skip_jacobian));
      prall.push_back(*pr);
    }
  }
  int nrow = 0;
  int nvar = 0;
  BOOST_FOREACH(const Spectrum& s, sall) {
    nrow += s.spectral_domain().data().rows();
    nvar = std::max(nvar, s.spectral_range().data_ad().number_variable());
  }
  Unit ud, ur;
  if(sall.size() > 0) {
    ud = sall[0].spectral_domain().units();
    ur = sall[0].spectral_range().units();
  }
  ArrayAd<double, 1> sr(nrow, nvar);
  Array<double, 1> sd(nrow);
  for(int i = 0; i < (int) sall.size(); ++i) {
    sd(prall[i]) = sall[i].spectral_domain().data() * 
      FullPhysics::conversion(sall[i].spectral_domain().units(), ud);
    sr.value()(prall[i]) = sall[i].spectral_range().data() *
      FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    if(nvar > 0)
      sr.jacobian()(prall[i], Range::all()) =
	sall[i].spectral_range().data_ad().jacobian() *
	FullPhysics::conversion(sall[i].spectral_range().units(), ur);
  }
  return Spectrum(SpectralDomain(sd, ud), SpectralRange(sr, ur));
}

//-----------------------------------------------------------------------
/// This is the measured radiance, all stacked together as one long 
/// spectrum (so band 0, followed by band 1, etc.).
//-----------------------------------------------------------------------

Spectrum ForwardModel::measured_radiance_all() const
{
  std::vector<Spectrum> sall;
  std::vector<Range> prall;
  for(int i = 0; i < number_spectrometer(); ++i) {
    boost::optional<Range> pr = pixel_range(i);
    if(pr) {
      sall.push_back(measured_radiance(i));
      prall.push_back(*pr);
    }
  }

  if(sall.size() == 0)
    throw Exception("Measured radiance Spectrum empty, pixel ranges must be empty");

  int nrow = 0;
  int nvar = 0;
  bool have_uncertainty = true;
  BOOST_FOREACH(const Spectrum& s, sall) {
    nrow += s.spectral_domain().data().rows();
    nvar = std::max(nvar, s.spectral_range().data_ad().number_variable());
    if(s.spectral_range().uncertainty().rows() == 0)
      have_uncertainty = false;
  }
  Unit ud, ur;
  if(sall.size() > 0) {
    ud = sall[0].spectral_domain().units();
    ur = sall[0].spectral_range().units();
  }
  ArrayAd<double, 1> sr(nrow, nvar);
  Array<double, 1> sd(nrow);
  Array<double, 1> uncer;
  if(have_uncertainty)
    uncer.resize(nrow);
  for(int i = 0; i < (int) sall.size(); ++i) {
    sd(prall[i]) = sall[i].spectral_domain().data() * 
      FullPhysics::conversion(sall[i].spectral_domain().units(), ud);
    sr.value()(prall[i]) = sall[i].spectral_range().data() *
      FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    if(have_uncertainty)
      uncer(prall[i]) = sall[i].spectral_range().uncertainty() *
	FullPhysics::conversion(sall[i].spectral_range().units(), ur);
    if(nvar > 0)
      sr.jacobian()(prall[i], Range::all()) =
	sall[i].spectral_range().data_ad().jacobian() *
	FullPhysics::conversion(sall[i].spectral_range().units(), ur);
  }
  return Spectrum(SpectralDomain(sd, ud), SpectralRange(sr, ur, uncer));
}
