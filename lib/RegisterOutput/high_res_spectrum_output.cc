#include "high_res_spectrum_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "oco_forward_model.h"
#include "lsi_rt.h"

boost::shared_ptr<RegisterOutputBase> hr_spec_as_register_output_base(const boost::shared_ptr<HighResSpectrumOutput>& hr_spec) {
  return boost::dynamic_pointer_cast<RegisterOutputBase>(hr_spec);
}

void hr_spec_add_as_oco_fm_observer(HighResSpectrumOutput& hr_spec, OcoForwardModel& oco_fm) {
  oco_fm.add_observer(hr_spec);
}

void hr_spec_add_as_lsi_rt_observer(HighResSpectrumOutput& hr_spec, LsiRt& rt) {
  rt.add_observer(hr_spec);
}

#include "register_lua.h"
REGISTER_LUA_CLASS(HighResSpectrumOutput)
.def(luabind::constructor<>())
.def("as_register_output_base", &hr_spec_as_register_output_base)
.def("add_as_observer", &hr_spec_add_as_oco_fm_observer)
.def("add_as_observer", &hr_spec_add_as_lsi_rt_observer)
REGISTER_LUA_END()
#endif

const Array<double, 1> HighResSpectrumOutput::saved_spectral_domain(const std::string& spectra_name)
{
  Array<double, 1> res(0);
  BOOST_FOREACH(const boost::shared_ptr<Spectrum>& s, saved_spectra[spectra_name]) {
    // Its possible to have empty spectra if a certain index is notified
    // but not the one before it
    if(s) {
      int nrow = s->spectral_domain().data().rows();
      Range save_range(res.rows(), res.rows() + nrow - 1);
      res.resizeAndPreserve(res.rows() + nrow);
      res(save_range) = s->spectral_domain().data();
    }
  }
  return res;
}

const Array<double, 1> HighResSpectrumOutput::saved_spectral_range(const std::string& spectra_name)
{
  Array<double, 1> res(0);
  BOOST_FOREACH(const boost::shared_ptr<Spectrum>& s, saved_spectra[spectra_name]) {
    // Its possible to have empty spectra if a certain index is notified
    // but not the one before it
    if(s) {
      int nrow = s->spectral_range().data().rows();
      Range save_range(res.rows(), res.rows() + nrow - 1);
      res.resizeAndPreserve(res.rows() + nrow);
      res(save_range) = s->spectral_range().data();
    }
  }
  return res;
}

//-----------------------------------------------------------------------
/// On notification this class registers a dataset based on the 
/// name in the NamedSpectrum. If the same name comes in again it
/// would be re-registered. Multiple bands will be saved together
/// an output as one array.
//-----------------------------------------------------------------------

void HighResSpectrumOutput::notify_update(const boost::shared_ptr<FullPhysics::NamedSpectrum>& Obs)
{
  register_named_spectrum(Obs, "grid_", "radiance_");
}

//-----------------------------------------------------------------------
/// Calls update for each NamedSpectum in a notification sending a 
/// vector of NamedSpectrum, such as for delivering stokes values.
//-----------------------------------------------------------------------

void HighResSpectrumOutput::notify_update(const std::vector<boost::shared_ptr<NamedSpectrum> >& Obs) {
  BOOST_FOREACH(boost::shared_ptr<NamedSpectrum> named_spec, Obs)
    register_named_spectrum(named_spec, "", "reflectance_");
}

void HighResSpectrumOutput::register_named_spectrum(const boost::shared_ptr<FullPhysics::NamedSpectrum>& named_spec,
                                                    const std::string& domain_prefix, const std::string& range_prefix) {
  std::vector<boost::shared_ptr<Spectrum> > &name_spec_vec = saved_spectra[named_spec->name()];
  if ((int) name_spec_vec.size() <= named_spec->index())
    name_spec_vec.resize(named_spec->index()+1);

  boost::shared_ptr<Spectrum> spec_copy( new Spectrum(named_spec->spectral_domain().clone(),
                                                      named_spec->spectral_range().clone()) );
  name_spec_vec[named_spec->index()] = boost::shared_ptr<Spectrum>(spec_copy);

  BOOST_FOREACH(boost::shared_ptr<Output> out, output_files) {
    // Only call registration if output object is still a valid pointer
    if (out) {
      boost::function<Array<double, 1>() > f;

      if (domain_prefix.length() > 0) {
        f = boost::bind(&HighResSpectrumOutput::saved_spectral_domain, this, named_spec->name());
        out->register_data_source("/HighResSpectra/" + domain_prefix + named_spec->name(), f);
      }

      if (range_prefix.length() > 0) {
        f = boost::bind(&HighResSpectrumOutput::saved_spectral_range, this, named_spec->name());
        out->register_data_source("/HighResSpectra/" + range_prefix + named_spec->name(), f);
      }
    }
  }
}
