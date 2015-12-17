#include "level_1b_heritage.h"
#include <boost/foreach.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bHeritage, Level1b)
.def(luabind::constructor<const std::string&, 
     const std::string&,
     const boost::shared_ptr<NoiseModel>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read the given Level 1B heritage file.
//-----------------------------------------------------------------------

Level1bHeritage::Level1bHeritage(const std::string& Sounding_info_file, 
				 const std::string& Spectrum_file,
				 const boost::shared_ptr<NoiseModel>& Nm)
: hmeta(Sounding_info_file), hspec(Spectrum_file), noise_model_(Nm)
{
  number_spectrometer_ = 
    hmeta.value<std::vector<double> >("sounding_latitude").size();
  // Not really a double, but reading of a std::vector<int> has different semantics 
  // (to handle things like the retrieval_indices).
  std::vector<double> sp(hspec.value<std::vector<double> >("start_pixels"));
  BOOST_FOREACH(double& f, sp)
    f -= 1;			// 1 based to zero based.
  for(int i = 0; i < ((int) sp.size()) - 1; ++i)
    r.push_back(Range((int) sp[i], (int) sp[i + 1] - 1));
  r.push_back(Range((int) sp.back(), hspec.data().rows() - 1));
}

SpectralRange Level1bHeritage::radiance(int Spec_index) const
{ 
  range_check(Spec_index, 0, number_spectrometer());
  Array<double, 1> rad(hspec.data("Radiance")(r[Spec_index]));
  Array<double, 1> uncer;
  if(noise_model_)
    uncer.reference(noise_model_->uncertainty(Spec_index, rad));
  return SpectralRange(rad, Unit("W / cm^2 / sr / cm^-1"), uncer);
}
