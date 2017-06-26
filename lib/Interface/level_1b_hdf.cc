#include "level_1b_hdf.h"
#include "ostream_pad.h"
#include "fp_exception.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<NoiseModel> level_1b_hdf_noise_model_get(const Level1bHdf& F)
{ return F.noise_model(); }
void level_1b_hdf_noise_model_set(Level1bHdf& F, const boost::shared_ptr<NoiseModel>& Nm)
{ F.noise_model(Nm); }
REGISTER_LUA_DERIVED_CLASS(Level1bHdf, Level1b)
.property("noise_model", 
	  &level_1b_hdf_noise_model_get,
	  &level_1b_hdf_noise_model_set)
REGISTER_LUA_END()
#endif

Level1bHdf::Level1bHdf(const std::string& Fname, 
		       const boost::shared_ptr<HdfSoundingId>& Sounding_id)
: file_name(Fname), hdf_sounding_id_(Sounding_id)
{
  hfile.reset(new HdfFile(Fname));
}

Level1bHdf::Level1bHdf(const boost::shared_ptr<HdfFile>& Hfile, 
		       const boost::shared_ptr<HdfSoundingId>& Sounding_id)
: file_name(Hfile->file_name()), hfile(Hfile), hdf_sounding_id_(Sounding_id)
{
  // All work done in initialization above
}

void Level1bHdf::print(std::ostream& Os) const
{
  Os << "Level 1B Hdf file:\n"
     << "  File name:   " << file_name << "\n";
  OstreamPad opad(Os, "     ");
  Os << "  Sounding:\n";
  opad << *hdf_sounding_id_;
  opad.strict_sync();
  Os << "  Noise Model:\n";
  if(noise_model_)
    opad << *noise_model_;
  else
    opad << "No noise model";
  opad.strict_sync();
}

// See base class for description.
SpectralRange Level1bHdf::radiance(int Spec_index) const
{
  SpectralRange r = radiance_no_uncertainty(Spec_index);
  if(noise_model_)
    return SpectralRange(r.data(), r.units(), 
			 noise_model_->uncertainty(Spec_index, r.data()));
  else
    return r;
}
