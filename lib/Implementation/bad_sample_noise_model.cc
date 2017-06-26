#include "bad_sample_noise_model.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(BadSampleNoiseModel, NoiseModel)
.def(luabind::constructor<const boost::shared_ptr<NoiseModel>&,const blitz::Array<double, 2>&, double>())
.def(luabind::constructor<const boost::shared_ptr<NoiseModel>&,const blitz::Array<bool, 2>&, double>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Calculate uncertainty, using underlying model except where we have
/// bad samples, where we set this to a fixed value.
//-----------------------------------------------------------------------

blitz::Array<double, 1> BadSampleNoiseModel::uncertainty
(int Spec_index, const blitz::Array<double, 1>& Radiance) const 
{
  range_check(Spec_index, 0, bad_sample_mask_.rows());
  blitz::Array<double, 1> uncer = 
    underlying_noise_model_->uncertainty(Spec_index, Radiance);
  blitz::Array<double, 1> res(uncer.shape());
  res(Range::all()) = where(bad_sample_mask_(Spec_index, Range::all()), bad_sample_uncer_, uncer);
  return res;
}

void BadSampleNoiseModel::print(std::ostream& Os) const
{
  Os << "BadSampleNoiseModel\n"
     << "  Bad sample uncertainty: " << bad_sample_uncer_ << "\n";
  OstreamPad opad(Os, "     ");
  Os << "  Underlying noise model:\n";
  opad << *underlying_noise_model_;
  opad.strict_sync();
}

