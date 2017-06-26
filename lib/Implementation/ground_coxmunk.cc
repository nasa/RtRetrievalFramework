#include "ground_coxmunk.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GroundCoxmunk, Ground)
.def(luabind::constructor<const double,
                          const bool&, 
                          const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

GroundCoxmunk::GroundCoxmunk(const double Windspeed,
                const bool& Ws_flag, 
                const blitz::Array<double, 1>& Refr_index)
    : SubStateVectorArray<Ground>(Windspeed, Ws_flag),
      refractive_index_(Refr_index)
{
}

ArrayAd<double, 1> GroundCoxmunk::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars;
    spars.resize(4, windspeed().number_variable());
    spars(0) = windspeed();
    spars(1) = refractive_index(spec_index);
    spars(2) = 0;
    // 0 = no shadowing, 1 = include shadowing
    spars(3) = 0;
    return spars;
}

const AutoDerivative<double> GroundCoxmunk::windspeed() const 
{ 
    return coefficient()(0); 
}

const double GroundCoxmunk::refractive_index(const int Spec_idx) const 
{ 
    return refractive_index_(Spec_idx); 
}

boost::shared_ptr<Ground> GroundCoxmunk::clone() const {
  return boost::shared_ptr<GroundCoxmunk>(new GroundCoxmunk(coefficient().value()(0), used_flag_value()(0), refractive_index_));
}

std::string GroundCoxmunk::state_vector_name_i(int i) const {
  return "Ground Coxmunk Windspeed";
}

void GroundCoxmunk::print(std::ostream& Os) const
{
  OstreamPad opad(Os, "    ");
  Os << "GroundCoxmunk:\n";
  opad << "Windspeed: " << windspeed().value() << '\n'  
       << "Refractive Index: " << refractive_index_ << '\n';
  opad.strict_sync();
}

void GroundCoxmunk::update_sub_state_hook()
{
  // Check that windspeed isn't negative.
  if(windspeed().value() < 0)
    throw Exception("Windspeed has gone negative");
}
