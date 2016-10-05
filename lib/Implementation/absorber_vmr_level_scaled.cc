#include "absorber_vmr_level_scaled.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevelScaled, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Pressure>&,
			  const blitz::Array<double, 1>&,
			  double, 
			  bool,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrLevelScaled::AbsorberVmrLevelScaled
(const boost::shared_ptr<Pressure>& Press,
 const blitz::Array<double, 1>& Vmr_profile,
 double Scale,                         
 bool Scale_flag,
 const std::string& Gas_name)
: AbsorberVmrScaled(Press, Scale, Scale_flag, Gas_name), vmr_profile_(Vmr_profile)
{
}

blitz::Array<double, 1> AbsorberVmrLevelScaled::vmr_profile() const
{ 
  return vmr_profile_;
}

blitz::Array<double, 1> AbsorberVmrLevelScaled::pressure_profile() const
{ 
  return press->pressure_grid().value.value();
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrLevelScaled::clone
(const boost::shared_ptr<Pressure>& Press) const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrLevelScaled(Press, vmr_profile_, coeff(0).value(),used_flag(0),
				gas_name()));
}

void AbsorberVmrLevelScaled::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrLevelScaled:\n"
     << "  Gas name:       " << gas_name() << "\n"
     << "  Scale:          " << scale_factor() << "\n"
     << "  Retrieval flag: " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  VMR Profile:\n";
  opad << vmr_profile_ << "\n";
  opad.strict_sync();
}
