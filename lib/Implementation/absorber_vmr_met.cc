#include "absorber_vmr_met.h"
#include "ostream_pad.h"
#include "old_constant.h"
#include <boost/make_shared.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrMet, AbsorberVmr)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Pressure>&,
			  double, 
			  bool,
			  const std::string&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

AbsorberVmrMet::AbsorberVmrMet
(const boost::shared_ptr<Meteorology>& Met_file,
 const boost::shared_ptr<Pressure>& Press,
 double Scale,                         
 bool Scale_flag,
 const std::string& Gas_name)
: AbsorberVmrScaled(Press, Scale, Scale_flag, Gas_name), met(Met_file)
{
  std::string gname = gas_name();
  boost::to_upper(gname);
  if (gname != "H2O" &&
      gname != "O3") {
    std::stringstream err_msg; 
    err_msg << "Only H2O and O3 is supported by AbsorberVmrMet, unknown absorber: " 
	    << gas_name();
    throw Exception(err_msg.str());
  }
}

blitz::Array<double, 1> AbsorberVmrMet::specific_humidity() const
{ 
  return met->specific_humidity();
}

blitz::Array<double, 1> AbsorberVmrMet::vmr_profile() const
{ 
  return met->vmr(gas_name());
}

blitz::Array<double, 1> AbsorberVmrMet::pressure_profile() const
{ 
  return met->pressure_levels();
}

boost::shared_ptr<AbsorberVmr> AbsorberVmrMet::clone
(const boost::shared_ptr<Pressure>& Press) const
{
  return boost::shared_ptr<AbsorberVmr>
    (new AbsorberVmrMet(met, Press, coeff(0).value(),used_flag(0),
			  gas_name()));
}

void AbsorberVmrMet::print(std::ostream& Os) const
{ 
  OstreamPad opad(Os, "    ");
  Os << "AbsorberVmrMet:\n"
     << "  Gas name:       " << gas_name() << "\n"
     << "  Scale:          " << scale_factor() << "\n"
     << "  Retrieval flag: " << (used_flag_value()(0) ? 
					"True\n" : "False\n")
     << "  Meteorology:\n";
  opad << *met << "\n";
  opad.strict_sync();
}
