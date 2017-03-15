#include "absorber_vmr_met_output.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_met_create
(const boost::shared_ptr<AbsorberVmr>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrMetOutput
     (boost::dynamic_pointer_cast<AbsorberVmrMet>(A)));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrMetOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_met_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AbsorberVmrMetOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrScaled> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrScaled>(a->clone());
  boost::shared_ptr<AbsorberVmrMet> afreeze_met = 
    boost::dynamic_pointer_cast<AbsorberVmrMet>(afreeze);
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor_apriori", 
     &AbsorberVmrScaled::scale_factor, afreeze);
  out->register_data_source("/RetrievalResults/specific_humidity_profile_met",
			   &AbsorberVmrMet::specific_humidity, afreeze_met);
}

void AbsorberVmrMetOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor", 
     &AbsorberVmrScaled::scale_factor, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor_uncert", 
     &AbsorberVmrScaled::scale_uncertainty, a);
}
