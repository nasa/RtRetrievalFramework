#include "absorber_vmr_fixed_level_scaled_output.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_fixed_level_scaled_create
(const boost::shared_ptr<AbsorberVmr>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrFixedLevelScaledOutput
     (boost::dynamic_pointer_cast<AbsorberVmrFixedLevelScaled>(A)));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrFixedLevelScaledOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_fixed_level_scaled_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AbsorberVmrFixedLevelScaledOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrFixedLevelScaled> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrFixedLevelScaled>(a->clone());
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor_apriori", 
     &AbsorberVmrFixedLevelScaled::scale_factor, afreeze);
}

void AbsorberVmrFixedLevelScaledOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor", 
     &AbsorberVmrFixedLevelScaled::scale_factor, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor_uncert", 
     &AbsorberVmrFixedLevelScaled::scale_uncertainty, a);
}
