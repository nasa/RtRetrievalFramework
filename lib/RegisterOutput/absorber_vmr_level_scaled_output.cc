#include "absorber_vmr_level_scaled_output.h"
#include <boost/algorithm/string/case_conv.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_level_scaled_create
(const boost::shared_ptr<AbsorberVmr>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrLevelScaledOutput
     (boost::dynamic_pointer_cast<AbsorberVmrLevelScaled>(A)));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevelScaledOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_level_scaled_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AbsorberVmrLevelScaledOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrScaled> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrScaled>(a->clone());
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_scale_factor_apriori", 
     &AbsorberVmrScaled::scale_factor, afreeze);
}

void AbsorberVmrLevelScaledOutput::register_output(const boost::shared_ptr<Output>& out) const
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
