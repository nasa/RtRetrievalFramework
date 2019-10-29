#include "absorber_vmr_log_level_output.h"
#include <boost/algorithm/string.hpp>
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_log_level_create
(const boost::shared_ptr<AbsorberVmr>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrLogLevelOutput
     (boost::dynamic_pointer_cast<AbsorberVmrLogLevel>(A)));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLogLevelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_log_level_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AbsorberVmrLogLevelOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrLogLevel> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrLogLevel>(a->clone());
  std::string gname = a->gas_name();
  boost::algorithm::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_log_profile_apriori", 
     &AbsorberVmrLogLevel::log_vmr_profile, afreeze);
}

void AbsorberVmrLogLevelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string gname = a->gas_name();
  boost::algorithm::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_log_profile", 
     &AbsorberVmrLogLevel::log_vmr_profile, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_log_profile_uncert", 
     &AbsorberVmrLogLevel::log_vmr_uncertainty, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_log_profile_covariance_matrix", 
     &AbsorberVmrLogLevel::log_vmr_covariance, a);
}
