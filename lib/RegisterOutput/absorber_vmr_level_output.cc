#include "absorber_vmr_level_output.h"
#include <boost/algorithm/string.hpp>
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_level_create
(const boost::shared_ptr<AbsorberVmr>& A)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrLevelOutput
     (boost::dynamic_pointer_cast<AbsorberVmrLevel>(A)));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrLevelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_level_create)
]
REGISTER_LUA_END()
#endif

// Helper class that get co2_grad_del
class AbsorberVmrLevelHelper {
public:
  AbsorberVmrLevelHelper(const boost::shared_ptr<AbsorberVmrLevel>& a,
			 const boost::shared_ptr<AbsorberVmrLevel>& afreeze)
    : a_(a), afreeze_(afreeze) { }
//-----------------------------------------------------------------------
/// This term is useful in the warn level determination. The fact that
/// this is level 20 and level 13 is pretty much empirical I think,
/// this originally came from Chris O'Dell's bias term
//-----------------------------------------------------------------------

  double co2_grad_del() const
  {
    blitz::Array<double, 1> p = a_->vmr_profile();
    blitz::Array<double, 1> ap = afreeze_->vmr_profile();
    double res = fill_value<double>();
    // Don't try to calculate this unless we have enough levels.
    if(p.rows() >= 20) 
      res = (p(19) - p(12)) - (ap(19) - ap(12));
    return res;
  }
private:
  boost::shared_ptr<AbsorberVmrLevel>  a_, afreeze_;
};

// See base class for description

void AbsorberVmrLevelOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrLevel> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrLevel>(a->clone());
  std::string gname = a->gas_name();
  boost::algorithm::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_profile_apriori", 
     &AbsorberVmrLevel::vmr_profile, afreeze);
  if(gname == "co2") {
    boost::shared_ptr<AbsorberVmrLevelHelper> h(new AbsorberVmrLevelHelper(a, afreeze));
    out->register_data_source
      ("/RetrievalResults/co2_vertical_gradient_delta",
       &AbsorberVmrLevelHelper::co2_grad_del, h);
  }
}

void AbsorberVmrLevelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string gname = a->gas_name();
  boost::algorithm::to_lower(gname);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_profile", 
     &AbsorberVmrLevel::vmr_profile, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_profile_uncert", 
     &AbsorberVmrLevel::vmr_uncertainty, a);
  out->register_data_source
    ("/RetrievalResults/" + gname + "_profile_covariance_matrix", 
     &AbsorberVmrLevel::vmr_covariance, a);
}
