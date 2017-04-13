#include "absorber_vmr_fixed_level_output.h"
#include "fill_value.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_fixed_level_create
(const boost::shared_ptr<AbsorberVmr>& A,
 const boost::shared_ptr<StateVector>& Sv)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new AbsorberVmrFixedLevelOutput
     (boost::dynamic_pointer_cast<AbsorberVmrFixedLevel>(A), Sv));
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrFixedLevelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &abs_vmr_fixed_level_create)
]
REGISTER_LUA_END()
#endif

class AbsorberVmrFixedLevelOutputHelper {
public:
  AbsorberVmrFixedLevelOutputHelper
  (const boost::shared_ptr<AbsorberVmrFixedLevel>& Absvmr,
   const boost::shared_ptr<StateVector>& Sv)
    : absvmr(Absvmr), sv (Sv) {}
  Array<double, 2> vmr_covariance() const
  {
    firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
    ArrayAd<double, 1> p = absvmr->pressure()->pressure_grid().value;
    Array<AutoDerivative<double>, 1> vmrv_t(p.rows());
    for(int i = 0; i < p.rows(); ++i)
      vmrv_t(i) = absvmr->volume_mixing_ratio(p(i));
    ArrayAd<double, 1> vmrv(vmrv_t);
    Array<double, 2> res(vmrv.rows(), vmrv.rows());
    if(vmrv.is_constant())
      res = 0;			// Normal only encounter this in
  				// testing, when we haven't yet set up
  				// a covariance matrix and state vector.
    else {
      Array<double, 2> dvmr_dstate(vmrv.jacobian());
      Array<double, 2> cov(sv->state_covariance());
      Array<double, 2> t(dvmr_dstate.rows(), cov.cols()) ;
      t = sum(dvmr_dstate(i1, i3) * cov(i3, i2), i3);
      res = sum(t(i1, i3) * dvmr_dstate(i2, i3), i3);
    }
    return res;
  }
  Array<double, 1> vmr_uncertainty() const
  {
    Array<double, 2> cov(vmr_covariance());
    Array<double, 1> res(cov.rows());
    for(int i = 0; i < res.rows(); ++i)
      res(i) = (cov(i, i) > 0 ? sqrt(cov(i, i)) : 0.0);
    return res;
  }
private:
  boost::shared_ptr<AbsorberVmrFixedLevel> absvmr;
  boost::shared_ptr<StateVector> sv;
};


// See base class for description

void AbsorberVmrFixedLevelOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the pressure state
  boost::shared_ptr<AbsorberVmrFixedLevel> afreeze = 
    boost::dynamic_pointer_cast<AbsorberVmrFixedLevel>(a->clone());
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  out->register_data_source_pad
    ("/RetrievalResults/" + gname + "_profile_apriori", 
     &AbsorberVmrFixedLevel::volume_mixing_ratio_active_level, 
     afreeze, num_level, fill_value<double>());
}

void AbsorberVmrFixedLevelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  std::string gname = a->gas_name();
  boost::to_lower(gname);
  boost::shared_ptr<AbsorberVmrFixedLevelOutputHelper> h
    (new AbsorberVmrFixedLevelOutputHelper(a, sv));
  out->register_data_source_pad
    ("/RetrievalResults/" + gname + "_profile", 
     &AbsorberVmrFixedLevel::volume_mixing_ratio_active_level, 
     a, num_level, fill_value<double>());
  out->register_data_source_pad
    ("/RetrievalResults/" + gname + "_profile_uncert", 
     &AbsorberVmrFixedLevelOutputHelper::vmr_uncertainty, h, num_level, fill_value<double>());
  out->register_data_source_pad
    ("/RetrievalResults/" + gname + "_profile_covariance_matrix", 
     &AbsorberVmrFixedLevelOutputHelper::vmr_covariance, h, num_level, fill_value<double>());
}
