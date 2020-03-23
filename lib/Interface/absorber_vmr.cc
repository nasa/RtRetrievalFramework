#include "absorber_vmr.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AbsorberVmr)
.def("vmr_grid_value", &AbsorberVmr::vmr_grid_value)
REGISTER_LUA_END()

// typedef to distinguish between copying value or moving value (C++11) push_back prototoypes 
typedef void(std::vector<boost::shared_ptr<AbsorberVmr> >::*pbt1)(
        const std::vector<boost::shared_ptr<AbsorberVmr> >::value_type&);

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AbsorberVmr> >, VectorAbsorberVmr)
.def(luabind::constructor<>())
.def("push_back", ((pbt1) &std::vector<boost::shared_ptr<AbsorberVmr> >::push_back))
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Return the vmr on the pressure grid.
//-----------------------------------------------------------------------

ArrayAd<double, 1> AbsorberVmr::vmr_grid(const Pressure& P) const
{
  ArrayAd<double, 1> pgrid = P.pressure_grid().convert(Unit("Pa")).value;
  blitz::Array<AutoDerivative<double>, 1> res(pgrid.rows());
  for(int i = 0; i < res.rows(); ++i)
    res(i) = volume_mixing_ratio(pgrid(i));
  return ArrayAd<double, 1>(res);
}
