#include "absorber_vmr.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AbsorberVmr)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AbsorberVmr> >,
			VectorAbsorberVmr)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<AbsorberVmr> >::push_back)
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
