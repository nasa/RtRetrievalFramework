#include "hdf_constant.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(HdfConstant, Constant)
.def(luabind::constructor<const boost::shared_ptr<HdfFile>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Read HDF file to get constants.
//-----------------------------------------------------------------------

HdfConstant::HdfConstant(const boost::shared_ptr<HdfFile>& Hdf_file)
{
  // Fill me in
}
