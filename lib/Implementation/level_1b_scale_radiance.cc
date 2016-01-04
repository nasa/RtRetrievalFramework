#include "level_1b_scale_radiance.h"
#include "fp_exception.h"
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Level1bScaleRadiance, Level1b)
.def(luabind::constructor<const boost::shared_ptr<Level1b> &, const blitz::Array<double, 1>&>())
REGISTER_LUA_END()
#endif

SpectralRange Level1bScaleRadiance::radiance(int Spec_index) const
{
    range_check(Spec_index, 0, (int) scaling.rows());
    SpectralRange rad = l1b->radiance(Spec_index);
    return SpectralRange(Array<double, 1>(rad.data() * scaling(Spec_index)), rad.units(), rad.uncertainty());
}

void Level1bScaleRadiance::print(std::ostream& Os) const
{
  Os << "Level1bScaleRadiance";
}
