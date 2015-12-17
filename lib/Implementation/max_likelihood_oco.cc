#include <max_likelihood_oco.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(MaxLikelihoodOCO, MaxLikelihood)
.def(luabind::constructor< const boost::shared_ptr<ForwardModel>& >())
REGISTER_LUA_END()
#endif



MaxLikelihoodOCO::MaxLikelihoodOCO(const boost::shared_ptr<ForwardModel>& fm)
  : ModelMeasure(
      fm->measured_radiance_all().spectral_range().data(),
      Array<double, 1>(sqr(fm->measured_radiance_all().spectral_range().uncertainty()))), 
    ModelMeasureOCO(fm)
{}
