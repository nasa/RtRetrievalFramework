#include "ground_coxmunk_scaled.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<Ground> ground_cm_scaled_output_create(boost::shared_ptr<Ground>& coxmunk, boost::shared_ptr<Ground>& brdf_weight)
{
    return boost::shared_ptr<GroundCoxmunkScaled>
        (new GroundCoxmunkScaled(boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk),
                                 boost::dynamic_pointer_cast<GroundBrdfWeight>(brdf_weight)));
}


REGISTER_LUA_DERIVED_CLASS(GroundCoxmunkScaled, Ground)
.def(luabind::constructor<boost::shared_ptr<GroundCoxmunk>&, boost::shared_ptr<GroundBrdfWeight>&>())
.scope
[
    luabind::def("create", &ground_cm_scaled_output_create)
]
REGISTER_LUA_END()
#endif


GroundCoxmunkScaled::GroundCoxmunkScaled(const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                         const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight)
    : coxmunk_(Coxmunk), brdf_weight_(Brdf_weight)
{
    std::vector<boost::shared_ptr<SubStateVectorObserver> > proxied;
    proxied.push_back(coxmunk_);
    proxied.push_back(brdf_weight_);
    initialize(proxied);
}

ArrayAd<double, 1> GroundCoxmunkScaled::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars = coxmunk_->surface_parameter(wn, spec_index);
    AutoDerivative<double> weight_param = brdf_weight_->surface_parameter(wn, spec_index)(0);
    spars(0) = weight_param;
    return spars;
}

boost::shared_ptr<Ground> GroundCoxmunkScaled::clone() const {
    return boost::shared_ptr<GroundCoxmunkScaled>(new GroundCoxmunkScaled(coxmunk_, brdf_weight_));
}

void GroundCoxmunkScaled::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundCoxmunkScaled:\n";
    coxmunk_->print(opad);
    brdf_weight_->print(opad);
    opad.strict_sync();
}
