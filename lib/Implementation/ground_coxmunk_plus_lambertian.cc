#include "ground_coxmunk_plus_lambertian.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<Ground> ground_cm_plus_lamb_output_create(boost::shared_ptr<Ground>& coxmunk, boost::shared_ptr<Ground>& lambertian)
{
    return boost::shared_ptr<GroundCoxmunkPlusLambertian>
        (new GroundCoxmunkPlusLambertian(boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk),
                                         boost::dynamic_pointer_cast<GroundLambertian>(lambertian)));
}


REGISTER_LUA_DERIVED_CLASS(GroundCoxmunkPlusLambertian, Ground)
.def(luabind::constructor<boost::shared_ptr<GroundCoxmunk>&, boost::shared_ptr<GroundLambertian>&>())
.scope
[
    luabind::def("create", &ground_cm_plus_lamb_output_create)
]
REGISTER_LUA_END()
#endif


GroundCoxmunkPlusLambertian::GroundCoxmunkPlusLambertian(const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
                                                         const boost::shared_ptr<GroundLambertian>& Lambertian)
    : coxmunk_(Coxmunk), lambertian_(Lambertian)
{
    std::vector<boost::shared_ptr<SubStateVectorObserver> > proxied;
    proxied.push_back(coxmunk_);
    proxied.push_back(lambertian_);
    initialize(proxied);
}

ArrayAd<double, 1> GroundCoxmunkPlusLambertian::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars = coxmunk_->surface_parameter(wn, spec_index);
    AutoDerivative<double> lamb_param = lambertian_->surface_parameter(wn, spec_index)(0);
    spars(3) = lamb_param;
    return spars;
}

boost::shared_ptr<Ground> GroundCoxmunkPlusLambertian::clone() const {
    return boost::shared_ptr<GroundCoxmunkPlusLambertian>(new GroundCoxmunkPlusLambertian(coxmunk_, lambertian_));
}

void GroundCoxmunkPlusLambertian::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundCoxmunkPlusLambertian:\n";
    coxmunk_->print(opad);
    lambertian_->print(opad);
    opad.strict_sync();
}
