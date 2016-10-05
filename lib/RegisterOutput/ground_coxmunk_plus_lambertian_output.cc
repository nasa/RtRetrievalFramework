#include "ground_coxmunk_plus_lambertian_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

boost::shared_ptr<RegisterOutputBase> ground_cm_plus_lamb_output_create(boost::shared_ptr<Ground>& coxmunk, boost::shared_ptr<Ground>& lambertian, const std::vector<std::string>& hdf_band_names)
{
    return boost::shared_ptr<RegisterOutputBase>
        (new GroundCoxmunkPlusLambertianOutput(boost::dynamic_pointer_cast<GroundCoxmunk>(coxmunk), 
                                               boost::dynamic_pointer_cast<GroundLambertian>(lambertian),
                                               hdf_band_names));
}

REGISTER_LUA_DERIVED_CLASS(GroundCoxmunkPlusLambertianOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<GroundCoxmunk>&, 
                          const boost::shared_ptr<GroundLambertian>&,
                          const std::vector<std::string>&>())
.scope
[
    luabind::def("create", &ground_cm_plus_lamb_output_create)
]
REGISTER_LUA_END()
#endif

GroundCoxmunkPlusLambertianOutput::GroundCoxmunkPlusLambertianOutput(
        const boost::shared_ptr<GroundCoxmunk>& Coxmunk,
        const boost::shared_ptr<GroundLambertian>& Lambertian,
        const std::vector<std::string>& Hdf_band_names)
: surface_type("Coxmunk,Lambertian")
{
    coxmunk_output.reset(new GroundCoxmunkOutput(Coxmunk));
    lambertian_output.reset(new GroundLambertianOutput(Lambertian, Hdf_band_names));
}

void GroundCoxmunkPlusLambertianOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    coxmunk_output->register_output_apriori(out);
    lambertian_output->register_output_apriori(out);
}

void GroundCoxmunkPlusLambertianOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    coxmunk_output->register_output(out);
    lambertian_output->register_output(out);

    out->register_data_source("/RetrievalResults/surface_type", surface_type.c_str()); 
}
