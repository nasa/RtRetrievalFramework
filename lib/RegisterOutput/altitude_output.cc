#include "altitude_output.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AltitudeOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<Altitude>&, const boost::shared_ptr<Pressure>&>())
REGISTER_LUA_END()
#endif

blitz::Array<double, 1> altitude_levels(boost::shared_ptr<Altitude>& Alt, boost::shared_ptr<Pressure>& Pres)
{
    blitz::Array<double, 1> pressure_levels(Pres->pressure_grid().convert(units::Pa).value.value());
    blitz::Array<double, 1> altitude(pressure_levels.shape());
    for(int lev_idx = 0; lev_idx < altitude.rows(); lev_idx++) {
        altitude(lev_idx) = Alt->altitude(AutoDerivativeWithUnit<double>(pressure_levels(lev_idx), units::Pa)).convert(units::m).value.value();
    }
    return altitude;
}

void AltitudeOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    // Freeze the pressure state
    boost::shared_ptr<Pressure> pfreeze = pres->clone();
    boost::shared_ptr<Altitude> afreeze = alt->clone();

    boost::function<blitz::Array<double, 1> ()> f = boost::bind(&altitude_levels, afreeze, pfreeze);
    out->register_data_source("/RetrievalResults/vector_altitude_levels_apriori", f);
}

void AltitudeOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    boost::function<blitz::Array<double, 1>()> f = boost::bind(&altitude_levels, alt, pres);
    out->register_data_source("/RetrievalResults/vector_altitude_levels", f);
}

