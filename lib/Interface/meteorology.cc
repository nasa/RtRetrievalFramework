#include "meteorology.h"
#include "log_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

typedef Array<double, 1> (Meteorology::*f_val)() const;
typedef Array<double, 1> (Meteorology::*f_interp)(const Array<double, 1>&) const;

typedef Array<double, 1> (Meteorology::*f_vmr)(const std::string&) const;
typedef Array<double, 1> (Meteorology::*f_vmr_interp)(const std::string&, const Array<double, 1>&) const;

REGISTER_LUA_CLASS(Meteorology)
.def("pressure_levels", &Meteorology::pressure_levels)
.def("specific_humidity", ((f_val) &Meteorology::specific_humidity))
.def("specific_humidity", ((f_interp) &Meteorology::specific_humidity))
.def("vmr", ((f_vmr) &Meteorology::vmr))
.def("vmr", ((f_vmr_interp) &Meteorology::vmr))
.def("temperature", ((f_val) &Meteorology::temperature))
.def("temperature", ((f_interp) &Meteorology::temperature))
.def("surface_pressure", &Meteorology::surface_pressure)
.def("windspeed", &Meteorology::windspeed)
REGISTER_LUA_END()
#endif

Array<double, 1> Meteorology::specific_humidity(const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(specific_humidity(), Pressure_level);
}

Array<double, 1> Meteorology::vmr(const std::string& Species, const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(vmr(Species), Pressure_level);
}

Array<double, 1> Meteorology::temperature(const Array<double, 1>& Pressure_level) const
{
    return interpolate_to_grid(temperature(), Pressure_level);
}

Array<double, 1> Meteorology::interpolate_to_grid(const Array<double, 1>& Profile, const Array<double, 1>& Dest_pressure_levels) const
{
    LogLogInterpolate<double, double> prof_interp(pressure_levels().begin(), pressure_levels().end(), Profile.begin());
    Array<double, 1> prof_res(Dest_pressure_levels.shape());
    for(int i = 0; i < Dest_pressure_levels.rows(); ++i) {
        prof_res(i) = prof_interp(Dest_pressure_levels(i));
    }
    return prof_res;
}

