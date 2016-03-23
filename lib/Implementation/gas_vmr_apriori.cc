#include "gas_vmr_apriori.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(GasVmrApriori)
.def(luabind::constructor<const boost::shared_ptr<Ecmwf>&,
                          const boost::shared_ptr<Level1b>&,
                          const boost::shared_ptr<Altitude>&,
                          const HdfFile&,
                          const std::string&,
                          const std::string&>())
// Expose version which requires a Pressure argument
.def("apriori_vmr", (const blitz::Array<double, 1>(GasVmrApriori::*)(const Pressure&) const) &GasVmrApriori::apriori_vmr)
REGISTER_LUA_END()
#endif

GasVmrApriori::GasVmrApriori(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
                             const boost::shared_ptr<Level1b>& L1b_file,
                             const boost::shared_ptr<Altitude>& Alt,
                             const HdfFile& Hdf_static_input,
                             const std::string& Hdf_group,
                             const std::string& Gas_name)
{
    blitz::Array<double, 1> model_temp;
    
    // Read pressure and temperature grids
    Ecmwf_file->temperature_grid(model_press, model_temp);

    blitz::Array<double, 1> model_alt(model_press.rows());
    for(int lev_idx = 0; lev_idx < model_press.rows(); lev_idx++) {
        model_alt(lev_idx) = Alt->altitude(AutoDerivativeWithUnit<double>(model_press(lev_idx), units::Pa)).convert(units::km).value.value();
    }

    // Get reference data from HDF file
    Array<double, 1> ref_altitude = Hdf_static_input.read_field<double, 1>(Hdf_group + "/altitude_grid");
    double ref_latitude = Hdf_static_input.read_field<double>(Hdf_group + "/latitude");
    Time ref_time = Time::parse_time(Hdf_static_input.read_field<std::string>(Hdf_group + "/time_string"));
    double ref_tropo_alt = Hdf_static_input.read_field<double>(Hdf_group + "/tropopause_altitude");
    ref_vmr.reference( Hdf_static_input.read_field<double, 1>(Hdf_group + "/Gas/" + Gas_name + "/vmr") );
    gas_name = Gas_name;

    // Allow values to be in either increasing altitude or increasing pressure order
    bool in_increasing_alt = Hdf_static_input.read_field<int>(Hdf_group + "/increasing_altitude_order");
    if (!in_increasing_alt) {
        // Reverse on a copy to make sure order is correct behind the scenes
        // to satisfy LinearInterpolate
        ref_altitude = ref_altitude.copy().reverse(firstDim);
        ref_vmr = ref_vmr.copy().reverse(firstDim);
    }

    // Get observation values from L1B file
    double obs_latitude = L1b_file->latitude(0).value;
    Time obs_time = L1b_file->time(0);

    ref_apriori.reset(new ReferenceVmrApriori(model_alt.reverse(firstDim), model_temp.reverse(firstDim), ref_altitude, ref_latitude, ref_time, ref_tropo_alt, obs_latitude, obs_time)); 
}

const blitz::Array<double, 1> GasVmrApriori::apriori_vmr() const
{
    Array<double, 1> ap_vmr = ref_apriori->apriori_vmr(ref_vmr, gas_name);
    ap_vmr.reverseSelf(firstDim);
    return ap_vmr;
}

const blitz::Array<double, 1> GasVmrApriori::apriori_vmr(const Pressure& pressure) const
{
    // Make a copy to avoid issues with memory access into the reversed view
    Array<double, 1> model_ap_vmr(model_press.shape());
    model_ap_vmr = apriori_vmr();

    Array<double, 1> press_levels(pressure.pressure_grid().convert(units::Pa).value.value());

    LinearInterpolate<double, double> mod_vmr_interp(model_press.begin(), model_press.end(), model_ap_vmr.begin());

    Array<double, 1> interp_vmr(press_levels.shape());
    for(int lev_idx = 0; lev_idx < interp_vmr.rows(); lev_idx++) {
        interp_vmr(lev_idx) = mod_vmr_interp(press_levels(lev_idx));
    }

    return interp_vmr;
}
