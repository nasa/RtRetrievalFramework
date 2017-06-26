#include "gas_vmr_apriori.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"

REGISTER_LUA_CLASS(GasVmrApriori)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
                          const boost::shared_ptr<Level1b>&,
                          const boost::shared_ptr<Altitude>&,
                          const HdfFile&,
                          const std::string&,
                          const std::string&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
                          const boost::shared_ptr<Level1b>&,
                          const boost::shared_ptr<Altitude>&,
                          const HdfFile&,
                          const std::string&,
                          const std::string&,
                          const int>())
// Expose version which requires a Pressure argument
.def("apriori_vmr", (const blitz::Array<double, 1>(GasVmrApriori::*)(const Pressure&) const) &GasVmrApriori::apriori_vmr)
.def("tropopause_pressure", &GasVmrApriori::tropopause_pressure)
REGISTER_LUA_END()
#endif
    
// temp_avg_window is the number of points before and after the current temperature value
// to average among

GasVmrApriori::GasVmrApriori(const boost::shared_ptr<Meteorology>& Met_file,
                             const boost::shared_ptr<Level1b>& L1b_file,
                             const boost::shared_ptr<Altitude>& Alt,
                             const HdfFile& Hdf_static_input,
                             const std::string& Hdf_group,
                             const std::string& Gas_name,
                             const int temp_avg_window)
{
    // Read pressure and temperature grids
    blitz::Array<double, 1> ecmwf_press = Met_file->pressure_levels();
    model_press.resize(ecmwf_press.rows());
    model_press = ecmwf_press;

    blitz::Array<double, 1> model_temp = Met_file->temperature();

    // Smooth the model temperature out with a simple moving average to reduce problems finding
    // tropopause altitude due to kinks in the data
    blitz::Array<double, 1> smoothed_model_temp(model_temp.rows());
    for (int out_idx = 0; out_idx < model_temp.rows(); out_idx++) {
        double sum = 0;
        int count = 0;
        int avg_beg = std::max(out_idx - temp_avg_window, 0);
        int avg_end = std::min(out_idx + temp_avg_window, model_temp.rows() - 1);
        for(int avg_idx = avg_beg; avg_idx <= avg_end; avg_idx++) {
            sum += model_temp(avg_idx);
            count++;
        }
        smoothed_model_temp(out_idx) = sum/count;
    }

    model_alt.resize(model_press.rows());
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

    ref_apriori.reset(new ReferenceVmrApriori(model_press.reverse(firstDim), model_alt.reverse(firstDim), smoothed_model_temp.reverse(firstDim), ref_altitude, ref_latitude, ref_time, ref_tropo_alt, obs_latitude, obs_time)); 
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

const double GasVmrApriori::tropopause_pressure() const
{
    // Make a copy and reverse to satisfy array increasing ordering and memory contiguous requirements of linear interpolator
    blitz::Array<double, 1> rev_alt(model_alt.shape());
    rev_alt = model_alt.copy().reverse(firstDim);
    blitz::Array<double, 1> rev_press(rev_alt.shape());
    rev_press = model_press.copy().reverse(firstDim);

    LinearInterpolate<double, double> alt_press_interp(rev_alt.begin(), rev_alt.end(), rev_press.begin());
    return alt_press_interp(tropopause_altitude().convert(units::km).value);
}
