#include "reference_vmr_apriori.h"
#include "wgs84_constant.h"
#include "old_constant.h"
#include "linear_interpolate.h"
#include "fp_logger.h"

using namespace FullPhysics;
using namespace blitz;
using namespace boost::posix_time;

/// Inter-hemispheric gradients
/// Positive values imply more in the NH than the SH than
/// would be expected based solely on the age of the air.
std::map<std::string, double> lat_gradient =
{
    {"CO2",     0.000 },
    {"O3",      0.200 },
    {"CO",      0.280 },
    {"CH4",     0.033 },
    {"NO2",     0.250 },
    {"NH3",     0.200 },
    {"HNO3",    0.100 },
    {"H2CO",    0.200 },
    {"HCN",     0.100 },
    {"CH3F",    0.200 },
    {"CH3Cl",   0.200 },
    {"CF4",     0.200 },
    {"CCl2F2",  0.200 },
    {"CCl3F",   0.200 },
    {"CH3CCl3", 0.200 },
    {"CCl4",    0.200 },
    {"C2H6",    0.300 },
    {"C2H4",    0.300 },
    {"C2H2",    0.300 },
    {"CHClF2",  0.200 },
    {"CH3Br",   0.200 },
    {"HCOOH",   0.200 },
    {"CHCl2F",  0.200 },
    {"SF6",     0.300 },
    {"F113",    0.300 },
    {"F142b",   0.200 },
    {"CH3OH",   0.200 },
    {"CH3CHO",  0.200 },
    {"CH3CN",   0.200 },
    {"NF3",     0.300 },
    {"CHF3",    0.200 },
    {"f141b",   0.200 },
    {"CH3COOH", 0.200 },
    {"C3H8",    0.500 }
};

std::map<std::string, double> secular_trend = 
{
    {"CO2",     0.0052},
    {"N2O",     0.001},
    {"CO",     -0.006},
    {"CH4",     0.003},
    {"HF",     -0.01}
};

// Amplitude at surface at 50N
std::map<std::string, double> seasonal_amplitude =
{
    {"CO2",     0.0081},
    {"CO",      0.25},
};
//
// From original code
const double lapse_rate_threshold = -2;

ReferenceVmrApriori::ReferenceVmrApriori(const blitz::Array<double, 1>& Model_altitude,
                                         const blitz::Array<double, 1>& Model_temperature,
                                         const blitz::Array<double, 1>& Ref_altitude,
                                         const double Ref_latitude,
                                         const Time& Ref_time,
                                         const double Ref_tropopause_altitude,
                                         const double Obs_latitude,
                                         const Time& Obs_time)
: model_altitude(Model_altitude), model_temperature(Model_temperature),
  ref_altitude(Ref_altitude), ref_latitude(Ref_latitude), ref_time(Ref_time), ref_tropopause_altitude(Ref_tropopause_altitude),
  obs_latitude(Obs_latitude), obs_time(Obs_time)
{
}

//-----------------------------------------------------------------------
/// Calculate the tropopause altitude for the model data used to 
/// initalize the class.
///
/// Which is the first instance where the lapse rate exceeds -2K/km
//-----------------------------------------------------------------------

double ReferenceVmrApriori::model_tropopause_altitude() const
{
    // Equatorial radius (km)
    double radius = OldConstant::wgs84_a.convert(units::km).value;

    double last_lr_alt = 0;
    double last_lapse_rate = 0;
    double ztrop = 0;
    for (int i = 1; i < model_altitude.rows(); i++) {
        // Lapse rate
        double lapse_rate = (model_temperature(i) - model_temperature(i - 1)) / 
            (model_altitude(i) - model_altitude(i - 1));

        // Altitude at which lapse rate == lr
        double lr_alt = 0.5 * (model_altitude(i) + model_altitude(i - 1));

        if(model_altitude(i - 1) > 5.0 && lapse_rate > lapse_rate_threshold) {
            ztrop = last_lr_alt + (lr_alt - last_lr_alt) * (-2 - last_lapse_rate) /(lapse_rate - last_lapse_rate);
            ztrop = ztrop / (1 - ztrop / radius);  // convert H to Z
            return ztrop;
        }

        last_lr_alt = lr_alt;
        last_lapse_rate = lapse_rate;
    }
    throw Exception("Could not determine tropopause altitude");
}

//-----------------------------------------------------------------------
/// Computes an altitude grid for resampling that takes into account
/// a difference in the tropopause altitude for the target from the
/// model.
//-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::effective_altitude() const 
{
    double mod_tropo_alt = model_tropopause_altitude();
    Array<double, 1> effective_altitudes(model_altitude.shape());
    for (int lev_idx = 0; lev_idx < model_altitude.rows(); lev_idx++) {
        double zeff;
        if(model_altitude(lev_idx) < mod_tropo_alt) {
            // troposphere
            zeff = model_altitude(lev_idx) * ref_tropopause_altitude / mod_tropo_alt;
        } else { 
            // stratosphere
            zeff = model_altitude(lev_idx) + (ref_tropopause_altitude - mod_tropo_alt) *
                exp(-(model_altitude(lev_idx) - mod_tropo_alt) / 10.0);
            zeff = zeff - exp(-std::pow(obs_latitude / 25.0, 4))
                * 3.5 * mod_tropo_alt * std::pow(model_altitude(lev_idx) / mod_tropo_alt - 1, 2)  
                * exp(-(model_altitude(lev_idx) - mod_tropo_alt) / 9.0);
        }
        
        // Don't let effective altitude exceeed limit of model altitude
        if(zeff > model_altitude(model_altitude.rows() - 1)) {
            zeff = model_altitude(model_altitude.rows() - 1);
        }

        effective_altitudes(lev_idx) = zeff;
    }
    return effective_altitudes;
}

//-----------------------------------------------------------------------
/// Computes the Age of Air at particular location  (Altitude, Latitude)
/// relative to the surface at 50N, where Age=0
//-----------------------------------------------------------------------

const double ReferenceVmrApriori::age_of_air(const double altitude) const
{
    double mod_tropo_alt = model_tropopause_altitude();
    
    double fl = obs_latitude / 22;
    double aoa = 0.313 - 0.085 * exp(-std::pow((obs_latitude - 49) / 18, 2))
        -0.268 * exp(-1.42 * altitude / (altitude + mod_tropo_alt)) * fl / sqrt(1 + std::pow(fl, 2));
    if(altitude > mod_tropo_alt) 
        aoa = aoa + 7.0 * (altitude - mod_tropo_alt) / altitude;
    return aoa;
}

//-----------------------------------------------------------------------
/// Resamples a VMR to the effective model altitude grid that accounts
/// for the difference in tropopause altitudes
//-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::resample_to_model_grid(const blitz::Array<double, 1>& vmr) const
{
    Array<double, 1> eff_altitude(effective_altitude());
    if(vmr.rows() != ref_altitude.rows()) {
        Exception err_msg;
        err_msg << "Altitude size: " << ref_altitude.rows() << ", does not match "
            << "vmr size: " << vmr.rows();
        throw err_msg;
    }

    LinearInterpolate<double, double> mod_grid_interp(ref_altitude.begin(), ref_altitude.end(), vmr.begin());
    Array<double, 1> mod_grid_vmr(eff_altitude.shape());
    for(int lev_idx = 0; lev_idx < eff_altitude.rows(); lev_idx++) {
        mod_grid_vmr(lev_idx) = mod_grid_interp(eff_altitude(lev_idx));
    }
    return mod_grid_vmr;
}

//-----------------------------------------------------------------------
/// Modifies the vmr profiles to account for the difference in 
/// latitude between the observation latitude and the reference latitude
///
/// In the middle stratosphere, gas distributions are assumed symmetrical
/// about equator. At the surface, gas distributions are assumed 
/// anti-symmetric about equator. At intermediate altitudes the profiles 
/// are interpolated between these limiting behaviors.
//-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::apply_latitude_gradient(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const
{
    double mod_tropo_alt = model_tropopause_altitude();

    double gas_gradient;
    try {
        gas_gradient = lat_gradient.at(gas_name);
    } catch(const std::out_of_range& exc) {
        Logger::warning() << "apply_latitude_gradients: No latitude gradient found for: " << gas_name;
        gas_gradient = 0.0;
    }
    double xref = gas_gradient * (ref_latitude / 15) / sqrt(1 + std::pow(ref_latitude/15, 2));
    double xobs = gas_gradient * (obs_latitude / 15) / sqrt(1 + std::pow(obs_latitude/15, 2));

    Array<double, 1> grad_vmr(vmr.shape());
    for(int lev_idx = 0; lev_idx < vmr.rows(); lev_idx++) {
        double frac = 1.0 / (1.0 + std::pow(model_altitude(lev_idx) / mod_tropo_alt, 2));
        grad_vmr(lev_idx) = vmr(lev_idx) * (1 + frac * xobs) / (1 + frac * xref);
    }
    return grad_vmr;
}

///-----------------------------------------------------------------------
/// Modifies the a priori profiles on a gas-by-gas basis to account for 
/// the difference in time between the observation/model and the reference 
/// vmrs. This includes the secular trend but not the seasonal cycle.
///-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::apply_secular_trend(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const
{
    double time_diff = obs_time.frac_year() - ref_time.frac_year();

    double trend;
    try {
        trend = secular_trend.at(gas_name);
    } catch(const std::out_of_range& exc) {
        Logger::warning() << "apply_secular_trend: No secular trend found for: " << gas_name;
        trend = 0.0;
    }

    Array<double, 1> secular_vmr(vmr.shape());
    for(int lev_idx = 0; lev_idx < vmr.rows(); lev_idx++) {
        double aoa = age_of_air(model_altitude(lev_idx));
        secular_vmr(lev_idx) = vmr(lev_idx) * (1 + trend * (time_diff - aoa));
        if(gas_name == "HF") 
            secular_vmr(lev_idx) = secular_vmr(lev_idx) / (1.0 + exp((-time_diff + aoa - 16) / 5.0));
        if(gas_name == "SF6")
            secular_vmr(lev_idx) = 1.5 * secular_vmr(lev_idx) / (1.0 + exp((-time_diff + aoa - 4) / 9.0));
    }

    return secular_vmr;
}

///-----------------------------------------------------------------------
/// Modifies the a priori vmr profile to account for the season of 
/// the observation/model.
///-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::apply_seasonal_cycle(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const
{
    double twopi = 2.0 * OldConstant::pi;

    // Remove year just to get how far into the year
    double obs_year_frac = obs_time.frac_year() - ptime(obs_time).date().year();

    double amplitude;
    try {
        amplitude = seasonal_amplitude.at(gas_name);
    } catch(const std::out_of_range& exc) {
        Logger::warning() << "apply_seasonal_cycle: No seasonal cycle amplitude found for: " << gas_name;
        amplitude = 0.0;
    }

    Array<double, 1> seasonal_vmr(vmr.shape());
    for(int lev_idx = 0; lev_idx < vmr.rows(); lev_idx++) {
        double zobs = model_altitude(lev_idx);
        double aoa = age_of_air(zobs);
        double sca;
        if (gas_name == "CO2") {
            // seasonal variation
            double sv = std::sin(twopi * (obs_year_frac - 0.834 - aoa));
            // seasonal variation
            double svnl = sv + 1.80 * exp(-std::pow((obs_latitude - 74) / 41, 2)) * (0.5 - std::pow(sv, 2));
            sca = svnl * exp(-aoa/0.20) * (1 + 1.33 * exp(-std::pow((obs_latitude - 76) / 48, 2))*(zobs + 6.0) / (zobs + 1.4));
        } else {
            // basic seasonal variation
            double sv = std::sin(twopi * (obs_year_frac - 0.89));
            // latitude dependence 
            double svl = sv * (obs_latitude / 15) / sqrt(1 + std::pow(obs_latitude / 15, 2));
            // altitude dependence
            sca = svl * exp(-aoa / 1.60);
        }
        seasonal_vmr(lev_idx) = vmr(lev_idx) * (1.0 + sca * amplitude);
    }
    return seasonal_vmr;
}

///-----------------------------------------------------------------------
/// Creates the a priori VMR using the various transformation methods
/// of the class.
///-----------------------------------------------------------------------

const blitz::Array<double, 1> ReferenceVmrApriori::apriori_vmr(const blitz::Array<double, 1>& vmr, const std::string& gas_name) const
{
    Array<double, 1> mod_grid_vmr = resample_to_model_grid(vmr);
    Array<double, 1> lat_grad_vmr = apply_latitude_gradient(mod_grid_vmr, gas_name);
    Array<double, 1> secular_trend_vmr = apply_secular_trend(lat_grad_vmr, gas_name);
    Array<double, 1> seasonal_cycle_vmr = apply_seasonal_cycle(secular_trend_vmr, gas_name);
    return seasonal_cycle_vmr;
}
