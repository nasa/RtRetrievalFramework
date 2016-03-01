#include "reference_vmr_apriori.h"
#include "linear_interpolate.h"

using namespace FullPhysics;
using namespace blitz;

ReferenceVmrApriori::ReferenceVmrApriori(const blitz::Array<double, 1>& Model_altitude,
                                         const blitz::Array<double, 1>& Model_temperature,
                                         const blitz::Array<double, 1>& Ref_altitude,
                                         const double Ref_latitude,
                                         const double Ref_tropopause_altitude,
                                         const double Obs_latitude)
: model_altitude(Model_altitude), model_temperature(Model_temperature),
  ref_altitude(Ref_altitude), ref_latitude(Ref_latitude), ref_tropopause_altitude(Ref_tropopause_altitude),
  obs_latitude(Obs_latitude)
{
    /*Array<double, 1> mod_grid_vmr = resample_ref_vmrs_to_model_grid(Ref_vmr);
    ap_vmr.reference(mod_grid_vmr);
    resample_at_effective_altitudes();
    apply_latitude_gradients();
    apply_secular_trends();
    apply_seasonal_cycle();*/
}


//-----------------------------------------------------------------------
/// Calculate the tropopause altitude for the model data used to 
/// initalize the class.
///
/// Which is the first instance where the lapse rate exceeds -2K/km
//-----------------------------------------------------------------------

double ReferenceVmrApriori::model_tropopause_altitude() const
{
    // From original code
    const double lapse_rate_threshold = -2;
    const double radius = 6378.137;  // Equatorial radius (km)

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

void ReferenceVmrApriori::apply_latitude_gradients()
{
}

void ReferenceVmrApriori::apply_secular_trends()
{
}

void ReferenceVmrApriori::apply_seasonal_cycle()
{
}

