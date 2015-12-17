#include "chapman_boa_rt.h"
#include "wgs84_constant.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

boost::shared_ptr<RadiativeTransfer> chapman_boa_rt_create
(const boost::shared_ptr<RtAtmosphere>& Atm,
 const blitz::Array<double, 1>& Sza, 
 const SpectralBound& Spec_bound)
{
  const boost::shared_ptr<AtmosphereOco> atm_oco(boost::dynamic_pointer_cast<AtmosphereOco>(Atm));
  return boost::shared_ptr<RadiativeTransfer>(new ChapmanBoaRT(atm_oco, Sza, Spec_bound));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ChapmanBoaRT, RadiativeTransfer)
.scope
[
 luabind::def("create", &chapman_boa_rt_create)
]
REGISTER_LUA_END()
#endif

ChapmanBoaRT::ChapmanBoaRT(const boost::shared_ptr<AtmosphereOco>& Atm,
			   const blitz::Array<double, 1>& Sza) 
  : atm(Atm), sza(Sza)
{
  // Watch atmosphere for changes, so we clear cache if needed.
  atm->add_observer(*this);

  // Cache is state since we have just been initialized
  chapman_cache_stale.resize(Sza.rows());
  chapman_cache_stale = true;

  chapman_boa.resize(Sza.rows());
}

ChapmanBoaRT::ChapmanBoaRT(const boost::shared_ptr<AtmosphereOco>& Atm,
			   const blitz::Array<double, 1>& Sza, 
			   const SpectralBound& Spec_bound) 
  : spec_bound(Spec_bound), atm(Atm), sza(Sza)
{
  // Watch atmosphere for changes, so we clear cache if needed.
  atm->add_observer(*this);

  // Cache is state since we have just been initialized
  chapman_cache_stale.resize(spec_bound.number_spectrometer());
  chapman_cache_stale = true;

  chapman_boa.resize(spec_bound.number_spectrometer());
}

void ChapmanBoaRT::compute_chapman_factors(const int spec_idx) const {

  if(chapman_cache_stale(spec_idx)) {
    // Constants needed by ChapmanBoa
    double rearth = OldConstant::wgs84_a.convert(units::km).value;
    double rfindex_param = 0.000288;

    // Can not use OCO refracrive index if CO2 and H2O are not
    // present in the state structure
    bool can_use_oco_refr = atm->absorber_ptr()->gas_index("CO2") != -1 and 
      atm->absorber_ptr()->gas_index("H2O") != -1;

    // Atmospheric values
    Array<AutoDerivative<double>, 1> height_grid( atm->altitude(spec_idx).convert(units::km).value.to_array() );
    Array<AutoDerivative<double>, 1> press_grid( atm->pressure_ptr()->pressure_grid().convert(units::Pa).value.to_array() );
    Array<AutoDerivative<double>, 1> temp_grid(press_grid.rows());
    Array<AutoDerivative<double>, 1> co2_vmr(press_grid.rows());
    Array<AutoDerivative<double>, 1> h2o_vmr(press_grid.rows());
    for(int i = 0; i < temp_grid.rows(); ++i) {
      temp_grid(i) =
	atm->temperature_ptr()->temperature(AutoDerivativeWithUnit<double>(press_grid(i), units::Pa)).convert(units::K).value;
      
      if(can_use_oco_refr) {
	co2_vmr(i) = atm->absorber_ptr()->absorber_vmr("CO2")->
	  volume_mixing_ratio(press_grid(i));
	h2o_vmr(i) = atm->absorber_ptr()->absorber_vmr("H2O")->
	  volume_mixing_ratio(press_grid(i));
      }
    }

    // Calculate reference wavelengths, if we can..
    double ref_wavelength;
    if(spec_bound.number_spectrometer() > 0)
      ref_wavelength = spec_bound.center(spec_idx, units::micron).value;

    Logger::debug() << "Generating Chapman Factors\n";
    boost::shared_ptr<AtmRefractiveIndex> refr_index;
    Logger::debug() << "Band " << (spec_idx+1) << " using ";
    if(can_use_oco_refr && spec_bound.number_spectrometer() > 0 && 
       OcoRefractiveIndex::wavelength_in_bounds(ref_wavelength)) {
      Logger::debug() << "OCO";
      refr_index.reset(new OcoRefractiveIndex(ref_wavelength, press_grid, temp_grid, co2_vmr, h2o_vmr));
    } else {
      Logger::debug() << "simple";
      refr_index.reset(new SimpleRefractiveIndex(rfindex_param, press_grid, temp_grid));
    }
    Logger::debug() << " refractive index class.\n";

    chapman_boa[spec_idx].reset( new ChapmanBOA(rearth, sza(spec_idx), height_grid, refr_index) ); 

    // Do not need recomputing until updates from Atmosphere class
    chapman_cache_stale(spec_idx) = false;
  }

}

Spectrum ChapmanBoaRT::reflectance(const SpectralDomain& Spec_domain, 
				int Spec_index, bool Skip_jacobian) const
{
  ArrayAd<double, 1> res;
  if(Skip_jacobian)
    res.reference(ArrayAd<double, 1>(stokes(Spec_domain, Spec_index)
				     (Range::all(), 0)));
  else 
    res.reference(stokes_and_jacobian(Spec_domain, Spec_index)
		  (Range::all(), 0));
  return Spectrum(Spec_domain, SpectralRange(res, units::inv_sr));
}

blitz::Array<double, 2> ChapmanBoaRT::stokes(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));
  /// \todo : Make this faster to where it doesn't compute jacobians, right now it just
  /// fails to return the jacobian part and it is computed unncecessarily

  // Compute chapman factors if needed
  compute_chapman_factors(Spec_index);

  Array<double, 1> wn(Spec_domain.wavenumber());
  Array<double, 2> trans(wn.extent(firstDim), number_stokes());
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);

  for(int i = 0; i < wn.rows(); ++i) {
    trans(i, 0) = value( chapman_boa[Spec_index]->transmittance(atm->optical_depth_wrt_state_vector(wn(i), Spec_index).to_array(), 0) );
    if(disp) *disp += 1;
  }

  return trans;
}

ArrayAd<double, 2> ChapmanBoaRT::stokes_and_jacobian(const SpectralDomain& Spec_domain, int Spec_index) const
{
  FunctionTimer ft(timer.function_timer(true));

  // Compute chapman factors if needed
  compute_chapman_factors(Spec_index);

  Array<double, 1> wn(Spec_domain.wavenumber());
  ArrayAd<double, 2> trans_jac(wn.extent(firstDim), number_stokes(), 1);
  boost::shared_ptr<boost::progress_display> disp = progress_display(wn);

  for(int i = 0; i < wn.extent(firstDim); ++i) {
    AutoDerivative<double> trans_wn = chapman_boa[Spec_index]->transmittance(atm->optical_depth_wrt_state_vector(wn(i), Spec_index).to_array(), 0);
    trans_jac.resize_number_variable(trans_wn.number_variable());
    trans_jac(i, 0) = trans_wn;
    if(disp) *disp += 1;
  }

  return trans_jac;
}

void ChapmanBoaRT::print(std::ostream& Os, bool Short_form) const
{
  Os << "ChapmanBoaRT";
  OstreamPad opad(Os, "  ");
  if(!Short_form) {
    Os << "\nAtmosphere:\n";
    opad << *atm << "\n";
    opad.strict_sync();
  }
}
