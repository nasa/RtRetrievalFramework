#include "atmosphere_oco.h"
#include "old_constant.h"
#include "linear_algebra.h"
#include "pressure_fixed_level.h"
#include "temperature_fixed_level.h"
#include "ground.h"
#include "ostream_pad.h"
#include <boost/foreach.hpp>
#include <cmath>
using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(AtmosphereOco, RtAtmosphere)
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
	     const boost::shared_ptr<Pressure>&,
	     const boost::shared_ptr<Temperature>&,
	     const boost::shared_ptr<Aerosol>&,
	     const boost::shared_ptr<RelativeHumidity>&,
	     const boost::shared_ptr<Ground>&,
	     const std::vector<boost::shared_ptr<Altitude> >&,
	     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
	     const boost::shared_ptr<Pressure>&,
	     const boost::shared_ptr<Temperature>&,
	     const boost::shared_ptr<RelativeHumidity>&,
	     const boost::shared_ptr<Ground>&,
	     const std::vector<boost::shared_ptr<Altitude> >&,
	     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
	     const boost::shared_ptr<Pressure>&,
	     const boost::shared_ptr<Temperature>&,
	     const boost::shared_ptr<Aerosol>&,
	     const boost::shared_ptr<RelativeHumidity>&,
	     const std::vector<boost::shared_ptr<Altitude> >&,
	     const boost::shared_ptr<Constant>&>())
.def(luabind::constructor<const boost::shared_ptr<Absorber>&,
	     const boost::shared_ptr<Pressure>&,
	     const boost::shared_ptr<Temperature>&,
	     const boost::shared_ptr<RelativeHumidity>&,
	     const std::vector<boost::shared_ptr<Altitude> >&,
	     const boost::shared_ptr<Constant>&>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Index numbers used in Intermediate variables.
//-----------------------------------------------------------------------

const int taug_index = 0;
const int taur_index = 1;
const int taua_0_index = 2;

//-----------------------------------------------------------------------
/// Create an Atmosphere with the given objects.
///
/// Note that it is ok for aerosolv and/or groundv to be null. For
/// up looking mode (e.g, TCCON FTS), we don't include a ground
/// portion. For a pure Rayleigh atmosphere, we don't include any
/// aerosols.
//-----------------------------------------------------------------------

AtmosphereOco::AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
	     const boost::shared_ptr<Pressure>& pressurev,
	     const boost::shared_ptr<Temperature>& temperaturev,
	     const boost::shared_ptr<Aerosol>& aerosolv,
	     const boost::shared_ptr<RelativeHumidity>& rhv,
             const boost::shared_ptr<Ground>& groundv,
	     const std::vector<boost::shared_ptr<Altitude> >& altv,
	     const boost::shared_ptr<Constant>& C)
  : absorber(absorberv), pressure(pressurev), temperature(temperaturev),
    aerosol(aerosolv), rh(rhv), ground_ptr(groundv), constant(C), alt(altv), 
    sv_size(0),
    wn_tau_cache(-1),
    spec_index_tau_cache(-1),
    nlay(-1)
{
  if(!absorber)
    throw Exception("Absorber is not allowed to be null in AtmosphereOco");
  if(!pressure)
    throw Exception("Pressure is not allowed to be null in AtmosphereOco");
  if(!temperature)
    throw Exception("Temperature is not allowed to be null in AtmosphereOco");
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    if(!a)
      throw Exception("Altitude is not allowed to be null in AtmosphereOco");
  // Aerosol and ground actually is allowed to be null, so we don't check it.
  
  initialize();
}

//-----------------------------------------------------------------------
/// Create an Atmosphere with the given objects.
///
/// This is redundant with the other constructors, but Lua doesn't 
/// support passing nulls, so we have versions that leave particular
/// arguments off.
//-----------------------------------------------------------------------

AtmosphereOco::AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
	     const boost::shared_ptr<Pressure>& pressurev,
	     const boost::shared_ptr<Temperature>& temperaturev,
	     const boost::shared_ptr<Aerosol>& aerosolv,
	     const boost::shared_ptr<RelativeHumidity>& rhv,
	     const std::vector<boost::shared_ptr<Altitude> >& altv,
	     const boost::shared_ptr<Constant>& C)
  : absorber(absorberv), pressure(pressurev), temperature(temperaturev),
    aerosol(aerosolv), rh(rhv), constant(C), alt(altv), 
    sv_size(0),
    wn_tau_cache(-1),
    spec_index_tau_cache(-1),
    nlay(-1)
{
  if(!absorber)
    throw Exception("Absorber is not allowed to be null in AtmosphereOco");
  if(!pressure)
    throw Exception("Pressure is not allowed to be null in AtmosphereOco");
  if(!temperature)
    throw Exception("Temperature is not allowed to be null in AtmosphereOco");
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    if(!a)
      throw Exception("Altitude is not allowed to be null in AtmosphereOco");
  // Aerosol and ground actually is allowed to be null, so we don't check it.
  
  initialize();
}

AtmosphereOco::AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
	     const boost::shared_ptr<Pressure>& pressurev,
	     const boost::shared_ptr<Temperature>& temperaturev,
	     const boost::shared_ptr<RelativeHumidity>& rhv,
	     const boost::shared_ptr<Ground>& groundv,
             const std::vector<boost::shared_ptr<Altitude> >& altv,
	     const boost::shared_ptr<Constant>& C)
  : absorber(absorberv), pressure(pressurev), temperature(temperaturev),
    rh(rhv), ground_ptr(groundv), constant(C), alt(altv), 
    sv_size(0),
    wn_tau_cache(-1),
    spec_index_tau_cache(-1),
    nlay(-1)
{
  if(!absorber)
    throw Exception("Absorber is not allowed to be null in AtmosphereOco");
  if(!pressure)
    throw Exception("Pressure is not allowed to be null in AtmosphereOco");
  if(!temperature)
    throw Exception("Temperature is not allowed to be null in AtmosphereOco");
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    if(!a)
      throw Exception("Altitude is not allowed to be null in AtmosphereOco");
  // Aerosol and ground actually is allowed to be null, so we don't check it.
  
  initialize();
}

AtmosphereOco::AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
	     const boost::shared_ptr<Pressure>& pressurev,
	     const boost::shared_ptr<Temperature>& temperaturev,
	     const boost::shared_ptr<RelativeHumidity>& rhv,
	     const std::vector<boost::shared_ptr<Altitude> >& altv,
	     const boost::shared_ptr<Constant>& C)
  : absorber(absorberv), pressure(pressurev), temperature(temperaturev),
    constant(C), rh(rhv), alt(altv), sv_size(0), wn_tau_cache(-1), 
    spec_index_tau_cache(-1),
    nlay(-1)
{
  if(!absorber)
    throw Exception("Absorber is not allowed to be null in AtmosphereOco");
  if(!pressure)
    throw Exception("Pressure is not allowed to be null in AtmosphereOco");
  if(!temperature)
    throw Exception("Temperature is not allowed to be null in AtmosphereOco");
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    if(!a)
      throw Exception("Altitude is not allowed to be null in AtmosphereOco");
  // Aerosol and ground actually is allowed to be null, so we don't check it.
  
  initialize();
}

void AtmosphereOco::initialize()
{
  rayleigh.reset(new Rayleigh(pressure, alt, *constant));
  if(aerosol)
    aerosol->add_observer(*this);
  pressure->add_observer(*this);

  column_optical_depth_cache().reset( new ArrayAdMapCache<double, double, 1>() );
}

void AtmosphereOco::notify_add(StateVector& Sv)
{
  // Attach each of the classes we use directly. We then don't need to
  // do anything in this class to handle state vector updates, we rely
  // on absorber, pressure, temperature and aerosol to handle their
  // own state.
  Sv.add_observer(*absorber);
  Sv.add_observer(*pressure);
  Sv.add_observer(*temperature);
  if(aerosol) {
    Sv.add_observer(*aerosol);
  }
  if(ground_ptr) {
    Sv.add_observer(*ground_ptr);
  }
  sv_size = Sv.state().size();
}

void AtmosphereOco::notify_remove(StateVector& Sv)
{
  // Attach each of the classes we use directly. We then don't need to
  // do anything in this class to handle state vector updates, we rely
  // on absorber, pressure, temperature and aerosol to handle their
  // own state.
  Sv.remove_observer(*absorber);
  Sv.remove_observer(*pressure);
  Sv.remove_observer(*temperature);
  if(aerosol) {
    Sv.remove_observer(*aerosol);
  }
  if(ground_ptr) {
    Sv.remove_observer(*ground_ptr);
  }
}

//-----------------------------------------------------------------------
/// Changes the aerosol class used. Notifies the relevant obsevers 
/// and invalidates the cache
//-----------------------------------------------------------------------

void AtmosphereOco::set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv) {
    // Remove observers from old aerosol instance
    Sv.remove_observer(*aerosol);

    // Switch to new instance
    aerosol = new_aerosol;

    // Register observers
    aerosol->add_observer(*this);
    Sv.add_observer(*aerosol);

    // Run notify update so gradient is set up correctly
    aerosol->notify_update(Sv);

    // Invalidate caches
    wn_tau_cache = -1;
    spec_index_tau_cache = -1;
}

//-----------------------------------------------------------------------
/// Most of the calculation of each variation of scattering_moment is
/// the same. This does the common part of the calculation, picking up
/// after frac_aer and frac_ray have been filled in.
//-----------------------------------------------------------------------

ArrayAd<double, 3> 
AtmosphereOco::scattering_moment_common(double wn, int nummom, 
					int numscat) const
{
  FunctionTimer ft(timer.function_timer());
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  Array<double, 2> coefsr(RayleighGreekMoment::array());
  Range ra(Range::all());

  // If we don't have any Aerosol (i.e., a Rayleigh atmosphere only),
  // then calculating the scattering matrix is considerably
  // simpler. 
  //
  // This special case is simple enough that it is worth while
  // handling separately. 
  if(rayleigh_only_atmosphere()) {
    Range r1, r2;
    int s2;
    if(nummom == -1 || nummom > coefsr.rows() - 1)
      r1 = Range(0, coefsr.rows() - 1);
    else
      r1 = Range(0, nummom);
    if(numscat == -1) {
      r2 = ra;
      s2 = coefsr.cols();
    } else {
      r2 = Range(0, numscat - 1);
      s2 = numscat;
    }
    int s1 = (nummom == -1 ? coefsr.rows() : nummom + 1);

    // number of variables in rayleigh only case should
    // be at least 2 to match how tau and omega jacobians
    // are set up, see resizing of intermediate_v variable
    ArrayAd<double, 3> res(s1, number_layer(), s2, 2);
    if(s1 > coefsr.rows())
      res = 0;
    for(int i = 0; i < number_layer(); ++i)
      res(r1, i, r2) = coefsr(r1, r2);
    return res;
  } else {
    // Handle case where we do have aerosol.
    Range r1(0, coefsr.rows() - 1);
    Range ra = Range::all();
    ArrayAd<double, 3> 
      pf(aerosol->pf_mom(wn, frac_aer, nummom, numscat));
    Range r2(0, coefsr.rows() - 1);
    pf.value()(r2,ra,ra) += frac_ray.value()(i2) * coefsr(i1,i3);
    pf.jacobian()(r2,ra,ra,ra) += frac_ray.jacobian()(i2,i4) * coefsr(i1,i3);
    pf(0, ra, 0) = 1;
    return pf;
  } // End handling of atmosphere with aerosol.
}

void AtmosphereOco::reset_timer()
{ 
  RtAtmosphere::reset_timer(); 
  Aerosol::timer.reset_elapsed();
  Absorber::timer.reset_elapsed();
}

std::string AtmosphereOco::timer_info() const
{
  std::ostringstream os;
  os << RtAtmosphere::timer_info() << "\n"
     << "   " << Absorber::timer << "\n"
     << "   " << Aerosol::timer << "\n";
  return os.str();
}

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the aerosol is changed, so we observe aerosol
/// and mark the cache when it changes. 
//-----------------------------------------------------------------------

void AtmosphereOco::notify_update(const Aerosol& A)
{
  wn_tau_cache = -1;
  notify_update_do(*this);
}

//-----------------------------------------------------------------------
/// This handles filling the various cached variables. We first check
/// to see if these values are already cached, if so then we skip the
/// calculation. 
/// Returns true if cache is filled, otherwise returns false
//-----------------------------------------------------------------------

bool AtmosphereOco::fill_cache(double wn, int spec_index) const
{
  if(fabs(wn - wn_tau_cache) < 1e-6 &&
     spec_index == spec_index_tau_cache)
    return true;
  
  // If spectrometer changes then erase totaltaug cache
  if(spec_index != spec_index_tau_cache) {
    if (totaltaug_cache)
      totaltaug_cache->clear();
    spec_index_tau_cache = spec_index;
  }

  wn_tau_cache = wn;
  FunctionTimer ft(timer.function_timer());
  calc_intermediate_variable(wn, spec_index);
  calc_rt_parameters(wn, intermediate_v);
  
  return false;
}

//-----------------------------------------------------------------------
/// This calculates intermediate_v. We also fill in totaltaug, since
/// it naturally falls out of this calculation.
//-----------------------------------------------------------------------

void AtmosphereOco::calc_intermediate_variable(double wn, int spec_index) 
const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  range_check(spec_index, 0, number_spectrometer());

//-----------------------------------------------------------------------
// We set up the intermediate variables to be taug, taur, and taua_i
// for each aerosol. To keep things straight, we create the
// variables taug, taur, and taua_i which just point to the right
// parts of the Intermediate variables.
//-----------------------------------------------------------------------

  intermediate_v.resize(number_layer(), 
			2 + (aerosol ? aerosol->number_particle() : 0),
			sv_size);
  ArrayAd<double, 1> taug(intermediate_v(ra, taug_index));
  ArrayAd<double, 1> taur(intermediate_v(ra, taur_index));
  ArrayAd<double, 2> taua_i;
  if(aerosol) {
    Range taua_r(taua_0_index, taua_0_index + 
		 (aerosol ? aerosol->number_particle() - 1 : 0));
    taua_i.reference(intermediate_v(ra, taua_r));
  }
  
//-----------------------------------------------------------------------
// Now, fill in taur and taug. We also fill in totaltaug which is the
// sum over all the layers for each gas absorber. This is needed for
// calculating column_optical_depth_main_gas() and
// column_optical_depth_water_vapor().
//-----------------------------------------------------------------------

  if(taur.is_constant())
    taur.value() = rayleigh->optical_depth_each_layer(wn, spec_index).value();
  else
    taur = rayleigh->optical_depth_each_layer(wn, spec_index);
  ArrayAd<double, 2> taug_i =
    absorber->optical_depth_each_layer(wn, spec_index);

  totaltaug.resize(taug_i.cols(), taug_i.number_variable());
  totaltaug.value() = sum(taug_i.value()(i2, i1), i2);
  if(!taug_i.is_constant())
    totaltaug.jacobian() = sum(taug_i.jacobian()(i3, i1, i2), i3);
  else
    totaltaug.jacobian() = 0;

  // If cache exists, then store the computed value
  if(totaltaug_cache)
    totaltaug_cache->insert(wn, totaltaug);

  taug.value() = sum(taug_i.value()(i1, i2), i2);
  if(!taug.is_constant() && !taug_i.is_constant())
    taug.jacobian() = sum(taug_i.jacobian()(i1, i3, i2), i3);
  else
    taug.jacobian() = 0;

//-----------------------------------------------------------------------
/// Add in aerosol, if we have any.
//-----------------------------------------------------------------------

  if(aerosol) {
    if(taua_i.is_constant())
      taua_i.value() = aerosol->optical_depth_each_layer(wn).value();
    else
      taua_i = aerosol->optical_depth_each_layer(wn);
  }
}

//-----------------------------------------------------------------------
/// This calculates the optical depth, single scatter albedo, and
/// fractions used to calculate the scattering_moment for the given
/// wave number. 
///
/// These calculations are coupled enough that it makes sense to do them 
/// at the same time.
///
/// \param wn The wave number to calculate parameters for.
/// \param iv The intermediate variables (i.e., taug, taur,
///    taua_i). This can either be the cached data intermediate_v, or
///    something passed into this class.
//-----------------------------------------------------------------------

void AtmosphereOco::calc_rt_parameters
(double wn, const ArrayAd<double, 2>& iv) const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Range ra(Range::all());
  tau.resize(number_layer(), iv.cols());
  omega.resize(tau.shape(), tau.number_variable());
  frac_ray.resize(tau.shape(), tau.number_variable());

  ArrayAd<double, 1> taug(iv(ra, taug_index));
  ArrayAd<double, 1> taur(iv(ra, taur_index));
  tau.value() = taur.value() + taug.value();
  tau.jacobian() = 1;		// We just add all the taus together
				// to get total tau.

  // We need taur with derivatives with respect to the intermediate
  // variables in a few different places, so store here
  Array<double, 2> scratch(taur.rows(), iv.cols());
  scratch = 0;
  scratch(ra, taur_index) = 1;
  ArrayAd<double, 1> taur_wrt_iv(taur.value(), scratch.copy());
  scratch(ra, taur_index) = 0;

//-----------------------------------------------------------------------
/// Add in aerosol, if we have any, and use to finish up tau, omega,
/// frac_aer and frac_ray.
//-----------------------------------------------------------------------

  if(aerosol) {
    Range taua_r(taua_0_index, taua_0_index + 
		 (aerosol ? aerosol->number_particle() - 1 : 0));
    ArrayAd<double, 2> taua_i(iv(ra, taua_r));
    tau.value() += sum(taua_i.value(), i2);
    ArrayAd<double, 2> aersc(taua_i.shape(), iv.cols());
    for(int i = 0; i < taua_i.cols(); ++i) {
      scratch(ra, taua_0_index + i) = 1;
      ArrayAd<double, 1> t(taua_i.value()(Range::all(), i), scratch);
      aersc(ra, i) = aerosol->ssa_each_layer(wn, i, t);
      scratch(ra, taua_0_index + i) = 0;
    }
    frac_aer.resize(aersc.shape(), iv.cols());
    ArrayAd<double, 1> ssasum(taur_wrt_iv.copy());

    // Because this is a bottle neck, we have explicit expressions for
    // the Jacobians here. This is just the simple chain rule 

    ssasum.value() += sum(aersc.value(), i2);
    ssasum.jacobian() += sum(aersc.jacobian()(i1, i3, i2), i3);

    frac_aer.value() = aersc.value() / ssasum.value()(i1);
    frac_aer.jacobian() = aersc.jacobian() / ssasum.value()(i1) -
      aersc.value()(i1,i2) / (ssasum.value()(i1) * ssasum.value()(i1)) * 
      ssasum.jacobian()(i1, i3);

    frac_ray.value() = taur_wrt_iv.value() / ssasum.value();
    frac_ray.jacobian() = taur_wrt_iv.jacobian() / ssasum.value()(i1) -
      taur_wrt_iv.value()(i1) / (ssasum.value()(i1) * ssasum.value()(i1)) *
      ssasum.jacobian();

    omega.value() = ssasum.value() / tau.value();
    omega.jacobian() = ssasum.jacobian() / tau.value()(i1) -
      ssasum.value()(i1) / (tau.value()(i1) * tau.value()(i1)) *
      tau.jacobian();
  } else {
    // Note this is ok that frac_ray and frac_aer don't get updated on
    // this branch of the logic, if you look at scattering_moment you'll see 
    // these value won't get used.
    omega.value() = taur_wrt_iv.value() / tau.value();
    omega.jacobian() = taur_wrt_iv.jacobian() / tau.value()(i1) -
      taur_wrt_iv.value()(i1) / (tau.value()(i1) * tau.value()(i1)) *
      tau.jacobian();
  }
}

//-----------------------------------------------------------------------
/// This clones a Atmosphere object. This is a deep copy, all of the
/// objects that are part of this are cloned also (e.g., Pressure,
/// Temperature). 
///
/// This cloned copy will *not* be attached to any StateVector, nor
/// will any Observer<Atmosphere> objects be attached (although the
/// original object is unchanged). You can attach the clone to any
/// objects you wish to.
///
/// This property is particularly useful to "freeze" the state. For
/// example, if the StateVector was set the apriori state and the
/// Atmosphere attached to the StateVector is cloned, then the cloned
/// version will continue to be the "Apriori atmosphere", even if the
/// StateVector is subsequently changed thus updating the original
/// object.
//-----------------------------------------------------------------------

boost::shared_ptr<AtmosphereOco> AtmosphereOco::clone() const
{
  boost::shared_ptr<Pressure> pressure_clone = pressure->clone();
  boost::shared_ptr<Temperature> temperature_clone = 
    temperature->clone(pressure_clone);
  boost::shared_ptr<Ground> ground_clone;
  if(ground_ptr)
    ground_clone = ground_ptr->clone();
  std::vector<boost::shared_ptr<Altitude> > alt_clone;
  BOOST_FOREACH(const boost::shared_ptr<Altitude>& a, alt)
    alt_clone.push_back(a->clone(pressure_clone, temperature_clone));
  boost::shared_ptr<Absorber> absorber_clone =
    absorber->clone(pressure_clone, temperature_clone, alt_clone);
  boost::shared_ptr<RelativeHumidity> rh_clone =
    rh->clone(absorber_clone, temperature_clone, pressure_clone);
  boost::shared_ptr<Aerosol> aerosol_clone;
  if(aerosol)
    aerosol_clone = aerosol->clone(pressure_clone, rh_clone);

  boost::shared_ptr<AtmosphereOco> res
    (new AtmosphereOco(absorber_clone, pressure_clone, temperature_clone,
		       aerosol_clone, rh_clone, ground_clone, alt_clone, 
		       constant));
  return res;
}

//-----------------------------------------------------------------------
/// For unit test purposes, it is useful to be able to directly change
/// the surface pressure. This is intended just for testing
/// purposes. This only works if the Pressure is a
/// PressureFixedLevel, otherwise it will fail.
//-----------------------------------------------------------------------

void AtmosphereOco::set_surface_pressure_for_testing(double x)
{
  PressureFixedLevel& p = dynamic_cast<PressureFixedLevel&>(*pressure);
  p.set_surface_pressure(x);
}

void AtmosphereOco::print(std::ostream& Os) const
{
  Os << "AtmosphereOco:\n";
  OstreamPad opad(Os, "    ");
  Os << "  Constant:\n";
  opad << *constant << "\n";
  opad.strict_sync();
  Os << "  Absorber:\n";
  opad << *absorber << "\n";
  opad.strict_sync();
  Os << "  Pressure:\n";
  opad << *pressure << "\n";
  opad.strict_sync();
  Os << "  Temperature:\n";
  opad << *temperature << "\n";
  opad.strict_sync();
  Os << "  Aerosol:\n";
  if(aerosol)
    opad << *aerosol << "\n";
  else
    opad << "Rayleigh only, no aerosols\n";
  opad.strict_sync();
  Os << "  Ground:\n";
  if(ground_ptr)
    opad << *ground_ptr << "\n";
  else
    opad << "No ground, uplooking mode\n";
  opad.strict_sync();
  for(int i = 0; i < (int) alt.size(); ++i) {
    Os << "  Altitude[" << i << "]:\n";
    opad << *(alt[i]) << "\n";
    opad.strict_sync();
  }
}

// See bass class for description
ArrayAdWithUnit<double, 1> AtmosphereOco::altitude(int spec_index) const
{
  range_check(spec_index, 0, number_spectrometer());
  ArrayAdWithUnit<double, 1> p = pressure->pressure_grid();
  Array<AutoDerivative<double>, 1> res(p.rows());
  Unit u = alt[spec_index]->altitude(p(0)).units;
  for(int i = 0; i < res.rows(); ++i)
    res(i) = alt[spec_index]->altitude(p(i)).convert(u).value;
  return ArrayAdWithUnit<double, 1>(ArrayAd<double, 1>(res), u);
}

