#ifndef ATMOSPHERE_OCO_H
#define ATMOSPHERE_OCO_H
#include "rt_atmosphere.h"
#include "absorber.h"
#include "pressure.h"
#include "temperature.h"
#include "aerosol.h"
#include "ground.h"
#include "relative_humidity.h"
#include "rayleigh.h"
#include "rayleigh_greek_moment.h"
#include "array_ad_cache.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the atmosphere portion of the state, and uses
  this to set up the atmosphere and ground information needed to run
  the Radiative transfer code.

  This particular implementation forwards most of the work to other
  classes such as Absorber and Aerosol. This class then coordinates
  these other classes, and provides the calculations needed to set up
  the RT code.

  For some set ups, aerosol_ptr and/or ground_ptr may be null. For a
  Rayleigh only atmosphere, we don't have any aerosol to include. For
  up looking (e.g., TCCON FTS), there is no ground portion included in
  the radiative transfer.

  To speed up the calculation of the Jacobian in LIDORT, we make use
  of "intermediate" variables instead of directly using state vector
  variables. A description of this in more detail can be found in 
  doc/LIDORT_Jacobian.pdf
*******************************************************************/
class AtmosphereOco : public RtAtmosphere,
		      public Observer<Aerosol>, 
		      public Observer<Pressure> {
public:
  AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
		const boost::shared_ptr<Pressure>& pressurev,
		const boost::shared_ptr<Temperature>& temperaturev,
		const boost::shared_ptr<Aerosol>& aerosolv,
		const boost::shared_ptr<RelativeHumidity>& rhv,
		const boost::shared_ptr<Ground>& groundv,
		const std::vector<boost::shared_ptr<Altitude> >& altv,
		const boost::shared_ptr<Constant>& C);
  AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
		const boost::shared_ptr<Pressure>& pressurev,
		const boost::shared_ptr<Temperature>& temperaturev,
		const boost::shared_ptr<Aerosol>& aerosolv,
		const boost::shared_ptr<RelativeHumidity>& rhv,
		const std::vector<boost::shared_ptr<Altitude> >& altv,
		const boost::shared_ptr<Constant>& C);
  AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
		const boost::shared_ptr<Pressure>& pressurev,
		const boost::shared_ptr<Temperature>& temperaturev,
		const boost::shared_ptr<RelativeHumidity>& rhv,
		const boost::shared_ptr<Ground>& groundv,
		const std::vector<boost::shared_ptr<Altitude> >& altv,
		const boost::shared_ptr<Constant>& C);
  AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
		const boost::shared_ptr<Pressure>& pressurev,
		const boost::shared_ptr<Temperature>& temperaturev,
		const boost::shared_ptr<RelativeHumidity>& Rh,
		const std::vector<boost::shared_ptr<Altitude> >& altv,
		const boost::shared_ptr<Constant>& C);
  virtual ~AtmosphereOco() {}
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) const;
  virtual int number_spectrometer() const {return alt.size();}
  virtual int number_layer() const 
  {
    if(nlay < 0)
      nlay = pressure->number_layer();
    return nlay;
  }
  virtual AutoDerivative<double> 
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const
  {
    if (totaltaug_cache) {
      if (not totaltaug_cache->is_valid(wn))
        fill_cache(wn, spec_index);

      return (*totaltaug_cache)[wn](absorber->gas_index(Gas_name));
    } else {
      fill_cache(wn, spec_index);
      return totaltaug(absorber->gas_index(Gas_name));
    }
  }
  virtual boost::shared_ptr<ArrayAdCache<double, double, 1> >& column_optical_depth_cache()
  {
    return totaltaug_cache;
  };
  virtual ArrayAd<double, 1> 
    optical_depth_wrt_iv(double wn, int spec_index) const
  {
    fill_cache(wn, spec_index);
    return tau;
  }
  virtual ArrayAd<double, 1> 
    single_scattering_albedo_wrt_iv(double wn, int spec_index) const
  {
    fill_cache(wn, spec_index);
    return omega;
  }
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, int nummom = -1, 
		    int numscat = -1) const
  {
    fill_cache(wn, spec_index);
    return scattering_moment_common(wn, nummom, numscat);
  }
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_iv(double wn, int spec_index,
		const ArrayAd<double, 2>& iv) const
  {
    FunctionTimer ft(timer.function_timer());
    calc_rt_parameters(wn, iv);
    return tau;
  }
  virtual ArrayAd<double, 1> 
  single_scattering_albedo_wrt_iv(double wn, int spec_index,
			   const ArrayAd<double, 2>& iv) const
  {
    FunctionTimer ft(timer.function_timer());
    calc_rt_parameters(wn, iv);
    return omega;
  }
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, 
		    const ArrayAd<double, 2>& iv,
		    int nummom = -1, 
		    int numscat = -1) const
  {
    FunctionTimer ft(timer.function_timer());
    calc_rt_parameters(wn, iv);
    return scattering_moment_common(wn, nummom, numscat);
  }

  virtual ArrayAd<double, 2>
    intermediate_variable(double wn, int spec_index) const
  {
    fill_cache(wn, spec_index);
    return intermediate_v;
  }

  virtual const boost::shared_ptr<Ground> ground() const {return ground_ptr;}
  virtual bool uplooking() const {return !ground_ptr;}

  // We use attach_notify to directly attach the various object that
  // Atmosphere contains. This means we don't need to do anything with
  // changes to the StateVector in this class, it is already handled
  // by the objects we contain.
  virtual void notify_update(const StateVector& Sv) 
  { notify_update_do(*this); sv_size = (int) Sv.state().size();}

  virtual void print(std::ostream& Os) const;
  virtual void notify_update(const Aerosol& A);
  virtual void notify_update(const Pressure& P)
  { nlay = -1; }

  virtual void reset_timer();
  virtual std::string timer_info() const;

  const boost::shared_ptr<Pressure>& pressure_ptr() const {return pressure;}
  const boost::shared_ptr<Absorber>& absorber_ptr() const {return absorber;}
  const boost::shared_ptr<Temperature>& temperature_ptr() const 
  {return temperature;}
  const boost::shared_ptr<Aerosol>& aerosol_ptr() const {return aerosol;}
  void set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv);
  const boost::shared_ptr<RelativeHumidity>& relative_humidity_ptr() const
  { return rh; }
  const boost::shared_ptr<Constant>& constant_ptr() const {return constant;}
  const boost::shared_ptr<Rayleigh>& rayleigh_ptr() const {return rayleigh;}
  const std::vector<boost::shared_ptr<Altitude> >& altitude_ptr() const
  { return alt;}
  const boost::shared_ptr<Altitude>& altitude_ptr(int Spec_index) const
  { range_check(Spec_index, 0, number_spectrometer()); return alt[Spec_index];}
  boost::shared_ptr<AtmosphereOco> clone() const;

//-----------------------------------------------------------------------
/// Indicate we have rayleigh only atmosphere, i.e., we don't have any
/// aerosol content.
//-----------------------------------------------------------------------

  bool rayleigh_only_atmosphere() const {return !aerosol;}

  void set_surface_pressure_for_testing(double x);
private:

  boost::shared_ptr<Absorber> absorber;
  boost::shared_ptr<Pressure> pressure;
  boost::shared_ptr<Temperature> temperature;
  boost::shared_ptr<Aerosol> aerosol;
  boost::shared_ptr<RelativeHumidity> rh;
  boost::shared_ptr<Ground> ground_ptr;
  boost::shared_ptr<Rayleigh> rayleigh;
  boost::shared_ptr<Constant> constant;

  // The Altitude and Gravity constants depend on the specific
  // spectrometer we are using, because they see different ground
  // locations and this location enters into the hydrostatic equations
  // that determine the gravity and altitude constants.
  std::vector<boost::shared_ptr<Altitude> > alt;
  int sv_size;

  // We cache the last calculation of these values
  // Keeping these around to reuse helps since
  // its almost certain each subsequent calculation
  // will have the same size, meaning the .resize operation
  // has no work to do and no new memory has to be allocated
  mutable double wn_tau_cache;
  mutable int spec_index_tau_cache;
  mutable ArrayAd<double, 1> omega, taur, frac_ray, tau;
  mutable ArrayAd<double, 1> totaltaug;
  mutable ArrayAd<double, 2> frac_aer;
  mutable ArrayAd<double, 2> intermediate_v;
  mutable int nlay;

  // Items that might need to be cached for access
  // outside of say the radiative transfer loop, without
  // causing everything else to be recalculated
  boost::shared_ptr<ArrayAdCache<double, double, 1> > totaltaug_cache;

  void initialize();
  void calc_rt_parameters(double wn, const ArrayAd<double, 2>& iv) const;
  void calc_intermediate_variable(double wn, int spec_index) const;
  bool fill_cache(double wn, int spec_index) const;
  ArrayAd<double, 3> scattering_moment_common(double wn,
    int nummom, int numscat) const;
};
}
#endif
