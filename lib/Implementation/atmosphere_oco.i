// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "atmosphere_oco.h"
#include "sub_state_vector_array.h"
#include "rayleigh.h"
%}
%fp_shared_ptr(FullPhysics::AtmosphereOco);
%base_import(rt_atmosphere)
%base_import(aerosol)
%import "absorber.i"
%import "temperature.i"
%import "rayleigh.i"
%import "ground.i"
%import "relative_humidity.i"
namespace FullPhysics {
class AtmosphereOco;
}

%fp_shared_ptr(FullPhysics::Observer<FullPhysics::AtmosphereOco>);
%fp_shared_ptr(FullPhysics::Observable<FullPhysics::AtmosphereOco>);

namespace FullPhysics {
%template(ObserverAtmosphereOco) FullPhysics::Observer<AtmosphereOco>;

class AtmosphereOco : public RtAtmosphere,
    public Observer<Aerosol> {
public:
  AtmosphereOco(const boost::shared_ptr<Absorber>& absorberv,
	     const boost::shared_ptr<Pressure>& pressurev,
	     const boost::shared_ptr<Temperature>& temperaturev,
	     const boost::shared_ptr<Aerosol>& aerosolv,
	     const boost::shared_ptr<RelativeHumidity>& rhv,
	     const boost::shared_ptr<Ground>& groundv,
	     const std::vector<boost::shared_ptr<Altitude> >& altv,
	     const boost::shared_ptr<Constant>& C);
  virtual void notify_add(StateVector& Sv);
  virtual void notify_remove(StateVector& Sv);
  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) const;
  %python_attribute_derived(number_spectrometer, int)
  %python_attribute_derived(number_layer, int)
  virtual AutoDerivative<double> 
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const;
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_iv(double wn, int spec_index) const;
  virtual ArrayAd<double, 1> 
  single_scattering_albedo_wrt_iv(double wn, int spec_index) const;
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, int nummom = -1, 
		    int numscat = -1) const;
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_iv(double wn, int spec_index,
		const ArrayAd<double, 2>& iv) const;
  virtual ArrayAd<double, 1> 
  single_scattering_albedo_wrt_iv(double wn, int spec_index,
			   const ArrayAd<double, 2>& iv) const;
  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, 
			   const ArrayAd<double, 2>& iv,
			   int nummom = -1, int numscat = -1) const;
  virtual ArrayAd<double, 2>
    intermediate_variable(double wn, int spec_index) const;
  %python_attribute_derived(ground, boost::shared_ptr<Ground>)
  %python_attribute_derived(uplooking, bool)
  virtual void notify_update(const StateVector& Sv);
  virtual void notify_update(const Aerosol& A);
  virtual void reset_timer();
  %python_attribute_derived(timer_info, std::string)
  %python_attribute2(pressure, pressure_ptr, boost::shared_ptr<Pressure>)
  %python_attribute2(absorber, absorber_ptr, boost::shared_ptr<Absorber>)
  %python_attribute2(aerosol, aerosol_ptr, boost::shared_ptr<Aerosol>)
  void set_aerosol(boost::shared_ptr<Aerosol>& new_aerosol, StateVector& Sv);
  %python_attribute2(temperature, temperature_ptr, boost::shared_ptr<Temperature>)
  %python_attribute2(relative_humidity, relative_humidity_ptr, boost::shared_ptr<RelativeHumidity>)
  %python_attribute2(constant, constant_ptr, boost::shared_ptr<Constant>)
  %python_attribute2(rayleigh, rayleigh_ptr, boost::shared_ptr<Rayleigh>)
  %python_attribute2(altitude_obj, altitude_ptr, std::vector<boost::shared_ptr<Altitude> >)
  boost::shared_ptr<AtmosphereOco> clone() const;
  %python_attribute(rayleigh_only_atmosphere, bool)
  void set_surface_pressure_for_testing(double x);
};
}
