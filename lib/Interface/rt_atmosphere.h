#ifndef RT_ATMOSPHERE_H
#define RT_ATMOSPHERE_H
#include "state_vector.h"
#include "observer.h"
#include "array_ad_with_unit.h"
#include "accumulated_timer.h"
#include "ground.h"

namespace FullPhysics {
/****************************************************************//**
  This class is responsible for setting up the atmosphere and ground
  information needed to run the Radiative transfer code.

  There are many, many properties associated with the atmosphere. This
  class is not meant to model these properties, it is really the very
  limited information needed to run the Radiative transfer code.

  Note that this includes both the atmosphere and surface parameters
  needed by the RT code.

  The calculation of the Jacobians in LIDORT takes a time directly 
  proportional to the number of variables we are taking the Jacobian 
  with respect to, we use an "intermediate" set of variables for some
  of the reported gradients (e.g., AtmosphereOco uses taur, taug, and 
  tau for each of the aerosol). To support future Atmosphere classes, we
  are purposely vague on exactly what these intermediate variables are, at
  least through the RtAtmosphere interface. The
  "intermediate_variable"  
  function can be used to get the value of these intermediate
  variables and Jacobian with the state vector variables.

  A description of the intermediate variables can be found in 
  doc/LIDORT_Jacobian.pdf.

  Note that it is assumed by the LSI that averaging these intermediate
  variables to get average optical properties makes sense. This is
  true if the variables are taur etc., but might not be true in
  general. If we add a class derived from RtAtmosphere where this
  doesn't make sense, we will need to rework this interface.

  Other objects may depend on the RtAtmosphere, and should be updated
  when the RtAtmosphere is updated. To facilitate that, this class is
  an Oberverable, and objects can add themselves as Observers to be
  notified when the RtAtmosphere is updated.

  Because the absorber calculation tends to be a bottle neck, we keep
  a timer in this class. This class keeps track of the time used in
  the atmosphere calculations. Other classes can make use of
  this information for logging if desired.
*******************************************************************/
class RtAtmosphere : virtual public StateVectorObserver,
	             public Observable<RtAtmosphere> {
public:
  virtual ~RtAtmosphere() {}
  static AccumulatedTimer timer;
  virtual void add_observer(Observer<RtAtmosphere>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<RtAtmosphere>& Obs) 
  { remove_observer_do(Obs, *this);}
  virtual std::string timer_info() const;

//-----------------------------------------------------------------------
/// Number of layers we currently have.
//-----------------------------------------------------------------------

  virtual int number_layer() const = 0;

//-----------------------------------------------------------------------
/// Number of spectrometers we have.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// Altitude grid for current pressure grid. 
//-----------------------------------------------------------------------

  virtual ArrayAdWithUnit<double, 1> altitude(int spec_index) 
    const = 0;

//-----------------------------------------------------------------------
/// Total column optical depth for the given gas. This is 0 if the band
/// isn't one that sees that gas.
//-----------------------------------------------------------------------

  virtual AutoDerivative<double>
  column_optical_depth(double wn, int spec_index, const std::string& Gas_name) const = 0;

//-----------------------------------------------------------------------
/// The optical depth for each layer, for the given wave number.
///
/// The derivatives of the optical depth are with respect to the 
/// intermediate variables, rather than the state vector variables (see
/// description of RtAtmosphere class for discussion of this).
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \return Optical depth for each layer. This is number_layer() in size
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> 
    optical_depth_wrt_iv(double wn, int spec_index) const = 0;

//-----------------------------------------------------------------------
/// The single scattering albedo for each layer, for the given wave number.
///
/// The derivatives of the optical depth are with respect to the 
/// intermediate variables, rather than the state vector variables (see
/// description of RtAtmosphere class for discussion of this).
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \return Single scattering albedo for each layer. This is
/// number_layer() in size 
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> 
    single_scattering_albedo_wrt_iv(double wn, int spec_index) const = 0;

//-----------------------------------------------------------------------
/// The scattering moments for for each layer, for the given wave
/// number.
///
/// The scattering moments use the de Rooij convention for the 6
/// scattering matrix element.
///
/// The derivatives of the optical depth are with respect to the 
/// intermediate variables, rather than the state vector variables (see
/// description of RtAtmosphere class for discussion of this).
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \param nummom Number of moments to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \param numscat Number of scattering matrix elements to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \return Scattering moments for each layer. This is 
///         number_moment + 1 x number_layer() x number scattering
///         matrix elements
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, int nummom = -1, 
			   int numscat = -1) const = 0;

//-----------------------------------------------------------------------
/// This is a variation of optical_depth that takes the supplied value
/// for the intermediate variables rather than calculating it own
/// value. This is used by the LSI to get "average optical properties".
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \param iv Intermediate variable values to use.
/// \return Optical depth for each layer. This is number_layer() in size
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> 
  optical_depth_wrt_iv(double wn, int spec_index, const ArrayAd<double, 2>& iv) 
    const = 0;

//-----------------------------------------------------------------------
/// This is a variation of single_scattering_albedo that takes the
/// supplied value  for the intermediate variables rather than
/// calculating it own value. This is used by the LSI to get "average
/// optical properties". 
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \param iv Intermediate variable values to use.
/// \return Single scattering albedo for each layer. This is
/// number_layer() in size 
//-----------------------------------------------------------------------
  
  virtual ArrayAd<double, 1> 
    single_scattering_albedo_wrt_iv(double wn, int spec_index, 
				    const ArrayAd<double, 2>& iv) const = 0;

//-----------------------------------------------------------------------
/// This is a variation of scattering_moment that takes the
/// supplied value  for the intermediate variables rather than
/// calculating it own value. This is used by the LSI to get "average
/// optical properties". 
///
/// \param wn The wave number to calculate parameters for.
/// \param spec_index The spectrometer index
/// \param iv Intermediate variable values to use.
/// \param nummom Number of moments to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \param numscat Number of scattering matrix elements to include in
///           scatt_mom_each_layer, the default it to include all of
///           them.
/// \return Scattering moments for each layer. This is 
///         number_moment + 1 x number_layer() x number scattering
///         matrix elements
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3>
  scattering_moment_wrt_iv(double wn, int spec_index, 
			   const ArrayAd<double, 2>& iv,
			   int nummom = -1, 
			   int numscat = -1) const = 0;


  ArrayAd<double, 1> 
    optical_depth_wrt_state_vector(double wn, int spec_index) const;
  ArrayAd<double, 1> 
    single_scattering_albedo_wrt_state_vector(double wn, int spec_index) const;
  ArrayAd<double, 3>
  scattering_moment_wrt_state_vector(double wn, int spec_index, 
     int nummom = -1, int numscat = -1) const;

//-----------------------------------------------------------------------
/// This gives the values of the intermediate variables and the
/// Jacobian with respect to the state vector. This is
/// number_layer() x number variables
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2>
    intermediate_variable(double wn, int spec_index) const = 0;

//-----------------------------------------------------------------------
/// Object that represents the ground surface. If null then there
/// is no surface for this atmosphere.
//-----------------------------------------------------------------------

  virtual const boost::shared_ptr<Ground> ground() const = 0;

//-----------------------------------------------------------------------
/// Return true if  we have an atmosphere for uplooking mode, i.e., we
/// don't have a ground defined.
//-----------------------------------------------------------------------

  virtual bool uplooking() const = 0;

//-----------------------------------------------------------------------
/// Reset timer
//-----------------------------------------------------------------------

  virtual void reset_timer() { timer.reset_elapsed(); }
};
}
#endif
