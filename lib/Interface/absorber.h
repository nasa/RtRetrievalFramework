#ifndef ABSORBER_H
#define ABSORBER_H
#include "absorber_vmr.h"
#include "state_vector.h"
#include "observer.h"
#include "auto_derivative.h"
#include "accumulated_timer.h"
#include <vector>

namespace FullPhysics {
  class Pressure;
  class Temperature;
  class Altitude;
/****************************************************************//**
  This class maintains the absorber portion of the state.

  Other objects may depend on the absorber, and should be updated
  when the absorber is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the absorber is updated.

  Because the absorber calculation tends to be a bottle neck, we keep
  a timer in this class. This class keeps track of the time used in
  the optical_depth_each_layer function. Other classes can make use of
  this information for logging if desired.
*******************************************************************/
class Absorber : virtual public StateVectorObserver,
		    public Observable<Absorber> {
public:
  virtual ~Absorber() {}
  static AccumulatedTimer timer;
  virtual void add_observer(Observer<Absorber>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Absorber>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Number of species.
//-----------------------------------------------------------------------

  virtual int number_species() const 
  { 
    // Not sure of what the issue is, but SWIG 2.0.9 isn't happy with this
    // being a pure virtual. Just return a default value, this function gets
    // overridden by any "real" class derived from Absorber.
    return 0; 
  }

//-----------------------------------------------------------------------
/// Name of gases, in the order that optical_depth_each_layer returns
/// them.
//-----------------------------------------------------------------------

  virtual std::string gas_name(int Species_index) const = 0;

  virtual int gas_index(const std::string& Name) const;

//-----------------------------------------------------------------------
/// This gives the optical depth for each layer, for the given wave
/// number. Note this only includes the Absorbers portion of this,
/// Atmosphere class combines this with Rayleigh and Aerosol
/// scattering.
///
/// This has size of pres->number_active_layer() x number_species()
///
/// We include the derivative of this with respect to the state vector.
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> 
  optical_depth_each_layer(double wn, int spec_index) const = 0;

//-----------------------------------------------------------------------
/// This calculates the gas column, e.g., XCO2. This is the dry air
/// mole fraction of the gas, see section 3.5.4 of the ATB
///
/// We include the derivative of this with respect to the state vector.
//-----------------------------------------------------------------------
  
  virtual AutoDerivative<double> xgas(const std::string& Gas_name) const = 0;

//-----------------------------------------------------------------------
/// Returns the AbsorberVmr object for a given species index.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AbsorberVmr> absorber_vmr(const std::string& gas_name) const = 0;

//-----------------------------------------------------------------------
/// Clone an Absorber object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Absorber>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Absorber object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Absorber> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a Pressure, Altitude and Temperature
/// to use. The intent is that the Pressure, Altitude and Temperature
/// has been cloned from the original Pressure, Altitude and
/// Temperature (although this class has no way to verify this). This
/// allows sets of objects to be cloned using a common Pressure,
/// Altitude and Temperature clones, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Absorber> clone
    (const boost::shared_ptr<Pressure>& Press,
     const boost::shared_ptr<Temperature>& Temp,
     const std::vector<boost::shared_ptr<Altitude> >& Alt) const = 0;

};
}
#endif


