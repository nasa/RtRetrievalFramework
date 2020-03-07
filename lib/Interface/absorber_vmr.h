#ifndef ABSORBER_VMR_H
#define ABSORBER_VMR_H
#include "state_vector.h"
#include "auto_derivative.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This gives the Gas Absorber Volumn mixing ratio for a single gas. 
  This gets used by AbsorberAbsco class.

  When implementing a new class, you almost always will want to derive
  from AbsorberVmrImpBase rather than from this class. See that class for
  a description.
*******************************************************************/

class AbsorberVmr : virtual public StateVectorObserver,
	            public Observable<AbsorberVmr> {
public:
  virtual ~AbsorberVmr() {}
  virtual void add_observer(Observer<AbsorberVmr>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<AbsorberVmr>& Obs) 
  { remove_observer_do(Obs, *this);}
  
//-----------------------------------------------------------------------
/// Clone a AbsorberVmr object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<AbsorberVmr>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" AbsorberVmr object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AbsorberVmr> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a pressure to use. The intent is that
/// the pressure has been cloned from the original pressure (although
/// this class has no way to verify this). This allows sets of objects
/// to be cloned using a common Pressure clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;

//-----------------------------------------------------------------------
/// This indicates the name of this particular Absorber. The naming
/// convention is free form but recommended to use the short form
/// often used by HITRAN
//-----------------------------------------------------------------------

  virtual std::string gas_name() const = 0;

//-----------------------------------------------------------------------
/// This returns the volume mixing ratio at the given pressure level.
/// This is dimensionless, and the pressure is in Pascals
//-----------------------------------------------------------------------

  virtual AutoDerivative<double> 
  volume_mixing_ratio(const AutoDerivative<double>& P) const = 0;

  virtual ArrayAd<double, 1> vmr_grid(const Pressure& P) const;

//-----------------------------------------------------------------------
/// Lua doesn't like ArrayAd, so have a version of vmr_grid that just
/// returns the value.
//-----------------------------------------------------------------------

  blitz::Array<double,1> vmr_grid_value(const Pressure& P) const
  { return vmr_grid(P).value(); }
  
//-----------------------------------------------------------------------
/// Indicate what portion of the state vector is used to calculate the
/// VMR. 
//-----------------------------------------------------------------------

  virtual blitz::Array<bool, 1> state_used() const = 0;

};
}
#endif
