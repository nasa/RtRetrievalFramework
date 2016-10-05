#ifndef AEROSOL_EXTINCTION_H
#define AEROSOL_EXTINCTION_H
#include "state_vector.h"
#include "auto_derivative.h"
#include "pressure.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the aerosol extinction on each
  level.

  Other objects may depend on the AerosolExtinction, and should be updated
  when the AerosolExtinction is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the AerosolExtinction is updated.

  When implementing a new class, you almost always will want to derive
  from AerosolExtinctionImpBase rather than from this class. See that
  class for a description.
*******************************************************************/
class AerosolExtinction : virtual public StateVectorObserver,
			  public Observable<AerosolExtinction> {
public:
  virtual ~AerosolExtinction() {}
  virtual void add_observer(Observer<AerosolExtinction>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<AerosolExtinction>& Obs) 
  { remove_observer_do(Obs, *this);}
  
//-----------------------------------------------------------------------
/// Clone a AerosolExtinction object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<AerosolExtinction>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" AerosolExtinction object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AerosolExtinction> clone() const = 0;

//-----------------------------------------------------------------------
/// This version of clone takes a pressure to use. The intent is that
/// the pressure has been cloned from the original pressure (although
/// this class has no way to verify this). This allows sets of objects
/// to be cloned using a common Pressure clone, e.g. Atmosphere.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<AerosolExtinction> 
  clone(const boost::shared_ptr<Pressure>& Press) const = 0;

//-----------------------------------------------------------------------
/// Extinction for given layer.
//-----------------------------------------------------------------------

  virtual AutoDerivative<double> extinction_for_layer(int i) const = 0;

//-----------------------------------------------------------------------
/// Name of aerosol. 
//-----------------------------------------------------------------------
  virtual std::string aerosol_name() const = 0;
};
}
#endif
