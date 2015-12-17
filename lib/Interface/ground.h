#ifndef GROUND_H
#define GROUND_H
#include "state_vector.h"
#include "observer.h"
#include "auto_derivative.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the ground portion of the state.

  Other objects may depend on the ground, and should be updated
  when the ground is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the ground is updated.

  This class is unfortunately a bit hard coded. The surface types are
  one of a set of enumerations. The surface parameters depend on
  exactly what the surface type is. These types and parameters map to
  types and parameters found in the LIDORT and LRAD code, so the hard
  coding is intrinsic. 
 *******************************************************************/
class Ground : virtual public StateVectorObserver,
               public Observable<Ground> {
public:
  virtual ~Ground() {}
  virtual void add_observer(Observer<Ground>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Ground>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Surface parmeters.
///
/// What exactly these parameters mean is determined by the surface
/// type, see the discussion in the comments before the Ground class.
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> surface_parameter
    (const double wn, const int spec_index) const = 0;

//-----------------------------------------------------------------------
/// Clone a Ground object. Note that the cloned version will *not*
/// be attached to a StateVector, although you can of course attach
/// them after receiving the cloned object. 
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Ground object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Ground> clone() const = 0;

};
}
#endif
