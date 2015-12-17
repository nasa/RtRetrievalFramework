#ifndef STOKES_COEFFICIENT_H
#define STOKES_COEFFICIENT_H
#include "state_vector.h"
#include "observer.h"
#include "array_ad.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the stokes coefficient portion of the state.
*******************************************************************/

class StokesCoefficient : virtual public StateVectorObserver, 
			  public Observable<StokesCoefficient> {
public:
  virtual ~StokesCoefficient() {}
  virtual void add_observer(Observer<StokesCoefficient>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<StokesCoefficient>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Return Stokes coefficients used to go from Stokes vector to scalar
/// reflectance. This is number_spectrometer() x 4, and is unit less.
///
/// Note that is simple a matter of convenience that we have "4"
/// rather than just number_stokes(). This happens to be how the
/// stokes coefficients are given the Level 1 file. We only actually
/// use the first number_stokes() coefficients.
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> stokes_coefficient() const = 0;

//-----------------------------------------------------------------------
/// Clone a StokesCoefficient object. Note that the cloned version will *not*
/// be attached to a StateVector or Observer<StokesCoefficient>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" StokesCoefficient object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<StokesCoefficient> clone() const = 0;
};
}
#endif
