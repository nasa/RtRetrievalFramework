#ifndef DISPERSION_H
#define DISPERSION_H
#include "state_vector.h"
#include "observer.h"
#include "spectral_domain.h"

namespace FullPhysics {
/****************************************************************//**
  This class calculates the wavenumber for each pixel in a single band
  of an Instrument.
*******************************************************************/
class Dispersion: virtual public StateVectorObserver,
		  public Observable<Dispersion> {
public:
  virtual ~Dispersion() {}
  virtual void add_observer(Observer<Dispersion>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Dispersion>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Clone an Dispersion object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Dispersion>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Dispersion object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Dispersion> clone() const = 0;

//-----------------------------------------------------------------------
/// Returns as list of grid points for each instrument pixel, and the
/// gradient of the points wrt the state vector. This is for the
/// full instrument pixels, i.e., any windowing etc. happens in later
/// processing.
//-----------------------------------------------------------------------

  virtual SpectralDomain pixel_grid() const = 0;
};
}
#endif
