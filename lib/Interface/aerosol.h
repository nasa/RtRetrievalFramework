#ifndef AEROSOL_H
#define AEROSOL_H
#include "state_vector.h"
#include "accumulated_timer.h"
#include "pressure.h"
#include "relative_humidity.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the aerosol portion of the state.

  Other objects may depend on the aerosol, and should be updated
  when the aerosol is updated. To facilitate that, this class in
  an Oberverable, and objects can add themselves as Observers to be
  notified when the aerosol is updated.

  I'm not really sure what the interface for this class should be.
  Right now it is used only by AtmosphereOco, and there is only one
  instance AerosolOptical, so the functions are what AtmosphereOco
  needs. But we may perhaps want to modify this in the future to be
  more general. 
*******************************************************************/

class Aerosol: public StateVectorObserver, public Observable<Aerosol> {
public:
  virtual ~Aerosol() {}
  // Used as a convenience to collect timing information to report in 
  // logging
  static AccumulatedTimer timer;

  virtual void add_observer(Observer<Aerosol>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Aerosol>& Obs)
  { remove_observer_do(Obs, *this);}

  virtual boost::shared_ptr<Aerosol> clone() const = 0;
  virtual boost::shared_ptr<Aerosol> 
  clone(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<RelativeHumidity>& Rh) const = 0;

//-----------------------------------------------------------------------
/// This calculates the portion of the phase function moments that
/// come from the aerosol.
/// \param wn The wave number.
/// \param frac_aer This is number_active_layer() x number_particle()
/// \param nummom Number of moments to fill in
/// \param numscat Number of scatters to fill in
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 3> pf_mom(double wn, 
         const ArrayAd<double, 2>& frac_aer,
         int nummom = -1, int numscat = -1) const = 0;

//-----------------------------------------------------------------------
/// Number of aerosol particles
//-----------------------------------------------------------------------

  virtual int number_particle() const = 0;

//-----------------------------------------------------------------------
/// This gives the optical depth for each layer, for the given wave
/// number. Note this only includes the aerosol portion of this,
/// Atmosphere class combines this with Absorbers and rayleigh
/// scattering.
///
/// This calculates the derivatives with respect to the state vector.
///
/// This has size of number_active_layer() x number_particle().
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 2> optical_depth_each_layer(double wn) 
    const = 0;

//-----------------------------------------------------------------------
/// This gives the single scatter albedo for each layer, for the given wave
/// number, for the given particle. Note this only includes the
/// aerosol portion of this, 
/// Atmosphere class combines this with Rayleigh scattering.
///
/// We take in the optical depth of each layer. This is just what is
/// returned by optical_depth_each_layer(), we take this in because
/// we can change what the derivative of optical_depth_each_layer is
/// respect to, e.g. in AtmosphereOco we use taua_i.
///
/// This calculates the derivative with respect to whatever variables
/// Od is relative to.
///
/// This has size of number_active_layer()
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> 
  ssa_each_layer(double wn, int particle_index,
		 const ArrayAd<double, 1>& Od) const = 0;
};
}
#endif
