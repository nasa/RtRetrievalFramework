#ifndef RADIANT_DRIVER_H
#define RADIANT_DRIVER_H
#include "atmosphere_oco.h"
#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"
#include <boost/noncopyable.hpp>
#include <blitz/array.h>

using namespace FullPhysics;
using namespace blitz;

namespace FullPhysics {
/****************************************************************//**
  This class drives the RADIANT Radiative Transfer code.
*******************************************************************/
class RadiantDriver : public RadiativeTransferSingleWn,
                      public Observer<RtAtmosphere>,
                      public boost::noncopyable {
public:
  RadiantDriver(const boost::shared_ptr<RtAtmosphere>& Atm,
		const blitz::Array<double, 1>& Sza);
  virtual ~RadiantDriver();
  
  //-----------------------------------------------------------------------
  /// For performance, we cache some data as we calculate it. This
  /// becomes stale when the Atmosphere is changed, so we observe atm
  /// and mark the cache when it changes. 
  //-----------------------------------------------------------------------
  
  void notify_update(const RtAtmosphere& atm) { alt_spec_index_cache = -1; }
  
  //-----------------------------------------------------------------------
  /// Number of moments for scattering matrix.
  //-----------------------------------------------------------------------
  
  virtual int number_stream() const {return 0;}

  virtual void print(std::ostream& Os) const;

  virtual blitz::Array<double, 1> stokes_single_wn
    (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn
    (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;

private:
  static constexpr double SOLID_ANGLE = 6.79929414e-5;
  virtual int number_stokes() const {return 1;}

  blitz::Array<double, 1> sza;

  mutable blitz::Array<double, 1> tau_in;
  mutable blitz::Array<double, 2> l_tau_in;

  ArrayAd<double, 1> stokes_and_maybe_jacobian
    (double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;

  // Last index we updates the altitude for.
  mutable int alt_spec_index_cache;

};
}
#endif

