#ifndef SPURR_RT_H
#define SPURR_RT_H

#include "spurr_driver.h"
#include "radiative_transfer_single_wn.h"
#include "rt_atmosphere.h"
#include <boost/noncopyable.hpp>

namespace FullPhysics {

/****************************************************************//**
  Abstract Interface for Rt classes based on Spurr driver 
  implementations
 *******************************************************************/

class SpurrRt : public RadiativeTransferSingleWn,
		public Observer<RtAtmosphere>,
		public boost::noncopyable {
public:

  SpurrRt(const boost::shared_ptr<RtAtmosphere>& Atm,
	  const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
	  const blitz::Array<double, 1>& Sza, 
	  const blitz::Array<double, 1>& Zen, 
	  const blitz::Array<double, 1>& Azm);
  
  //-----------------------------------------------------------------------
  /// For performance, we cache some data as we calculate it. This
  /// becomes stale when the Atmosphere is changed, so we observe atm
  /// and mark the cache when it changes. 
  //-----------------------------------------------------------------------
  void notify_update(const RtAtmosphere& atm) { alt_spec_index_cache = -1; }

  /// Number of stokes in returned stokes values
  /// Note that LIDORT will only ever calculate the first stoke index for I,
  virtual int number_stokes() const { return stokes_coef->stokes_coefficient().cols(); }

  /// Number of quadtature streams in the cosine half space
  virtual int number_stream() const = 0;

  /// Number of moments for scattering matrix.
  virtual int number_moment() const = 0;

  /// Integer representing the surface type using the LIDORT indexing nomenclature
  virtual int surface_type() const { return surface_type_int; }

  virtual void print(std::ostream& Os, bool Short_form = false) const;

  virtual blitz::Array<double, 1> stokes_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;
  virtual ArrayAd<double, 1> stokes_and_jacobian_single_wn(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const;

protected:

  int surface_type_int;

  blitz::Array<double, 1> sza, zen, azm;

  boost::shared_ptr<SpurrRtDriver> rt_driver_;

  // Last index we updates the altitude/geometry for.
  mutable int alt_spec_index_cache, geo_spec_index_cache;
  virtual void update_altitude(int spec_index) const;
  virtual void update_geometry(int spec_index) const;
};
}
#endif

