#ifndef CHAPMAN_BOA_RT_H
#define CHAPMAN_BOA_RT_H

#include "radiative_transfer.h"
#include "atmosphere_oco.h"
#include "chapman_boa.h"
#include "spectral_bound.h"

namespace FullPhysics {

/****************************************************************//**
*******************************************************************/

class ChapmanBoaRT : public RadiativeTransfer,
		     public Observer<RtAtmosphere>,
                     public boost::noncopyable  {
public:
  ChapmanBoaRT(const boost::shared_ptr<AtmosphereOco>& Atm,
  	       const blitz::Array<double, 1>& Sza);

  ChapmanBoaRT(const boost::shared_ptr<AtmosphereOco>& Atm,
	       const blitz::Array<double, 1>& Sza, 
	       const SpectralBound& Spec_bound);

  virtual ~ChapmanBoaRT() {}

  //-----------------------------------------------------------------------
  /// Regenerate chapman factors when Atmosphere changes
  //-----------------------------------------------------------------------
  void notify_update(const RtAtmosphere& updated_atm) { chapman_cache_stale = true; }

  virtual int number_stokes() const { return 1; }

  virtual int number_spectrometer() const
  { return sza.extent(blitz::firstDim);}

  //-----------------------------------------------------------------------
  /// Pointer to the Atmosphere class we are using.
  //-----------------------------------------------------------------------

  const boost::shared_ptr<AtmosphereOco>& atmosphere_ptr() const {return atm;}

  // See description in base class
  virtual Spectrum reflectance
  (const SpectralDomain& Spec_domain, int Spec_index, 
   bool Skip_jacobian = false) const;
  virtual blitz::Array<double, 2> stokes(const SpectralDomain& Spec_domain, int Spec_index) const;
  virtual ArrayAd<double, 2> stokes_and_jacobian(const SpectralDomain& Spec_domain, int Spec_index) const;
 
  virtual void print(std::ostream& Os, bool Short_form = false) const;

private:

  /// Computes Chapman factors on demand
  void compute_chapman_factors(const int Spec_idx) const;

  /// Handles flagging whether or not chapman factors need to be regenerated
  mutable blitz::Array<bool, 1> chapman_cache_stale;

  /// Chapman BOA geometry class
  mutable std::vector<boost::shared_ptr<ChapmanBOA> > chapman_boa;

  /// Window range for use when picking refraction method
  SpectralBound spec_bound;

  /// Atmosphere object we are using.
  boost::shared_ptr<AtmosphereOco> atm;

  /// Solar zenith angles per spectrometer
  blitz::Array<double, 1> sza;

};
}

#endif
