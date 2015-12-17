#ifndef INSTRUMENT_CORRECTION_H
#define INSTRUMENT_CORRECTION_H
#include "state_vector.h"
#include "spectral_range.h"
#include "spectral_domain.h"

namespace FullPhysics {
/****************************************************************//**
  This class models an Instrument correction. This is used by
  IlsConvolution, and it applies zero or more corrections that allow
  the radiance results to be modified. Examples are a zero level
  offset correction, and a continuum correction.
*******************************************************************/

class InstrumentCorrection : virtual public StateVectorObserver,
			     public Observable<InstrumentCorrection> {
public:
  virtual ~InstrumentCorrection() {}

  virtual void add_observer(Observer<InstrumentCorrection>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<InstrumentCorrection>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Apply correction to radiance values, in place. If Radiance
/// includes a Jacobian, then we include the Jacobian
/// calculation. Otherwise we don't include the Jacobian in the
/// calculation. 
///
/// \param Pixel_grid - The grid point of each pixel. We only
///    use a subset of these points, but the full list is passed
///    in for use by the class. 
/// \param Pixel_list - List of pixels that actually appear in
///    Radiance, in the order that they appear.
/// \param Radiance - Radiance values, that will be corrected in
///    place. 
//-----------------------------------------------------------------------

  virtual void apply_correction
  (const SpectralDomain& Pixel_grid,
   const std::vector<int>& Pixel_list,
   SpectralRange& Radiance) const = 0;

//-----------------------------------------------------------------------
/// Clone an InstrumentCorrection object. Note that the cloned version
/// will *not* be attached to and StateVector or
/// Observer<InstrumentCorrection>, although you can of course attach
/// them after receiving the cloned object. 
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" InstrumentCorrection object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<InstrumentCorrection> clone() const = 0;
};
}
#endif
