#ifndef SPECTRUM_EFFECT_H
#define SPECTRUM_EFFECT_H
#include "state_vector.h"
#include "spectrum.h"

namespace FullPhysics {
  class ForwardModelSpectralGrid;

/****************************************************************//**
  This class models models any effects that need to be applied
  to high resolution spectra after the radiative transfer model
  has finished its work.
*******************************************************************/

class SpectrumEffect : virtual public StateVectorObserver,
                               public Observable<SpectrumEffect> {
public:
  virtual ~SpectrumEffect() {}

  virtual void add_observer(Observer<SpectrumEffect>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<SpectrumEffect>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Apply correction to spectrum in place. 
/// We pass in the forward model grids used. A class can use this to 
/// optimize its calculation, see for example FluorescenceEffect.
//-----------------------------------------------------------------------

  virtual void apply_effect(Spectrum& Spec, 
	    const ForwardModelSpectralGrid& Forward_model_grid) const = 0;

//-----------------------------------------------------------------------
/// Clone a SpectrumEffect object. Note that the cloned version
/// will *not* be attached to and StateVector or
/// Observer<SpectrumEffect>, although you can of course attach
/// them after receiving the cloned object. 
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" SpectrumEffect object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<SpectrumEffect> clone() const = 0;

//-----------------------------------------------------------------------
/// Name of spectrum effect, for use when outputting effects of effect
//-----------------------------------------------------------------------

  virtual std::string name() const = 0;

};
}
#endif
