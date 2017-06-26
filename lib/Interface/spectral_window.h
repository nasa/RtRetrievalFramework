#ifndef SPECTRAL_WINDOW_H
#define SPECTRAL_WINDOW_H
#include "printable.h"
#include <blitz/array.h>
#include <vector>
#include "spectral_bound.h"
#include "double_with_unit.h"
#include "spectrum.h"

namespace FullPhysics {
/****************************************************************//**
  This class represents a the spectral window.

  The definition of a spectral window is purposely fuzzy, we want to
  support things like excluding certain wavenumbers. So the interface
  simple takes a list of potential wavenumbers (e.g., the wavenumbers
  measured to the GOSAT oxygen A spectrometer) and returns the list of
  values that fall within the window. For example, if the window is
  just a wavenumber range, then all the wavenumbers that fall within
  that range are returned.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class SpectralWindow : public Printable<SpectralWindow> {
public:
  virtual ~SpectralWindow() {}

  SpectralDomain apply(const SpectralDomain& Grid, int Spec_index) const;
  Spectrum apply(const Spectrum& Spec, int Spec_index) const;

//-----------------------------------------------------------------------
/// Given a list of wavenumbers, this returns the indices that fall
/// within the window.
//-----------------------------------------------------------------------

  virtual std::vector<int> grid_indexes(const SpectralDomain& Grid, 
     int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Number of spectrometers.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const  = 0;

//-----------------------------------------------------------------------
/// Bounds of spectral window.
//-----------------------------------------------------------------------

  virtual SpectralBound spectral_bound() const = 0;

  virtual void print(std::ostream& Os) const { Os << "SpectralWindow";}
};
}
#endif
