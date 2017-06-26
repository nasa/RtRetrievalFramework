#ifndef FORWARD_MODEL_SPECTRAL_GRID_H
#define FORWARD_MODEL_SPECTRAL_GRID_H
#include "printable.h"
#include "spectrum.h"
#include "instrument.h"
#include "spectral_window.h"
#include "spectrum_sampling.h"
#include "spectral_domain.h"

namespace FullPhysics {

/****************************************************************//**
  This is the Forward Model spectral grid. This is in a separate class
  because this is a bit complicated. We have 3 grids to worry about

  1. The low resolution grid, which is the ultimate output of the 
     ForwardModel.
  2. The high resolution grid, which is where we calculate the RT
     on. This is the spectrum before we have convolved it with Ils.
     Depending on the options used, this grid might be nonuniform.
  3. The high resolution grid we interpolate to. If we are *not*
     doing nonuniform sampling, then this is the same as the high
     resolution grid.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class ForwardModelSpectralGrid : public Printable<ForwardModelSpectralGrid> {
public:
  ForwardModelSpectralGrid(
   const boost::shared_ptr<Instrument>& Inst,
   const boost::shared_ptr<SpectralWindow>& Spectral_window,
   const boost::shared_ptr<SpectrumSampling>& Spectrum_sampling) 
      : inst(Inst), spectral_window(Spectral_window), spectrum_sampling(Spectrum_sampling) {};

  ForwardModelSpectralGrid() {}
  virtual ~ForwardModelSpectralGrid() {}
  virtual void print(std::ostream& Os) const {Os << "ForwardModelSpectralGrid";}

  //-----------------------------------------------------------------------
  /// Number of spectrometer.
  //-----------------------------------------------------------------------

  int number_spectrometer() const { return spectral_window->number_spectrometer(); }

  //-----------------------------------------------------------------------
  /// The low resolution grid.
  //-----------------------------------------------------------------------

  const SpectralDomain low_resolution_grid(int Spec_index) const;

  //-----------------------------------------------------------------------
  /// The high resolution grid, possibly nonuniform
  //-----------------------------------------------------------------------

  const SpectralDomain high_resolution_grid(int Spec_index) const;

  //-----------------------------------------------------------------------
  /// The high resolution grid, interpolated to be uniform.
  //-----------------------------------------------------------------------

  const SpectralDomain high_resolution_interpolated_grid(int Spec_index) const;

  Spectrum interpolate_spectrum(const Spectrum& Spec_in, int Spec_index) const;

  //-----------------------------------------------------------------------
  /// Pixel indexes to use for low resolution grid.
  //-----------------------------------------------------------------------
  const std::vector<int> pixel_list(int Spec_index) const;

private:
  boost::shared_ptr<Instrument> inst;
  boost::shared_ptr<SpectralWindow> spectral_window;
  boost::shared_ptr<SpectrumSampling> spectrum_sampling;
};
}
#endif
