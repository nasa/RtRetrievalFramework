#ifndef FORWARD_MODEL_SPECTRAL_GRID_H
#define FORWARD_MODEL_SPECTRAL_GRID_H
#include "printable.h"
#include "spectrum.h"
#include <vector>

namespace FullPhysics {
  class Instrument;
  class SpectralWindow;
  class SpectrumSampling;
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
   const Instrument& Inst,
   const SpectralWindow& Spectral_window,
   const SpectrumSampling& Spectrum_sampling);
  ForwardModelSpectralGrid() {}
  virtual ~ForwardModelSpectralGrid() {}
  virtual void print(std::ostream& Os) const {Os << "ForwardModelSpectralGrid";}

//-----------------------------------------------------------------------
/// Number of spectrometer.
//-----------------------------------------------------------------------

  int number_spectrometer() const { return (int) lgrid.size(); }

//-----------------------------------------------------------------------
/// The low resolution grid.
//-----------------------------------------------------------------------

  const SpectralDomain& low_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return lgrid[Spec_index];
  }

//-----------------------------------------------------------------------
/// The high resolution grid, possibly nonuniform
//-----------------------------------------------------------------------

  const SpectralDomain& high_resolution_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return hgrid[Spec_index]; 
  }

//-----------------------------------------------------------------------
/// The high resolution grid, interpolated to be uniform.
//-----------------------------------------------------------------------

  const SpectralDomain& high_resolution_interpolated_grid(int Spec_index) const
  { 
    range_check(Spec_index, 0, number_spectrometer());
    return hgrid_inter[Spec_index]; 
  }
  Spectrum interpolate_spectrum(const Spectrum& Spec_in, int Spec_index) const;

//-----------------------------------------------------------------------
/// Pixel indexes to use for low resolution grid.
//-----------------------------------------------------------------------
  const std::vector<int>& pixel_list(int Spec_index) const
  {
    range_check(Spec_index, 0, number_spectrometer());
    return plist[Spec_index];
  }
private:
  std::vector<SpectralDomain> lgrid, hgrid, hgrid_inter;
  std::vector<std::vector<int> > plist;
  std::vector<bool> need_interpolation;
};
}
#endif
