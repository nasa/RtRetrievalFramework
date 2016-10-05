#ifndef SPECTRUM_SAMPLING_FIXED_SPACING_H
#define SPECTRUM_SAMPLING_FIXED_SPACING_H
#include "spectrum_sampling.h"

namespace FullPhysics {
/****************************************************************//**
  This generates a spectrum sampling that covers all the high
  resolution points needed to create the spectral domain covered by
  the given Instrument, subject to the given low resolution grid. For
  each range in the spectrum, we produce equally spaced points.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/

class SpectrumSamplingFixedSpacing : public SpectrumSampling {
public:
  SpectrumSamplingFixedSpacing(const ArrayWithUnit<double, 1>& Spec_spacing)
  : spec_spacing(Spec_spacing) { }

  virtual ~SpectrumSamplingFixedSpacing() {}
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const;
  virtual void print(std::ostream& Os) const 
  { Os << "SpectrumSamplingFixedSpacing\n"
       << "  Spacing: " << spec_spacing << "\n";}
private:
  ArrayWithUnit<double, 1> spec_spacing;
};
}
#endif
