#ifndef UNIFORM_SPECTRUM_SAMPLING_H
#define UNIFORM_SPECTRUM_SAMPLING_H
#include "spectrum_sampling.h"
#include "fp_exception.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This is a simple SpectrumSampling that is just a uniform
  sampling. This is useful for doing unit testing.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class UniformSpectrumSampling : public SpectrumSampling {
public:
  UniformSpectrumSampling(double wn_start, double wn_end, double wn_step);
  UniformSpectrumSampling(double wn_start1, double wn_end1, double wn_step1,
			  double wn_start2, double wn_end2, double wn_step2,
			  double wn_start3, double wn_end3, double wn_step3);

  virtual ~UniformSpectrumSampling() {}

//-----------------------------------------------------------------------
/// Wave numbers to use for the given spectrometer. 
//-----------------------------------------------------------------------

  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const
  { range_check(spec_index, 0, number_spectrometer()); 
    return spec_domain[spec_index];}
  virtual void print(std::ostream& Os) const;
private:
  std::vector<SpectralDomain> spec_domain;
};
}
#endif
