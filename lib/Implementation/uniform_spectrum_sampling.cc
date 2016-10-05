#include "uniform_spectrum_sampling.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor. This creates a fixed spacing grid. wn_end may be
/// adjusted slightly if it isn't exactly wn_start + wn_step * i for
/// some integer i.
//-----------------------------------------------------------------------

UniformSpectrumSampling::UniformSpectrumSampling(double wn_start, 
						 double wn_end, double wn_step)
: SpectrumSampling(1)
{
  // 0.05 is a bit arbitrary, it is just "a little more than wn_end"
  int nspec = (int) floor((wn_end + 0.05 * wn_step - wn_start) / wn_step) + 1;
  Array<double, 1> wn(nspec);
  for(int i = 0; i < nspec; ++i)
    wn(i) = wn_start + i * wn_step;
  spec_domain.push_back(SpectralDomain(wn));
}

//-----------------------------------------------------------------------
/// Constructor for 3 spectrum (this is a useful case, because it
/// matches GOSAT)
//-----------------------------------------------------------------------

UniformSpectrumSampling::UniformSpectrumSampling(
                          double wn_start1, double wn_end1, double wn_step1,
			  double wn_start2, double wn_end2, double wn_step2,
			  double wn_start3, double wn_end3, double wn_step3)
: SpectrumSampling(3)
{
  int nspec1 = (int) floor((wn_end1 + 0.05 * wn_step1 - wn_start1) / wn_step1) + 1;
  int nspec2 = (int) floor((wn_end2 + 0.05 * wn_step2 - wn_start2) / wn_step2) + 1;
  int nspec3 = (int) floor((wn_end3 + 0.05 * wn_step3 - wn_start3) / wn_step3) + 1;
  Array<double, 1> wn(nspec1);
  for(int i = 0; i < nspec1; ++i)
    wn(i) = wn_start1 + i * wn_step1;
  spec_domain.push_back(SpectralDomain(wn.copy()));
  wn.resize(nspec2);
  for(int i = 0; i < nspec2; ++i)
    wn(i) = wn_start2 + i * wn_step2;
  spec_domain.push_back(SpectralDomain(wn.copy()));
  wn.resize(nspec3);
  for(int i = 0; i < nspec3; ++i)
    wn(i) = wn_start3 + i * wn_step3;
  spec_domain.push_back(SpectralDomain(wn.copy()));
}

//-----------------------------------------------------------------------
/// Print to stream.
//-----------------------------------------------------------------------

void UniformSpectrumSampling::print(std::ostream& Os) const 
{ 
  Os << "UniformSpectrumSampling";
  for(int i = 0; i < number_spectrometer(); ++i)
    Os << "  " << i + 1 << ":\n"
       << "     wn_start: " << spec_domain[i].data()(0) << "\n"
       << "     wn_end:   " << spec_domain[i].data()(spec_domain[i].data().rows() - 1) << "\n"
       << "     wn_space: " << spec_domain[i].data()(1) - spec_domain[i].data()(0) << "\n";
}
