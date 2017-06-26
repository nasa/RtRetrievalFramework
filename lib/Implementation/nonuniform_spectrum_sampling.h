#ifndef NONUNIFORM_SPECTRUM_SAMPLING_H
#define NONUNIFORM_SPECTRUM_SAMPLING_H
#include "heritage_file.h"
#include "spectrum_sampling.h"
#include "fp_exception.h"
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This is a simple SpectrumSampling that is just a nonuniform
  sampling.

  We have an underlying Interpolated_sampling that we will 
  interpolate (by interpolate_spectrum).

  The input Grid might not overlap the Interpolated_sampling grid. 
  It might cover a larger area, or it might not exactly line up with
  the Interpolated_sampling. To reduce interpolation errors and also
  to avoid doing unnecessary calculations (e.g., the Grid is larger
  than Interpolated_sampling), we have spectrum_domain return *only* 
  a subset of the Interpolated_sampling. We take each Grid point and 
  find the closest Interpolated_sampling and return that. Note that 
  more than one Grid point might be closest to a single 
  Interpolated_sampling, which means that we may return a spectral
  domain smaller than the Interpolated_sampling.

  Note that there are a few closely related classes, with similar 
  sounding names. See \ref spectrum_doxygen for a description of each
  of these.
*******************************************************************/
class NonuniformSpectrumSampling : public SpectrumSampling {
public:
  NonuniformSpectrumSampling(const SpectralDomain& Grid,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);
  NonuniformSpectrumSampling(const SpectralDomain& Grid1,
     const SpectralDomain& Grid2,
     const SpectralDomain& Grid3,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);

  NonuniformSpectrumSampling(const HeritageFile& Grid_file,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);
  NonuniformSpectrumSampling(const HeritageFile& Grid_file1,
     const HeritageFile& Grid_file2,
     const HeritageFile& Grid_file3,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);

  NonuniformSpectrumSampling(const std::string& Grid_file,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);
  NonuniformSpectrumSampling(const std::string& Grid_file1,
     const std::string& Grid_file2,
     const std::string& Grid_file3,
     const boost::shared_ptr<SpectrumSampling>& Interpolated_sampling);


  virtual ~NonuniformSpectrumSampling() {}

  virtual SpectralDomain spectral_domain_interpolated(int Spec_index, 
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const
  { return interpolated_sampling->
      spectral_domain(Spec_index, Lowres_grid, Ils_half_width); }
  virtual SpectralDomain spectral_domain(int spec_index,
		 const SpectralDomain& Lowres_grid, 
		 const DoubleWithUnit& Ils_half_width) const;

  virtual bool need_interpolation(int Spec_index) const
  { return (spec_domain[Spec_index].data().rows() > 0); }
  virtual void print(std::ostream& Os) const;
private:
  SpectralDomain sort_sd(const SpectralDomain& In) const;
  SpectralDomain sort_sd(const HeritageFile& In) const;
  boost::shared_ptr<SpectrumSampling> interpolated_sampling;
  std::vector<SpectralDomain> spec_domain;
};
}
#endif
