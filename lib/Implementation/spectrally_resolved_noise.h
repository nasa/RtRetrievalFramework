#ifndef SPECTRALLY_RESOLVED_NOISE_H
#define SPECTRALLY_RESOLVED_NOISE_H

#include "noise_model.h"
#include <boost/shared_ptr.hpp>
#include <vector>

namespace FullPhysics {

/****************************************************************//**
  Adds a spectrally varying modification to an existing noise 
  model. Noise coefficents scale the parent noise
  model's uncertainty at certain grid points.
*******************************************************************/

class SpectrallyResolvedNoise: public NoiseModel {
public:

  SpectrallyResolvedNoise(boost::shared_ptr<NoiseModel> Base_model);

  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
  virtual void set_full_noise_scaling(int Spec_index, const blitz::Array<double, 1> Noise_scaling);
  virtual void set_single_noise_scaling(int Spec_index, double Noise_scaling);

  virtual void print(std::ostream& Os) const;

private:
  boost::shared_ptr<NoiseModel> base_model_;
  mutable std::vector<bool> band_single_value_;
  mutable std::vector<blitz::Array<double, 1> > band_noise_coeffs_;
};
}
#endif


