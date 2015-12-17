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

  //-----------------------------------------------------------------------
  /// Create a new SpectrallyResolvedNoise modifying an existing noise
  /// model. The created object should be used in place of the base model.
  /// The noise coefficents for each band must be set individually since
  /// they might vary in size
  //-----------------------------------------------------------------------
  SpectrallyResolvedNoise(boost::shared_ptr<NoiseModel> Base_model) : base_model_(Base_model) {};

  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;
  virtual void set_noise_coefficients(int Spec_index, const blitz::Array<double, 1> Noise_coeff); 

  virtual void print(std::ostream& Os) const;

private:
  boost::shared_ptr<NoiseModel> base_model_;
  mutable std::vector<blitz::Array<double, 1> > band_noise_coeffs_;
};
}
#endif


