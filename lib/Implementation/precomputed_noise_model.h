#ifndef PRECOMPUTED_NOISE_MODEL_H
#define PRECOMPUTED_NOISE_MODEL_H
#include "heritage_file.h"
#include "noise_model.h"

namespace FullPhysics {
/****************************************************************//**
  This class creates a noise model which simply passes back values
  read from a file.
*******************************************************************/
class PrecomputedNoiseModel: public NoiseModel {
public:

  PrecomputedNoiseModel(const HeritageFile& Noise_file);

  virtual blitz::Array<double, 1> uncertainty(int Spec_index, const blitz::Array<double, 1>& Radiance) const;

  virtual void print(std::ostream& Os) const;

private:

  std::vector< blitz::Array<double, 1> > file_noise_values;
  
};
}
#endif


