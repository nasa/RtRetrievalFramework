#ifndef NOISE_MODEL_H
#define NOISE_MODEL_H
#include <blitz/array.h>
#include <stdint.h>
#include "printable.h"

namespace FullPhysics {
/****************************************************************//**
  Interface for calculating noise/uncertainty values from 
  radiance data given some internal representation of the noise 
  model
*******************************************************************/
class NoiseModel: public Printable<NoiseModel> {
public:
  virtual ~NoiseModel() {}

  //-----------------------------------------------------------------------
  /// Uncertainty on radiance, for given spectral band.
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> uncertainty(int Spec_index, 
      const blitz::Array<double, 1>& Radiance) const = 0;

  //-----------------------------------------------------------------------
  /// Print description of object.
  //-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "NoiseModel";}

};
}
#endif
