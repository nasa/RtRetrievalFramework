#ifndef BARD_ML_PROBLEM
#define BARD_ML_PROBLEM
#include <max_likelihood.h>
#include <model_measure_bard.h>


namespace FullPhysics {

class BardMLProblem : public MaxLikelihood, public ModelMeasureBard {
public:
  BardMLProblem(const blitz::Array<double, 1>& measurement, 
                const blitz::Array<double, 1>& measurement_error_cov);
  virtual ~BardMLProblem() {}
  virtual void print(std::ostream& Os) const
  { Os << "BardMLProblem"; }
};
}

#endif
