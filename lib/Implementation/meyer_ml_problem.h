#ifndef MEYER_ML_PROBLEM
#define MEYER_ML_PROBLEM
#include <max_likelihood.h>
#include <model_measure_meyer.h>


namespace FullPhysics {

class MeyerMLProblem : public MaxLikelihood, public ModelMeasureMeyer {
public:
  MeyerMLProblem(const blitz::Array<double, 1>& measurement, 
                 const blitz::Array<double, 1>& measurement_error_cov);
  virtual ~MeyerMLProblem() {}
  virtual void print(std::ostream& Os) const
  { Os << "MeyerMLProblem"; }
};
}

#endif
