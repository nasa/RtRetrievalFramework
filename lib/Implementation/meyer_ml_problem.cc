#include <meyer_ml_problem.h>


using namespace FullPhysics;
using namespace blitz;

MeyerMLProblem::MeyerMLProblem(
    const blitz::Array<double, 1>& measurement, 
    const blitz::Array<double, 1>& measurement_error_cov)
  : ModelMeasure(measurement, measurement_error_cov),
    ModelMeasureMeyer()
{
}
