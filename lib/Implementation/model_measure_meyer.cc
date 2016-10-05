#include <model_measure_meyer.h>


using namespace FullPhysics;
using namespace blitz;


void ModelMeasureMeyer::model_eval()
{
  if(M.size() <= 0) {
    assert_parameter_set_correctly();
    M.resize(measurement_size());

    for(int i=0; i<measurement_size(); i++)
      M(i) = X(0) * exp( X(1)/((45.0 + 5.0*(i+1))+X(2)) );
  }
}



void ModelMeasureMeyer::jacobian_eval()
{
  if(K.size() <= 0) {
    assert_parameter_set_correctly();
    K.resize(measurement_size(), expected_parameter_size());

    for(int i=0; i<measurement_size(); i++) {
      double t = (45.0 + 5.0*(i+1));
      double expf = exp( X(1)/(t+X(2)) );
      K(i,0) = expf;  K(i,1) = X(0)/(t+X(2))*expf;  K(i,2) = -X(0)*X(1)/((t+X(2))*(t+X(2)))*expf;
    }
  }
}


void ModelMeasureMeyer::model_jacobian_eval()
{
  model_eval();
  jacobian_eval();
}
