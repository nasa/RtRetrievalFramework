#ifndef MODEL_MEASURE_MEYER_H
#define MODEL_MEASURE_MEYER_H
#include <model_measure.h>


namespace FullPhysics {
/******************************************************************
  This class implements what is common to 
    - "maximum a posteriori" for Meyer, and
    - "maximum likelihood" for Meyer
*******************************************************************/
class ModelMeasureMeyer : 
  virtual public ModelMeasure {

public:

  virtual ~ModelMeasureMeyer() {}

  virtual void model_eval();

  virtual void jacobian_eval();

  virtual void model_jacobian_eval();

  virtual int expected_parameter_size() const { return 3; }

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelMeasureMeyer"; }


protected:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  ModelMeasureMeyer() {}

};
}
#endif
