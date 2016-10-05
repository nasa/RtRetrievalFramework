#ifndef MODEL_MEASURE_BARD_H
#define MODEL_MEASURE_BARD_H
#include <model_measure.h>


namespace FullPhysics {
/******************************************************************
  This class implements what is common to 
    - "maximum a posteriori" for Bard, and
    - "maximum likelihood" for Bard
*******************************************************************/
class ModelMeasureBard : 
  virtual public ModelMeasure {

public:

  virtual ~ModelMeasureBard() {}

  virtual void model_eval();

  virtual void jacobian_eval();

  virtual void model_jacobian_eval();

  virtual int expected_parameter_size() const { return 3; }

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelMeasureBard"; }


protected:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  ModelMeasureBard() {}

};
}
#endif
