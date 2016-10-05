#ifndef MODEL_MEASURE_OCO_H
#define MODEL_MEASURE_OCO_H
#include <model_measure.h>
#include <forward_model.h>

#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/******************************************************************
  This class implements what is common to 
    - "maximum a posteriori" for OCO, and
    - "maximum likelihood" for OCO
*******************************************************************/
class ModelMeasureOCO : 
  virtual public ModelMeasure {

public:

  virtual ~ModelMeasureOCO() {}

  virtual void model_eval();

  virtual void jacobian_eval();

  virtual void model_jacobian_eval();

  virtual int expected_parameter_size() const;


//-----------------------------------------------------------------------
/// Sets the problem at a new point in the parameter space.
/// 
/// \param x Input value
//-----------------------------------------------------------------------

  virtual void parameters(const blitz::Array<double, 1>& x);


//-----------------------------------------------------------------------
/// Just returns the current values of parameters.
/// This method is redefined here (see the root base
/// class) because of a compiler bug; otherwise, there
/// should be no need for its redefinition.
/// 
/// \return Current parameter values
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> parameters() const
  { return ModelMeasure::parameters(); }


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "ModelMeasureOCO"; }


protected:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  ModelMeasureOCO(const boost::shared_ptr<ForwardModel>& fm);

  ModelMeasureOCO() {}

  void radiance_from_fm(bool Skip_jacobian=false);

  //  TEMPORARY
  //
  // Should go away after we end support for 
  // fixed pressure level grid.
  virtual void vanishing_params_update();

  boost::shared_ptr<ForwardModel> FM;
  Unit meas_units;

};
}
#endif
