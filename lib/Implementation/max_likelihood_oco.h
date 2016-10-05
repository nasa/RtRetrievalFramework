#ifndef MAX_LIKELIHOOD_OCO_H
#define MAX_LIKELIHOOD_OCO_H
#include <max_likelihood.h>
#include <model_measure_oco.h>


namespace FullPhysics {
/******************************************************************
  This class implements "maximum likelihood" for OCO.
*******************************************************************/
class MaxLikelihoodOCO : 
  public MaxLikelihood, public ModelMeasureOCO {

public:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  MaxLikelihoodOCO(const boost::shared_ptr<ForwardModel>& fm);

  virtual ~MaxLikelihoodOCO() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxLikelihoodOCO"; }

};
}
#endif
