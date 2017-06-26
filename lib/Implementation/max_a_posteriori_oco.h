#ifndef MAX_A_POSTERIORI_OCO_H
#define MAX_A_POSTERIORI_OCO_H
#include <max_a_posteriori.h>
#include <model_measure_oco.h>


namespace FullPhysics {
/******************************************************************
  This class implements "maximum a posteriori" for OCO.
*******************************************************************/
class MaxAPosterioriOCO : 
  public MaxAPosteriori, public ModelMeasureOCO {

public:

//-----------------------------------------------------------------------
/// Constructor
//-----------------------------------------------------------------------

  MaxAPosterioriOCO(const boost::shared_ptr<ForwardModel>& fm,
                    const blitz::Array<double, 1> a_priori_params,
                    const blitz::Array<double, 2> a_priori_cov);

  virtual ~MaxAPosterioriOCO() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxAPosterioriOCO"; }

protected:

  //  TEMPORARY
  //
  // Should go away after we end support for 
  // fixed pressure level grid.
  virtual void vanishing_params_update();

};
}
#endif
