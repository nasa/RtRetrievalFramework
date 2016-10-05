#ifndef MAX_LIKELIHOOD_H
#define MAX_LIKELIHOOD_H
#include <model_measure.h>


namespace FullPhysics {

//-----------------------------------------------------------------------
/// \brief The base class for maximum likelihood estimation.
///
/// This class is the base class for all classes that use 
/// maximum likelihood estimation method to implement the problem
/// of estimating the parameters of a statistical model.
//-----------------------------------------------------------------------

class MaxLikelihood : 
    virtual public ModelMeasure {

public:


//-----------------------------------------------------------------------
/// \brief Default Constructor
//-----------------------------------------------------------------------

  MaxLikelihood() {}


  virtual ~MaxLikelihood() {}


//-----------------------------------------------------------------------
/// \brief Prints description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "MaxLikelihood"; }


protected:


};
}
#endif
