#ifndef GROUND_COXMUNK_PLUS_LAMB_H
#define GROUND_COXMUNK_PLUS_LAMB_H

#include "ground_lambertian.h"
#include "ground_coxmunk.h"
#include "auto_derivative.h"
#include "sub_state_vector_proxy.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Coxmunk plus Lambertian ground type. 
*******************************************************************/

class GroundCoxmunkPlusLambertian: public Ground, public SubStateVectorProxy {
public:
  GroundCoxmunkPlusLambertian(const boost::shared_ptr<GroundCoxmunk>& Coxmunk, const boost::shared_ptr<GroundLambertian>& Lambertian);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const boost::shared_ptr<GroundCoxmunk> coxmunk() const { return coxmunk_; }

  virtual const boost::shared_ptr<GroundLambertian> lambertian() const { return lambertian_; }

  virtual boost::shared_ptr<Ground> clone() const;

  virtual void print(std::ostream& Os) const;

  virtual std::string desc() const { return "GroundCoxmunkPlusLambertian"; }
  
private:
  boost::shared_ptr<GroundCoxmunk> coxmunk_;
  boost::shared_ptr<GroundLambertian> lambertian_;
};
}
#endif
