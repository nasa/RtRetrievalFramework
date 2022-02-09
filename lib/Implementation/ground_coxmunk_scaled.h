#ifndef GROUND_COXMUNK_SCALED_H
#define GROUND_COXMUNK_SCALED_H

#include "ground_brdf_weight.h"
#include "ground_coxmunk.h"
#include "auto_derivative.h"
#include "sub_state_vector_proxy.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Coxmunk ground type times a polynomial
  scale.
*******************************************************************/

class GroundCoxmunkScaled: public Ground, public SubStateVectorProxy {
public:
  GroundCoxmunkScaled(const boost::shared_ptr<GroundCoxmunk>& Coxmunk, const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const boost::shared_ptr<GroundCoxmunk> coxmunk() const { return coxmunk_; }

  virtual const boost::shared_ptr<GroundBrdfWeight> brdf_weight() const { return brdf_weight_; }

  virtual boost::shared_ptr<Ground> clone() const;

  virtual void print(std::ostream& Os) const;

  virtual std::string desc() const { return "GroundCoxmunkScaled"; }

private:
  boost::shared_ptr<GroundCoxmunk> coxmunk_;
  boost::shared_ptr<GroundBrdfWeight> brdf_weight_;
};
}
#endif
