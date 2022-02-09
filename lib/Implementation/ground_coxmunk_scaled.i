// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_coxmunk_scaled.h"
#include "sub_state_vector_proxy.h"
%}

%base_import(ground)
%base_import(sub_state_vector_proxy)
%import "ground_coxmunk.i"
%import "ground_brdf_weight.i"

%fp_shared_ptr(FullPhysics::GroundCoxmunkScaled);
namespace FullPhysics {
class GroundCoxmunkScaled: public Ground, public SubStateVectorProxy {
public:
  GroundCoxmunkScaled(const boost::shared_ptr<GroundCoxmunk>& Coxmunk, const boost::shared_ptr<GroundBrdfWeight>& Brdf_weight);
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  %python_attribute(coxmunk, boost::shared_ptr<GroundCoxmunk>)
  %python_attribute(brdf_weight, boost::shared_ptr<GroundBrdfWeight>)
  virtual boost::shared_ptr<Ground> clone() const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const { return "GroundCoxmunkScaled"; }
};
}
