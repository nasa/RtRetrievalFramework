// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_coxmunk_plus_lambertian.h"
#include "sub_state_vector_proxy.h"
%}

%base_import(ground)
%base_import(sub_state_vector_proxy)
%import "ground_coxmunk.i"
%import "ground_lambertian.i"

%fp_shared_ptr(FullPhysics::GroundCoxmunkPlusLambertian);
namespace FullPhysics {
class GroundCoxmunkPlusLambertian: public Ground, public SubStateVectorProxy {
public:
  GroundCoxmunkPlusLambertian(const boost::shared_ptr<GroundCoxmunk>& Coxmunk, const boost::shared_ptr<GroundLambertian>& Lambertian);
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  %python_attribute(coxmunk, boost::shared_ptr<GroundCoxmunk>)
  %python_attribute(lambertian, boost::shared_ptr<GroundLambertian>)
  virtual boost::shared_ptr<Ground> clone() const;
  virtual void print(std::ostream& Os) const;
  virtual std::string desc() const { return "GroundCoxmunkPlusLambertian"; }
};
}
