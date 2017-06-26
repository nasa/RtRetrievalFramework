// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "ground_coxmunk.h"
#include "sub_state_vector_array.h"
%}

%base_import(ground)
%base_import(sub_state_vector_array)

%fp_shared_ptr(FullPhysics::GroundCoxmunk);
namespace FullPhysics {
class GroundCoxmunk: public SubStateVectorArray<Ground> {
public:
  GroundCoxmunk(const double Windspeed,
                const bool& Ws_flag, 
                const blitz::Array<double, 1>& Refr_index);
  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;
  virtual const AutoDerivative<double> windspeed() const;
  virtual const double refractive_index(const int Spec_idx) const;
  virtual boost::shared_ptr<Ground> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  virtual void print(std::ostream& Os) const;
  virtual void update_sub_state_hook();
};
}
