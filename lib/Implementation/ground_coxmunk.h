#ifndef GROUND_COXMUNK_H
#define GROUND_COXMUNK_H

#include "ground.h"
#include "auto_derivative.h"
#include "sub_state_vector_array.h"

namespace FullPhysics {
/****************************************************************//**
  This class implements a Coxmunk ground type. 
*******************************************************************/
class GroundCoxmunk: public SubStateVectorArray<Ground> {
public:
  GroundCoxmunk(const double Windspeed,
                const bool& Ws_flag, 
                const blitz::Array<double, 1>& Refr_index);

  virtual ArrayAd<double, 1> surface_parameter(const double wn, const int spec_index) const;

  virtual const AutoDerivative<double> windspeed() const;
  virtual const double windspeed_uncert() const;

  virtual const double refractive_index(const int Spec_idx) const;

  blitz::Array<double, 1> kernel_value_params(const int Spec_index);
  virtual const double kernel_value(const int Spec_index, const double Sza,
                                    const double Vza, const double Azm);
  virtual const blitz::Array<double, 1> kernel_value_giss
                                   (const int Spec_index, const double Sza,
                                    const double Vza, const double Azm);
  virtual const blitz::Array<AutoDerivative<double>, 1> kernel_value_l_giss
                                   (const int Spec_index, const double Sza,
                                    const double Vza, const double Azm);

  virtual boost::shared_ptr<Ground> clone() const;

  virtual std::string state_vector_name_i(int i) const;

  virtual void print(std::ostream& Os) const;

  virtual std::string desc() const { return "GroundCoxmunk"; }
  
  virtual void update_sub_state_hook();

private:

  blitz::Array<double, 1> refractive_index_;

};
}
#endif
