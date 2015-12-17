#ifndef ABSORBER_VMR_FIXED_LEVEL_SCALED_H
#define ABSORBER_VMR_FIXED_LEVEL_SCALED_H
#include "absorber_vmr_imp_base.h"
#include "pressure_level_input.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This implementation has the VMR passed to the constructor, and 
  applies a scale factor from the state vector.
*******************************************************************/
class AbsorberVmrFixedLevelScaled : public AbsorberVmrImpBase {
public:
  AbsorberVmrFixedLevelScaled(const boost::shared_ptr<Pressure>& Press,
        const boost::shared_ptr<PressureLevelInput>& Press_level,
        const blitz::Array<double, 1>& Vmr,
	bool Used_flag,
        double Scale,                         
        const std::string& Gas_name);
  virtual ~AbsorberVmrFixedLevelScaled() {}
  virtual void print(std::ostream& Os) const;
  virtual std::string state_vector_name_i(int i) const
  { return gas_name() + " Scaling factor"; }
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

//-----------------------------------------------------------------------
/// Scale factor.
//-----------------------------------------------------------------------

  double scale_factor() const { return coeff(0).value(); }

//-----------------------------------------------------------------------
/// Uncertainty of scale factor.
//-----------------------------------------------------------------------

  double scale_uncertainty() const
  { return (sv_cov_sub.rows() > 0 ? sqrt(sv_cov_sub(0,0)) : 0); }
protected:
  virtual void calc_vmr() const;
private:
  boost::shared_ptr<PressureLevelInput> press_level;
  blitz::Array<double, 1> vmr0;
};
}
#endif
