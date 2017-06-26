#ifndef ABSORBER_VMR_FIXED_LEVEL_H
#define ABSORBER_VMR_FIXED_LEVEL_H
#include "absorber_vmr_imp_base.h"
#include "pressure.h"
#include "pressure_level_input.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This implementation just gets the VMR from each level from the 
  state vector. 
*******************************************************************/
class AbsorberVmrFixedLevel : public AbsorberVmrImpBase {
public:
  AbsorberVmrFixedLevel(const boost::shared_ptr<Pressure>& Press,
	const boost::shared_ptr<PressureLevelInput>& Press_level,	   
	const blitz::Array<bool, 1>& Used_flag, 
	const blitz::Array<double, 1>& Vmr,
	const std::string& Gas_name);
  virtual ~AbsorberVmrFixedLevel() {}
  virtual void print(std::ostream& Os) const;

  virtual std::string state_vector_name_i(int i) const
  { return gas_name() + " VMR for Press Lvl " + 
      boost::lexical_cast<std::string>(i + 1); }
  virtual boost::shared_ptr<AbsorberVmr> clone() const;
  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

//-----------------------------------------------------------------------
/// Volume mixing ratio on the fixed pressure levels.
//-----------------------------------------------------------------------

  blitz::Array<double, 1> volume_mixing_ratio_level() const 
  {return coeff.value();}

//-----------------------------------------------------------------------
/// Volume mixing ratio on the fixed pressure levels, restricted to the 
/// active levels.
//-----------------------------------------------------------------------

  blitz::Array<double, 1> volume_mixing_ratio_active_level() const 
  {return coeff.value()(blitz::Range(0, press->number_level() - 1));}
protected:
  virtual void calc_vmr() const;
private:
  boost::shared_ptr<PressureLevelInput> press_level;
};
}
#endif
