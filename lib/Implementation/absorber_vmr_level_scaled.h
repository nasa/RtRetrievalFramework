#ifndef ABSORBER_VMR_LEVEL_SCALED_H
#define ABSORBER_VMR_LEVEL_SCALED_H
#include "absorber_vmr_scaled.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses a passed vmr profile
  (interpolated to the current pressure grid), along with a scale factor.
*******************************************************************/
class AbsorberVmrLevelScaled : public AbsorberVmrScaled {
public:
  AbsorberVmrLevelScaled(const boost::shared_ptr<Pressure>& Press,
			 const blitz::Array<double, 1>& Vmr_profile,
			 double Scale,                         
			 bool Scale_flag,
			 const std::string& Gas_name);
  virtual ~AbsorberVmrLevelScaled() {}

  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

  //-----------------------------------------------------------------------
  /// VMR values passed in from input
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> vmr_profile() const;

  //-----------------------------------------------------------------------
  /// Pressure levels that vmr is on
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> pressure_profile() const;

private:
  blitz::Array<double, 1> vmr_profile_;
};
}
#endif
