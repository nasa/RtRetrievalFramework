#ifndef ABSORBER_VMR_MET_H
#define ABSORBER_VMR_MET_H
#include "absorber_vmr_scaled.h"
#include "meteorology.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses the value from MET file
  (interpolated to the current pressure grid), along with a scale factor.
*******************************************************************/
class AbsorberVmrMet : public AbsorberVmrScaled {
public:
  AbsorberVmrMet(const boost::shared_ptr<Meteorology>& Met_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Scale,                         
		   bool Scale_flag,
		   const std::string& Gas_name);
  virtual ~AbsorberVmrMet() {}

  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

  //-----------------------------------------------------------------------
  /// Humidity from MET, used to write to output file
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> specific_humidity() const;

  //-----------------------------------------------------------------------
  /// VMR values converted from specific humidity
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> vmr_profile() const;

  //-----------------------------------------------------------------------
  /// Pressure levels that humidity is on from MET, used to write
  /// to output file 
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> pressure_profile() const;

private:
  boost::shared_ptr<Meteorology> met;
};
}
#endif
