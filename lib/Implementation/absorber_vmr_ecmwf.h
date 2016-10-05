#ifndef ABSORBER_VMR_ECMWF_H
#define ABSORBER_VMR_ECMWF_H
#include "absorber_vmr_scaled.h"
#include "ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This class maps the state vector to the absorber VMR on each
  level.

  This particular implementation uses the value from ECMWF file
  (interpolated to the current pressure grid), along with a scale factor.
*******************************************************************/
class AbsorberVmrEcmwf : public AbsorberVmrScaled {
public:
  AbsorberVmrEcmwf(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Scale,                         
		   bool Scale_flag,
		   const std::string& Gas_name);
  virtual ~AbsorberVmrEcmwf() {}

  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<AbsorberVmr> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

  //-----------------------------------------------------------------------
  /// Humidity from ECMWF, used to write to output file
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> specific_humidity_ecmwf() const;

  //-----------------------------------------------------------------------
  /// VMR values converted from specific humidity
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> vmr_profile() const;

  //-----------------------------------------------------------------------
  /// Pressure levels that humidity is on from ECMWF, used to write
  /// to output file 
  //-----------------------------------------------------------------------
  virtual blitz::Array<double, 1> pressure_profile() const;

private:
  boost::shared_ptr<Ecmwf> ecmwf;
};
}
#endif
