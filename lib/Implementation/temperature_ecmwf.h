#ifndef TEMPERATURE_ECMWF_H
#define TEMPERATURE_ECMWF_H
#include "temperature_offset.h"
#include "ecmwf.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the temperature portion of the state. This
  particular implementation uses the temperature from ECMWF file
  (interpolated to the current pressure grid), along with an offset.
*******************************************************************/
class TemperatureEcmwf: public TemperatureOffset {
public:
  TemperatureEcmwf(const boost::shared_ptr<Ecmwf>& Ecmwf_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Temp_offset,
		   bool Temp_flag);
  virtual ~TemperatureEcmwf() {}
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Temperature> clone() const
  { return clone(press->clone()); }
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const;

//-----------------------------------------------------------------------
/// Temperature from ECMWF, used to write to output file
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> temperature_profile() const
  { 
    blitz::Array<double, 1> t1, p;
    ecmwf->temperature_grid(p, t1);
    return t1;
  }

//-----------------------------------------------------------------------
/// Pressure levels that temperature is on from ECMWF, used to write
/// to output file 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> pressure_profile() const
  { 
    blitz::Array<double, 1> t1, p;
    ecmwf->temperature_grid(p, t1);
    return p;
  }

private:
  boost::shared_ptr<Ecmwf> ecmwf;
};
}
#endif
