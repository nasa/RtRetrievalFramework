// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "temperature_met.h"
%}
%base_import(temperature_offset)
%import "meteorology.i"
%import "pressure.i"
%fp_shared_ptr(FullPhysics::TemperatureMet);

namespace FullPhysics {
class TemperatureMet: public TemperatureOffset {
public:
  TemperatureMet(const boost::shared_ptr<Meteorology>& Met_file,
		   const boost::shared_ptr<Pressure>& Press,
		   double Temp_offset,
		   bool Temp_flag);
  virtual boost::shared_ptr<Temperature> clone() const;
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const ;
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(temperature_offset, double)
  %python_attribute(temperature_offset_uncertainty, double)
  %python_attribute(temperature_profile, blitz::Array<double, 1>)
  %python_attribute(pressure_profile, blitz::Array<double, 1>)
protected:
  virtual void calc_temperature_grid() const;
};
}
