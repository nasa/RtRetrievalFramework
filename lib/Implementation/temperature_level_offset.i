// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "temperature_level_offset.h"
%}
%base_import(temperature_offset)

%fp_shared_ptr(FullPhysics::TemperatureLevelOffset);

namespace FullPhysics {
class TemperatureLevelOffset: public TemperatureOffset {
public:
  TemperatureLevelOffset(const boost::shared_ptr<Pressure>& Press,
		   const blitz::Array<double, 1>& Temp_levels,
		   double Temp_offset,
		   bool Temp_flag);
  virtual boost::shared_ptr<Temperature> clone() const;
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const ;
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(temperature_offset, double)
  %python_attribute(temperature_offset_uncertainty, double)
protected:
  virtual void calc_temperature_grid() const;
};
}
