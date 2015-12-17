// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "temperature_fixed_level.h"
%}
%base_import(temperature_imp_base)
%import "pressure.i"
%import "pressure_level_input.i"

%fp_shared_ptr(FullPhysics::TemperatureFixedLevel);

namespace FullPhysics {
class TemperatureFixedLevel: public TemperatureImpBase {
public:
  TemperatureFixedLevel(const blitz::Array<bool, 1>& Flag_temp, 
	bool Flag_offset, const blitz::Array<double, 1>& Temp,
	double T_offset, 
        const boost::shared_ptr<Pressure>& Press,
        const boost::shared_ptr<PressureLevelInput>& Press_level);
  %python_attribute(temperature_levels, ArrayAd<double, 1>);
  virtual boost::shared_ptr<Temperature> clone() const;
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const ;
  virtual std::string state_vector_name_i(int i) const;
protected:
  virtual void calc_temperature_grid() const;
};
}
