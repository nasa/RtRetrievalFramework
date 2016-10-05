// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"
%{
#include "temperature_offset.h"
%}
%base_import(temperature_imp_base)

%fp_shared_ptr(FullPhysics::TemperatureOffset);

namespace FullPhysics {
class TemperatureOffset: public TemperatureImpBase {
public:
  TemperatureOffset(const boost::shared_ptr<Pressure>& Press,
		    double Temp_offset,
		    bool Temp_flag);
  virtual boost::shared_ptr<Temperature> clone() const;
  virtual boost::shared_ptr<Temperature> 
  clone(const boost::shared_ptr<Pressure>& Press) const ;
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(temperature_offset, double);
  %python_attribute(temperature_offset_uncertainty, double);
protected:
  virtual void calc_temperature_grid() const = 0;
};
}
