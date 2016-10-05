// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "pressure_fixed_level.h"
%}

%base_import(pressure_imp_base)
%import "pressure_level_input.i"

%fp_shared_ptr(FullPhysics::PressureFixedLevel);

namespace FullPhysics {
  class PressureFixedLevel;

class PressureFixedLevel : public PressureImpBase {
public:
  PressureFixedLevel(bool Pressure_flag, 
		     const boost::shared_ptr<PressureLevelInput>& Press_level,
		     double Surface_pressure);
  %python_attribute(surface_pressure_uncertainty, double);
  void set_surface_pressure(const AutoDerivative<double>& Surface_pressure);
  %python_attribute(number_active_level, int)
  %python_attribute(number_active_layer, int)
  %python_attribute(max_number_level, int)
  %python_attribute(pressure_active_levels, blitz::Array<double, 1>)
  virtual boost::shared_ptr<Pressure> clone() const;
  virtual std::string state_vector_name_i(int i) const;
protected:
  virtual void calc_pressure_grid() const;
};
}
