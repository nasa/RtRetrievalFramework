// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "pressure_sigma.h"
%}

%base_import(pressure_imp_base)
%import "auto_derivative.i"
%fp_shared_ptr(FullPhysics::PressureSigma);

namespace FullPhysics {
class PressureSigma : public PressureImpBase {
public:
  PressureSigma(const blitz::Array<double, 1>& A,
		const blitz::Array<double, 1>& B,
		double Surface_pressure, bool Pressure_flag);
  PressureSigma(const blitz::Array<double, 1>& Pressure_grid,
		double Surface_pressure, bool Pressure_flag);
  %python_attribute(surface_pressure_uncertainty, double)
  void set_surface_pressure(const AutoDerivative<double>& Surface_pressure);
  void set_levels_from_grid(const blitz::Array<double, 1>& Pressure_grid);
  virtual boost::shared_ptr<Pressure> clone() const;
  virtual std::string state_vector_name_i(int i) const;
  %python_attribute(a, blitz::Array<double, 1>)
  %python_attribute(b, blitz::Array<double, 1>)
protected:
  virtual void calc_pressure_grid() const;
};
}
