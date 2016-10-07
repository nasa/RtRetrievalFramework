#ifndef PRESSURE_SIGMA_H
#define PRESSURE_SIGMA_H
#include "state_vector.h"
#include "observer.h"
#include "fp_exception.h"
#include "pressure_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the pressure portion of the state. This
  particular implementation use pressure sigma levels to determine 
  the pressure levels as the surface pressure changes.
*******************************************************************/
class PressureSigma : public PressureImpBase {
public:
  PressureSigma(const blitz::Array<double, 1>& A,
		const blitz::Array<double, 1>& B,
		double Surface_pressure, bool Pressure_flag);

  PressureSigma(const blitz::Array<double, 1>& Pressure_grid,
		double Surface_pressure, bool Pressure_flag);

  virtual ~PressureSigma() {}

//-----------------------------------------------------------------------
/// Return the current surface pressure uncertainty. This is in Pascals.
//-----------------------------------------------------------------------

  double surface_pressure_uncertainty() const 
  {return (cov(0, 0) < 0 ? 0.0 : sqrt(cov(0, 0))); }

//-----------------------------------------------------------------------
/// Set the surface pressure. This is in Pascals.
//-----------------------------------------------------------------------

  void set_surface_pressure(const AutoDerivative<double>& Surface_pressure) 
  {
    coeff(0) = Surface_pressure;
    cache_stale = true;
    Observable<Pressure>::notify_update_do(*this);
  }
  
  void set_levels_from_grid(const blitz::Array<double, 1>& Pressure_grid);
 
  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Pressure> clone() const;
  virtual std::string state_vector_name_i(int i) const
  { return "Surface Pressure (Pascals)"; }
  const blitz::Array<double, 1>& a() const {return a_;}
  const blitz::Array<double, 1>& b() const {return b_;}
protected:
  virtual void calc_pressure_grid() const;
private:
  blitz::Array<double, 1> a_, b_;
};
}
#endif
