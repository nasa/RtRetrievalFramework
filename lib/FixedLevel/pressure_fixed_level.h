#ifndef PRESSURE_FIXED_LEVEL_H
#define PRESSURE_FIXED_LEVEL_H
#include "state_vector.h"
#include "observer.h"
#include "fp_exception.h"
#include "pressure_level_input.h"
#include "pressure_imp_base.h"

namespace FullPhysics {
/****************************************************************//**
  This class maintains the pressure portion of the state. This
  particular implementation has a fixed set of pressure levels, with
  only the surface pressure changing. As the surface pressure changes,
  it may pass a pressure level, changing the number of levels that lie
  above the surface.
*******************************************************************/
class PressureFixedLevel : public PressureImpBase {
public:
  PressureFixedLevel(bool Pressure_flag, 
		     const boost::shared_ptr<PressureLevelInput>& Press_level,
		     double Surface_pressure)
    : press_level(Press_level)
    {
      blitz::Array<double, 1> coeff(1);
      blitz::Array<bool, 1> used(1);
      coeff(0) = Surface_pressure;
      used(0) = Pressure_flag;
      init(coeff, used);
      cov.resize(1, 1);
      cov(0,0) = 1;
    }
  virtual ~PressureFixedLevel() {}

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
    const blitz::Array<double, 1>& plev = press_level->pressure_level();
    range_check(Surface_pressure.value(), plev(0), 
		plev(plev.rows() - 1));
    Observable<Pressure>::notify_update_do(*this);
  }

//-----------------------------------------------------------------------
/// Number of active levels, this is just the size of pressure_grid
/// for the current surface pressure.
//-----------------------------------------------------------------------

  int number_active_level() const 
  { return pressure_grid().value.rows(); }

//-----------------------------------------------------------------------
/// Number of active layers. This is 1 less than the number of
/// levels, since the levels give the top an bottom of a layer.
//-----------------------------------------------------------------------

  int number_active_layer() const { return number_active_level() - 1; }

  static const double new_level_fractional_size;

//-----------------------------------------------------------------------
/// Maximum number of levels that we can have.
//-----------------------------------------------------------------------

  int max_number_level() const {return press_level->pressure_level().rows();}

//-----------------------------------------------------------------------
/// Return the pressure on the fixed levels, but only include the
/// "active" portion. This is all the pressure levels above the
/// surface, plus the last one at or below the surface. This is the
/// same size as pressure_grid, but differs in that the bottom level
/// typically lies below the surface
//-----------------------------------------------------------------------

  blitz::Array<double, 1> pressure_active_levels() const
  { return press_level->pressure_level()
      (blitz::Range(0, number_active_level() - 1)); }

  virtual void print(std::ostream& Os) const;

  virtual boost::shared_ptr<Pressure> clone() const;
  virtual std::string state_vector_name_i(int i) const
  { return "Surface Pressure (Pascals)"; }
protected:
  virtual void calc_pressure_grid() const;
private:
  boost::shared_ptr<PressureLevelInput> press_level;
};
}
#endif
