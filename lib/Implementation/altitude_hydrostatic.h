#ifndef ALTITUDE_HYDROSTATIC_H
#define ALTITUDE_HYDROSTATIC_H
#include "altitude.h"
#include "pressure.h"
#include "temperature.h"
#include "double_with_unit.h"
#include "linear_interpolate.h"

namespace FullPhysics {
/****************************************************************//**
  This class handles the calculation of the altitude an gravity
  constants, automatically updating with the surface pressure or
  temperature profile is updated.

  We do this by solving the hydrostatic equations.

  \todo reference for this? (ATB?)
*******************************************************************/
class AltitudeHydrostatic : public Observer<Temperature>,
                 public Observer<Pressure>,
                 public Altitude
{
public:
  AltitudeHydrostatic(const boost::shared_ptr<Pressure>& P,
                      const boost::shared_ptr<Temperature>& T,
                      const DoubleWithUnit& Latitude, 
                      const DoubleWithUnit& Surface_height,
                      const int Num_sublayer = 10);
  virtual ~AltitudeHydrostatic() {};
  using Observer<Temperature>::notify_update;
  using Observer<Pressure>::notify_update;

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the pressure is changed, so we observe press
/// and mark the cache when it changes. 
//-----------------------------------------------------------------------

  virtual void notify_update(const Pressure& P)
  {
    cache_is_stale = true;   
  }

//-----------------------------------------------------------------------
/// For performance, we cache some data as we calculate it. This
/// becomes stale when the temperature is changed, so we observe temperature
/// and mark the cache when it changes. 
//-----------------------------------------------------------------------

  virtual void notify_update(const Temperature& T)
  {
    cache_is_stale = true;   
  }

  virtual AutoDerivativeWithUnit<double> 
  altitude(const AutoDerivativeWithUnit<double>& P) const
  { calc_alt_and_grav(); 
    AutoDerivativeWithUnit<double> p_pas = P.convert(units::Pa);
    return AutoDerivativeWithUnit<double>((*alt)(p_pas.value), units::km); }
  virtual AutoDerivativeWithUnit<double> 
  gravity(const AutoDerivativeWithUnit<double>& P) const
  { calc_alt_and_grav(); 
    AutoDerivativeWithUnit<double> p_pas = P.convert(units::Pa);
    return AutoDerivativeWithUnit<double>((*grav)(p_pas.value), "m/s^2"); 
  }
  virtual void print(std::ostream& Os) const { Os << "AltitudeHydrostatic"; }
  virtual boost::shared_ptr<Altitude> clone() const
  { boost::shared_ptr<Pressure> pnew = p->clone();
    return clone(pnew, t->clone(pnew)); }
  virtual boost::shared_ptr<Altitude> 
  clone(const boost::shared_ptr<Pressure>& Press,
        const boost::shared_ptr<Temperature>& Temp) const;
private:
  typedef LinearInterpolate<AutoDerivative<double>, AutoDerivative<double> >
    lin_type;
  mutable bool cache_is_stale;
  mutable boost::shared_ptr<lin_type> alt;
  mutable boost::shared_ptr<lin_type> grav;
  DoubleWithUnit latitude, surface_height;
  boost::shared_ptr<Pressure> p;
  boost::shared_ptr<Temperature> t;
  int num_sublayer;
  void calc_alt_and_grav() const;
  AutoDerivative<double> 
  gravity_calc(double gdlat, AutoDerivative<double> altit) const;
  void altitude_calc
  (const ArrayAdWithUnit<double, 1>& press_grid,  
   const ArrayAdWithUnit<double, 1>& temp_grid) const;
};

}
#endif
