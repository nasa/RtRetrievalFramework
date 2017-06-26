#include "tccon_apriori.h"
#include "pressure_sigma.h"
#include "temperature_level_offset.h"
#include "altitude_hydrostatic.h"
using namespace FullPhysics;
using namespace blitz;
using namespace boost::posix_time;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(TcconApriori)
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double, const Time&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double, const Time&, const DoubleWithUnit&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double, const Time&, const DoubleWithUnit&,
			  const DoubleWithUnit&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double, const Time&, const DoubleWithUnit&,
			  const DoubleWithUnit&, 
			  const DoubleWithUnit&>())
.def(luabind::constructor<const boost::shared_ptr<Meteorology>&,
			  const boost::shared_ptr<Level1b>&,
			  double, const Time&, const DoubleWithUnit&,
			  const DoubleWithUnit&, 
			  const DoubleWithUnit&,
			  const DoubleWithUnit&>())
.def(luabind::constructor<const boost::shared_ptr<Level1b>&,
			  const blitz::Array<double, 1>&,
			  const blitz::Array<double, 1>&,
			  double>())
.def(luabind::constructor<const boost::shared_ptr<Level1b>&,
			  const boost::shared_ptr<Pressure>&,
			  const boost::shared_ptr<Temperature>&>())
.def("co2_vmr_grid", &TcconApriori::co2_vmr_grid)
.def("tropopause_pressure", &TcconApriori::tropopause_pressure)
.def("planetary_boundary_layer_pressure", &TcconApriori::planetary_boundary_layer_pressure)
.def("fractional_amplitude_seasonal_cycle", &TcconApriori::fractional_amplitude_seasonal_cycle)
.def("age_air", &TcconApriori::age_air)
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. Latitude should be in degrees, height in meters.
//-----------------------------------------------------------------------

TcconApriori::TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
			   const boost::shared_ptr<Level1b>& L1b_file,
			   double Co2_ref,
			   const Time& Ref_time,
			   DoubleWithUnit Rate_increase,
			   DoubleWithUnit Age_air_pbl,
			   DoubleWithUnit Age_air_tropopause,
			   DoubleWithUnit Age_air_upper_stratosphere
)
: l1b(L1b_file),
  co2_ref(Co2_ref), co2_ref_time(Ref_time), rate_increase(Rate_increase),
  age_air_pbl(Age_air_pbl), age_air_trop(Age_air_tropopause),
  age_air_strat(Age_air_upper_stratosphere)
{
  surface_pressure = Met_file->surface_pressure();

  Array<double, 1> met_press = Met_file->pressure_levels();
  press_profile.resize(met_press.rows());
  press_profile = met_press;

  Array<double, 1> met_temp = Met_file->temperature();
  temp_profile.resize(met_temp.rows());
  temp_profile = met_temp;

  trop_pres = tropopause_pressure();
  pbl_pres = planetary_boundary_layer_pressure();
  fasc = fractional_amplitude_seasonal_cycle();
}

TcconApriori::TcconApriori(const boost::shared_ptr<Level1b>& L1b_file,
			   const blitz::Array<double, 1>& Press_profile,
			   const blitz::Array<double, 1>& Temp_profile,
			   double Surface_pressure,
			   double Co2_ref,
			   const Time& Ref_time,
			   DoubleWithUnit Rate_increase,
			   DoubleWithUnit Age_air_pbl,
			   DoubleWithUnit Age_air_tropopause,
			   DoubleWithUnit Age_air_upper_stratosphere
)
: press_profile(Press_profile), temp_profile(Temp_profile), 
  surface_pressure(Surface_pressure),
  l1b(L1b_file), co2_ref(Co2_ref), co2_ref_time(Ref_time), 
  rate_increase(Rate_increase),
  age_air_pbl(Age_air_pbl), age_air_trop(Age_air_tropopause),
  age_air_strat(Age_air_upper_stratosphere)
{

  trop_pres = tropopause_pressure();
  pbl_pres = planetary_boundary_layer_pressure();
  fasc = fractional_amplitude_seasonal_cycle();
}

TcconApriori::TcconApriori(const boost::shared_ptr<Level1b>& L1b_file,
			   const boost::shared_ptr<Pressure>& Pressure,
			   const boost::shared_ptr<Temperature>& Temp,
			   double Co2_ref,
			   const Time& Ref_time,
			   DoubleWithUnit Rate_increase,
			   DoubleWithUnit Age_air_pbl,
			   DoubleWithUnit Age_air_tropopause,
			   DoubleWithUnit Age_air_upper_stratosphere
)
: l1b(L1b_file),
  co2_ref(Co2_ref), co2_ref_time(Ref_time), rate_increase(Rate_increase),
  age_air_pbl(Age_air_pbl), age_air_trop(Age_air_tropopause),
  age_air_strat(Age_air_upper_stratosphere)
{
  press_profile.resize(Pressure->pressure_grid().rows());
  temp_profile.resize(press_profile.shape());
  for(int lev = 0; lev < press_profile.rows(); lev++) {
    AutoDerivativeWithUnit<double> p_ad = 
      Pressure->pressure_grid()(lev).convert(units::Pa);
    AutoDerivativeWithUnit<double> t_ad = 
      Temp->temperature(p_ad).convert(units::K);
    press_profile(lev) = p_ad.value.value();
    temp_profile(lev) = t_ad.value.value();
  }
  surface_pressure = Pressure->surface_pressure_value();

  trop_pres = tropopause_pressure();
  pbl_pres = planetary_boundary_layer_pressure();
  fasc = fractional_amplitude_seasonal_cycle();
}

//-----------------------------------------------------------------------
/// Calculate the Tropopause pressure.
//-----------------------------------------------------------------------

double TcconApriori::tropopause_pressure() const
{
  const double lapse_rate_threshold = -2;
  Array<double, 1> b(press_profile.shape());
  b = 0;
  boost::shared_ptr<PressureSigma> p
    (new PressureSigma(press_profile, b, surface_pressure, false));
  boost::shared_ptr<TemperatureLevelOffset> t
    (new TemperatureLevelOffset(p, temp_profile, 0, false));
  AltitudeHydrostatic a(p, t, l1b->latitude(0), 
			l1b->altitude(0));
  for(int i = press_profile.rows() - 1; i > 0; --i)
    if(press_profile(i) < 50000) {
      double lapse_rate = (temp_profile(i - 1) - temp_profile(i)) / 
	(a.altitude(AutoDerivativeWithUnit<double>(press_profile(i - 1), units::Pa)).convert(units::km).value.value() - a.altitude(AutoDerivativeWithUnit<double>(press_profile(i), units::Pa)).convert(units::km).value.value());
      if(lapse_rate > lapse_rate_threshold)
	return press_profile(i);
    }
  throw Exception("Couldn't find Tropopause pressure");
}

//-----------------------------------------------------------------------
/// Calculate planetary boundary layer pressure.
///
/// This assume that the daytime high-latitude PBL varies
/// from 90 mb in the tropics to 280 mb at high latitudes
//-----------------------------------------------------------------------

double TcconApriori::planetary_boundary_layer_pressure() const
{
  double lat = l1b->latitude(0).convert(units::rad).value;
  double pbl_pressure_atm = 
    (0.70 - 0.15 * cos(2 * lat) - 0.1 * sin(lat) * 
     sin(2 * units::pi * (l1b->time(0).frac_day_of_year() - 110)/365.25));
  return pbl_pressure_atm * 101325;
}

//-----------------------------------------------------------------------
/// Calculate fractional amplitude of seasonal cycle at surface
//-----------------------------------------------------------------------

double TcconApriori::fractional_amplitude_seasonal_cycle() const
{
  double lat = l1b->latitude(0).convert(units::deg).value;
  return 0.01 * sin(2 * lat * (1 - lat/720) * units::pi / 180) * 
    exp(lat / 45);
}

//-----------------------------------------------------------------------
/// Calculate the age of air for the given pressure level (pressure is
/// in pascals).
//-----------------------------------------------------------------------

DoubleWithUnit TcconApriori::age_air(double P) const
{
  double surf = surface_pressure;
  if(P > surf)
    return DoubleWithUnit(0, "year");
  if(P > pbl_pres)
    return age_air_pbl * (surf - P) / (surf - pbl_pres);
  if(P > trop_pres)
    return age_air_pbl + (age_air_trop - age_air_pbl) * 
      (pbl_pres - P) / (pbl_pres - trop_pres);
  return (age_air_strat - (age_air_strat - age_air_trop) *
	  sqrt(P / trop_pres));
}

//-----------------------------------------------------------------------
/// Return the CO2 volume mixing ratio for the given pressure level.
/// The pressure should be in Pascals.
//-----------------------------------------------------------------------

double TcconApriori::co2_vmr(double P) const
{
  DoubleWithUnit age_of_air = age_air(P);
  DoubleWithUnit time_air_from_ref = 
    DoubleWithUnit(l1b->time(0) - co2_ref_time, "s") - age_of_air;
  DoubleWithUnit vmr0 = co2_ref * (1 + rate_increase * time_air_from_ref);
  if(!vmr0.units.is_commensurate(units::dimensionless))
    throw Exception("Units are mangled");
  double sdma = sin(2 * units::pi * ((l1b->time(0).frac_day_of_year() + 75) / 365.2 - 
				     age_of_air.value));
  // Not sure where these magic numbers come from.
  sdma = 1.45 - exp(-1.11 * sdma);
  return vmr0.value * 
    (1 + fasc * exp(-((age_of_air) / DoubleWithUnit(0.25, "year")).value) * 
     sdma);
}

//-----------------------------------------------------------------------
/// Return CO2 VMR for each pressure in the given pressure grid.
//-----------------------------------------------------------------------

blitz::Array<double, 1> TcconApriori::co2_vmr_grid(const Pressure& P) const
{
  Array<double, 1> p(P.pressure_grid().convert(units::Pa).value.value());
  Array<double, 1> res(p.shape());
  for(int i = 0; i < p.rows(); ++i)
    res(i) = co2_vmr(p(i));
  return res;
}
