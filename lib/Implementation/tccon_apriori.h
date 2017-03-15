#ifndef TCCON_APRIORI_H
#define TCCON_APRIORI_H
#include "meteorology.h"
#include "pressure.h"
#include "temperature.h"
#include "level_1b.h"

namespace FullPhysics {
/****************************************************************//**
  This class is used to calculate a CO2 apriori using the same 
  method TCCON does.
*******************************************************************/

class TcconApriori : public Printable<TcconApriori> {
public:
  TcconApriori(const boost::shared_ptr<Meteorology>& Met_file,
	       const boost::shared_ptr<Level1b>& L1b_file,
	       double Co2_ref = 0.000380,
	       const Time& Ref_time = Time::parse_time("2005-01-01T00:00:00Z"),
	       DoubleWithUnit Rate_increase = 
	       DoubleWithUnit(0.005, "year^-1"),
	       DoubleWithUnit Age_air_pbl = 
	       DoubleWithUnit(0.2, "year"),
	       DoubleWithUnit Age_air_tropopause = 
	       DoubleWithUnit(0.4, "year"),
	       DoubleWithUnit Age_air_upper_stratosphere = 
	       DoubleWithUnit(5.5, "year"));

  TcconApriori(const boost::shared_ptr<Level1b>& L1b_file,
	       const blitz::Array<double, 1>& Press_profile,
	       const blitz::Array<double, 1>& Temp_profile,
	       double Surface_pressure,
	       double Co2_ref = 0.000380,
	       const Time& Ref_time = Time::parse_time("2005-01-01T00:00:00Z"),
	       DoubleWithUnit Rate_increase = 
	       DoubleWithUnit(0.005, "year^-1"),
	       DoubleWithUnit Age_air_pbl = 
	       DoubleWithUnit(0.2, "year"),
	       DoubleWithUnit Age_air_tropopause = 
	       DoubleWithUnit(0.4, "year"),
	       DoubleWithUnit Age_air_upper_stratosphere = 
	       DoubleWithUnit(5.5, "year"));

  TcconApriori(const boost::shared_ptr<Level1b>& L1b_file,
	       const boost::shared_ptr<Pressure>& Pressure,
	       const boost::shared_ptr<Temperature>& Temp,
	       double Co2_ref = 0.000380,
	       const Time& Ref_time = Time::parse_time("2005-01-01T00:00:00Z"),
	       DoubleWithUnit Rate_increase = 
	       DoubleWithUnit(0.005, "year^-1"),
	       DoubleWithUnit Age_air_pbl = 
	       DoubleWithUnit(0.2, "year"),
	       DoubleWithUnit Age_air_tropopause = 
	       DoubleWithUnit(0.4, "year"),
	       DoubleWithUnit Age_air_upper_stratosphere = 
	       DoubleWithUnit(5.5, "year"));

  virtual ~TcconApriori() {}
  double co2_vmr(double P) const;
  blitz::Array<double, 1> co2_vmr_grid(const Pressure& P) const;
  double tropopause_pressure() const;
  double planetary_boundary_layer_pressure() const;
  double fractional_amplitude_seasonal_cycle() const;
  DoubleWithUnit age_air(double P) const;

  void print(std::ostream& Os) const { Os << "TcconApriori"; }

private:

  blitz::Array<double, 1> press_profile;
  blitz::Array<double, 1> temp_profile;
  double surface_pressure;
    
  double trop_pres, pbl_pres, fasc;
  boost::shared_ptr<Level1b> l1b;
  double co2_ref;
  Time co2_ref_time;
  DoubleWithUnit rate_increase;
  DoubleWithUnit age_air_pbl;
  DoubleWithUnit age_air_trop;
  DoubleWithUnit age_air_strat;
};
}

#endif

