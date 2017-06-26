#ifndef UNIT_H
#define UNIT_H
#include "printable.h"
#include <boost/rational.hpp>
#include <boost/array.hpp>
#include <boost/operators.hpp>

namespace FullPhysics {
/****************************************************************//**
  Libraries such as boost::units allow unit handling where we know the
  units at compile time. This class provide the same sort of handling,
  but for instances where we know the units at runtime rather than
  compile time (e.g., based on input read).

  We do dimensional analysis based on the SI base units. In order,
  these are meter, kilogram, second, Kelvin, Ampere, mole, candela,
  steradian, radian, photon, sample_index

  Note that steradian, radian and sample_index are actually
  dimensionless, but it is useful to track them. Also photon is a
  photon count, which doesn't really have units either. But it is
  useful to track because we can determine the photon count at a
  particular wavelength to convert to cm^-1.
*******************************************************************/
class Unit : public Printable<Unit>,
	     boost::equality_comparable<Unit>,
	     boost::multipliable<Unit>,
	     boost::multipliable<Unit, double>,
	     boost::dividable<Unit>,
	     boost::dividable<Unit, double> {
public:
  enum {number_base_unit = 11};
  Unit(const std::string& Name, double Conversion_to_si, 
	      const boost::rational<int>& Length_power,
	      const boost::rational<int>& Mass_power,
	      const boost::rational<int>& Time_power,
	      const boost::rational<int>& Current_power,
	      const boost::rational<int>& Temperature_power,
	      const boost::rational<int>& Amount_power,
	      const boost::rational<int>& Luminous_intensity_power,
	      const boost::rational<int>& Solid_angle_power,
	      const boost::rational<int>& Angle_power,
              const boost::rational<int>& Photon_count_power,
              const boost::rational<int>& Sample_index);
  Unit(const std::string& Name, const Unit& Dunit);
  Unit(const std::string& Name_to_parse);
  Unit();

//-----------------------------------------------------------------------
/// Array of the powers of the base units (so m^2 would return (1,0,0,0,0,0,0,0))
//-----------------------------------------------------------------------

  const boost::array<boost::rational<int>, number_base_unit>& base_unit_powers() const 
  {return base_unit_powers_;}

//-----------------------------------------------------------------------
/// Conversion factor to go to SI units.
//-----------------------------------------------------------------------

  double conversion_to_si() const {return conversion_to_si_;}

//-----------------------------------------------------------------------
/// Test if this set of units is commensurate with another set. If
/// this return true then conversion() would succeed, otherwise it
/// would fail.
//-----------------------------------------------------------------------

  bool is_commensurate(const Unit& Units) const
  { return base_unit_powers() == Units.base_unit_powers(); }

  static Unit parse(const std::string& S);

//-----------------------------------------------------------------------
/// Name of unit. May be an empty string if a name wasn't assigned.
//-----------------------------------------------------------------------

  const std::string& name() const {return name_;}

//-----------------------------------------------------------------------
/// Set name of unit. 
//-----------------------------------------------------------------------

  void name(const std::string& V) {name_ = V;}

//-----------------------------------------------------------------------
/// Basic math operators for units.
//-----------------------------------------------------------------------
  Unit& operator*=(const Unit& Dunit);
  Unit& operator*=(double Scale_factor);
  Unit& operator/=(double Scale_factor);
  Unit& operator/=(const Unit& Dunit);
  bool operator==(const Unit& U) const;
  void print(std::ostream& Os) const;
private:
  boost::array<boost::rational<int>, number_base_unit> base_unit_powers_;
  double conversion_to_si_;
  std::string name_;
};

Unit pow(const Unit& Dunit, const boost::rational<int>& Exponent);
double conversion(const Unit& Dunit_from, const Unit& Dunit_to);
inline Unit operator/(double Scale_factor, const Unit& Unit)
{
  return Scale_factor * pow(Unit, -1);
}

//-----------------------------------------------------------------------
/// Define various units. Note that if you add something here, you 
/// should also look in unit.cc in the base_unit_ defined
/// for the parser.
//-----------------------------------------------------------------------
  namespace units {
    // base units
    const Unit m("m",     1.0, 1,0,0,0,0,0,0,0,0,0,0);
    const Unit kg("kg",   1.0, 0,1,0,0,0,0,0,0,0,0,0);
    const Unit s("s",     1.0, 0,0,1,0,0,0,0,0,0,0,0);
    const Unit K("K",     1.0, 0,0,0,1,0,0,0,0,0,0,0);
    const Unit A("A",     1.0, 0,0,0,0,1,0,0,0,0,0,0);
    const Unit mol("mol", 1.0, 0,0,0,0,0,1,0,0,0,0,0);
    const Unit cd("cd",   1.0, 0,0,0,0,0,0,1,0,0,0,0);
    const Unit sr("sr",   1.0, 0,0,0,0,0,0,0,1,0,0,0);
    const Unit rad("rad", 1.0, 0,0,0,0,0,0,0,0,1,0,0);
    const Unit ph("ph",   1.0, 0,0,0,0,0,0,0,0,0,1,0);
    const Unit sample_index("sample_index", 
                          1.0, 0,0,0,0,0,0,0,0,0,0,1);
    const Unit dimensionless("dimensionless", 
                          1.0, 0,0,0,0,0,0,0,0,0,0,0);

//-----------------------------------------------------------------------
/// Pi
//-----------------------------------------------------------------------

    const double pi = 3.14159265358979323846;
    
    const Unit deg("deg", pi / 180.0 * rad);
    const Unit arcsec("arcsec", deg / (60.0 * 60));

    // Some useful names
    const Unit day("day", 24 * 60 * 60 * s);
    const Unit year("year", 365.25 * 24 * 60 * 60 * s);
    const Unit cm("cm", 1e-2 * m);
    const Unit km("km", 1e3 * m);
    const Unit nm("nm", 1e-9 * m);
    const Unit micron("micron", 1e-6 * m);
    const Unit g("g", 1e-3 * kg);
    const Unit inv_cm("cm^-1", pow(cm, -1));
    const Unit inv_sr("sr^-1", pow(sr, -1));
    // Derived units
    const Unit N("N", kg * m / (s*s));
    const Unit Pa("Pa", N / (m * m));
    const Unit J("J", N * m);
    const Unit W("W", J / s);
  }
}
#endif
