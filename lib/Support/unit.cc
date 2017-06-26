#include "unit.h"
#include "fp_exception.h"
#include <cmath>
using namespace FullPhysics;

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
// Don't think this is needed, but leave this commented out for a
// short bit in case we find things breaking and need to come back to
// this. We might need to have a boost version conditional here if
// older versions of boost break. Believe this is the same as 
// phoenix_bind.hpp in older versions of boost.
// #include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

//-----------------------------------------------------------------------
/// Create a unit. This takes string giving the name (e.g., "m",
/// "m^2"), a factor to convert to SI units, and the power of each of
/// the base units.
//-----------------------------------------------------------------------

Unit::Unit(const std::string& Name, double Conversion_to_si, 
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
	    const boost::rational<int>& Sample_index)
: conversion_to_si_(Conversion_to_si),
  name_(Name)
{
  base_unit_powers_[0] = Length_power;
  base_unit_powers_[1] = Mass_power;
  base_unit_powers_[2] = Time_power;
  base_unit_powers_[3] = Current_power;
  base_unit_powers_[4] = Temperature_power;
  base_unit_powers_[5] = Amount_power;
  base_unit_powers_[6] = Luminous_intensity_power;
  base_unit_powers_[7] = Solid_angle_power;
  base_unit_powers_[8] = Angle_power;
  base_unit_powers_[9] = Photon_count_power;
  base_unit_powers_[10] = Sample_index;
}

//-----------------------------------------------------------------------
/// Create a unit. This variation takes name and applies it to a
/// passed in unit. This allows you to do things like:
///  Unit meter_sq("m ^ 2", meter * meter);
//-----------------------------------------------------------------------

Unit::Unit(const std::string& Name, const Unit& Dunit)
: base_unit_powers_(Dunit.base_unit_powers()),
  conversion_to_si_(Dunit.conversion_to_si()),
  name_(Name)
{
}

//-----------------------------------------------------------------------
/// Default constructor. This creates a dimensionless unit.
//-----------------------------------------------------------------------
Unit::Unit()
  :conversion_to_si_(1.0)
{
  for(int i = 0; i < number_base_unit; ++i)
    base_unit_powers_[i] = 0;
}

//-----------------------------------------------------------------------
/// Variation of constructor that uses parse on the given string to
/// determine the units.
//-----------------------------------------------------------------------

Unit::Unit(const std::string& Name_to_parse)
{
  *this = parse(Name_to_parse);
}

//-----------------------------------------------------------------------
/// Multiple two units together.
//-----------------------------------------------------------------------

Unit& Unit::operator*=(const Unit& Dunit)
{
  name_ += " " + Dunit.name();
  conversion_to_si_ *= Dunit.conversion_to_si();
  for(int i = 0; i < number_base_unit; ++i)
    base_unit_powers_[i] += Dunit.base_unit_powers()[i];
  return *this;
}

//-----------------------------------------------------------------------
/// Scale a unit.
//-----------------------------------------------------------------------
Unit& Unit::operator*=(double Scale_factor)
{
  name_ = "";
  conversion_to_si_ *= Scale_factor;
  return *this;
}

Unit& Unit::operator/=(double Scale_factor)
{
  name_ = "";
  conversion_to_si_ /= Scale_factor;
  return *this;
}

//-----------------------------------------------------------------------
/// Divide two units.
//-----------------------------------------------------------------------

Unit& Unit::operator/=(const Unit& Dunit)
{
  name_ += " / (" + Dunit.name() + ")";
  conversion_to_si_ /= Dunit.conversion_to_si();
  for(int i = 0; i < number_base_unit; ++i)
    base_unit_powers_[i] -= Dunit.base_unit_powers()[i];
  return *this;
}

//-----------------------------------------------------------------------
/// Test for equality. This does *not* verify that the name is the
/// same, instead we verify that the underlying units described are
/// the same (i.e., they are commensurate and have same
/// conversion_to_si factor)
//-----------------------------------------------------------------------

bool Unit::operator==(const Unit& U) const 
{
  return (is_commensurate(U) &&
	  fabs(conversion_to_si() / U.conversion_to_si() - 1.0) < 1e-8);
}

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

void Unit::print(std::ostream& Os) const
{
  Os << "Unit:\n"
     << "  Name: " << name_ << "\n"
     << "  Conversion to SI: " << conversion_to_si_ << "\n"
     << "  Unit powers: (";
  for(int i = 0; i < number_base_unit - 1; ++i)
    Os << base_unit_powers_[i] << ", ";
  Os << base_unit_powers_[number_base_unit - 1] << ")\n";
}

//-----------------------------------------------------------------------
// Take a dynamic unit 
//-----------------------------------------------------------------------

Unit FullPhysics::pow(const Unit& Dunit, const boost::rational<int>& Exponent)
{
  return Unit("", 
	      ::pow(Dunit.conversion_to_si(), boost::rational_cast<double>(Exponent)),
	      Dunit.base_unit_powers()[0] * Exponent,
	      Dunit.base_unit_powers()[1] * Exponent,
	      Dunit.base_unit_powers()[2] * Exponent,
	      Dunit.base_unit_powers()[3] * Exponent,
	      Dunit.base_unit_powers()[4] * Exponent,
	      Dunit.base_unit_powers()[5] * Exponent,
	      Dunit.base_unit_powers()[6] * Exponent,
	      Dunit.base_unit_powers()[7] * Exponent,
	      Dunit.base_unit_powers()[8] * Exponent,
	      Dunit.base_unit_powers()[9] * Exponent,
	      Dunit.base_unit_powers()[10] * Exponent);
}

//-----------------------------------------------------------------------
/// Return conversion factor to go from one unit to another. This
/// throws an exception of the units aren't commensurate.
//-----------------------------------------------------------------------

double FullPhysics::conversion(const Unit& Dunit_from, const Unit& Dunit_to)
{
  if(Dunit_from.base_unit_powers() != Dunit_to.base_unit_powers()) {
    Exception e;
    e << "Attempt to convert between Units that are not commensurate.\n"
      << "From unit:\n" << Dunit_from << "To unit:\n" << Dunit_to << "\n";
    throw e;
  }
  return Dunit_from.conversion_to_si() / Dunit_to.conversion_to_si();
}

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

/****************************************************************//**
  Symbol table to parse a base unit.
*******************************************************************/

struct base_unit_ : qi::symbols<char, Unit>
{
  base_unit_()
  {
    using namespace units;
    add
      ("m", m)
      ("meter", m)
      ("Meters", m)
      ("micron", micron)
      ("Microns", micron)
      ("um", micron)
      ("g", g)
      ("gram", g)
      ("s", s)
      ("second", s)
      ("Second", s)
      ("sec", s)
      ("day", day)
      ("year", year)
      ("K", K)
      ("A", A)
      ("mol", mol)
      ("mole", mol)
      ("cd", cd)
      ("sr", sr)
      ("rad", rad)
      ("deg", deg)
      ("Degrees", deg)
      ("ph", ph)
      ("Ph", ph)
      ("sample_index", sample_index)
      ("dimensionless", dimensionless)
      ("N", N)
      ("Pa", Pa)
      ("J", J)
      ("W", W)
      ("Wavenumbers", inv_cm)
      ("Meters per second", m / s)
      ;
  }
} base_unit;

/****************************************************************//**
  Symbol table to parse a prefix.
*******************************************************************/

struct prefix_ : qi::symbols<char, double>
{
  prefix_()
  {
    using namespace units;
    add
      ("Y", 1e24)
      ("yotta", 1e24)
      ("Z", 1e21)
      ("zetta", 1e21)
      ("E", 1e18)
      ("exa", 1e18)
      ("P", 1e15)
      ("peta", 1e15)
      ("T", 1e12)
      ("tera", 1e12)
      ("G", 1e9)
      ("giga", 1e9)
      ("M", 1e6)
      ("mega", 1e6)
      ("k", 1000)
      ("kilo", 1000)
      ("h", 100)
      ("hecto", 100)
      ("da", 10)
      ("deca", 10)
      ("d", 1e-1)
      ("deci", 1e-1)
      ("c", 1e-2)
      ("centi", 1e-2)
      ("m", 1e-3)
      ("milli", 1e-3)
      ("micro", 1e-6)
      ("n", 1e-9)
      ("nano", 1e-9)
      ("p", 1e-12)
      ("pico", 1e-12)
      ("f", 1e-15)
      ("femto", 1e-15)
      ("a", 1e-18)
      ("atto", 1e-18)
      ("z", 1e-21)
      ("zepto", 1e-21)
      ("y", 1e-24)
      ("yocto", 1e-24)
      ;
  }
} prefix;

/****************************************************************//**
  The grammer for parsing units.
*******************************************************************/
typedef std::string::const_iterator iterator_type;
struct unit_parser : qi::grammar<iterator_type, Unit()>
{
  unit_parser() : unit_parser::base_type(exp)
  {
    using namespace qi;
    using namespace boost::phoenix;
    typedef Unit (*pow_f)(const Unit&, 
				const boost::rational<int>&);
    // Eat white space
    ws = *lit(' ');
    // This rule gets a base type, with a prefix.
    prefixunit = 
      lexeme[prefix [_val *= _1] >>
    	     base_unit [_val *= _1]
            ];
    // This gets either a unit, or a unit with a prefix
    unit = prefixunit [_val *= _1] || base_unit [_val *= _1];

    // A term is either a single unit, or a expresion between
    // "(" and ")". Possibly raised to a power
    term = (unit [_val *= _1] || ('(' >> ws >> exp [_val *= _1] >> ws >> ')'))
      // This complicated expression just matches things like
      // "^5" or "^{-1}".
      >> -((ws >> lit('^') >> ws >> int_ [_val = bind((pow_f) &pow,_val, _1)])
	   ||
	   (ws >> lit('^') >> ws >> lit('{') >> ws 
	    >> int_ [_val = bind((pow_f) &pow,_val, _1)]
	    >> ws >> lit('}'))
	   );
    
    // A subexpression is either a term, or a term divided by another term

    subexp = term [_val *= _1] >> *(ws >> lit('/') >> ws >> term [_val /= _1]);

    // A expression is any product of subexpressions, possibly with '*'
    
    exp = ws >> subexp [_val *= _1] % (( ws >> lit('*') >> ws) || +lit(' '))
	     >> ws;
  }
  qi::rule<iterator_type, Unit()> prefixunit;
  qi::rule<iterator_type, Unit()> unit;
  qi::rule<iterator_type, Unit()> term;
  qi::rule<iterator_type, Unit()> subexp;
  qi::rule<iterator_type, Unit()> exp;
  qi::rule<iterator_type> ws;
};

//-----------------------------------------------------------------------
/// Parse a string, and return a Unit that matches the given
/// units. 
//-----------------------------------------------------------------------

Unit Unit::parse(const std::string& S)
{
  unit_parser unitp;
  Unit result = units::dimensionless;
  std::string::const_iterator iter = S.begin();
  bool r = qi::parse(iter, S.end(), unitp, result);
  if(r && iter == S.end()) {
    result.name(S);
    return result;
  }
  Exception e;
  e << "Parsing of unit string '" << S << "' failed.\n"
    << "Parsing stopped at: '" << std::string(iter, S.end()) << "'\n";
  throw e;
}

