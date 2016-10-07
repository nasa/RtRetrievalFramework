#ifndef OLD_CONSTANT_H
#define OLD_CONSTANT_H
#include "double_with_unit.h"
namespace FullPhysics {
/****************************************************************//**
  Contains useful constants.
*******************************************************************/

///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// This will get replaced with the new constants class
///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  namespace OldConstant {
//-----------------------------------------------------------------------
/// Avogadro constant
//-----------------------------------------------------------------------

const double avogadro_constant = 6.02214179e23;

//-----------------------------------------------------------------------
/// Molar weight of dry air, in g mol^-1
//-----------------------------------------------------------------------

const double molar_weight_dry_air = 28.9644;

//-----------------------------------------------------------------------
/// Molar weight of water, in g mol^-1
//-----------------------------------------------------------------------

const double molar_weight_water = 18.01528;


//-----------------------------------------------------------------------
/// Molar weight of O3, in g mol^-1
//-----------------------------------------------------------------------

const double molar_weight_ozone = 47.9982;

//-----------------------------------------------------------------------
/// Gas constant from the ideal gas law, in J mol-1 K-1
//-----------------------------------------------------------------------

const double gas_constant = 8.3144621; 

// Bring pi into this namespace, but make sure it matches the value
// already defined in units namespace.
const double pi = units::pi;

//-----------------------------------------------------------------------
/// Speed of light, in m/s. Note that this is actually exact, since the
/// speed of light is a defined quantity.
//-----------------------------------------------------------------------

const DoubleWithUnit speed_of_light(299792458, 
				    units::m / units::s);

//-----------------------------------------------------------------------
/// Earth angular rotation frequency (1/sec)
//-----------------------------------------------------------------------

const DoubleWithUnit earth_rot_speed(360 / 86164.09054, 
				     units::deg / units::s);
				// 360 degrees / sidereal day in
				// seconds.

//-----------------------------------------------------------------------
/// 1 AU in meters. This comes from
/// http://neo.jpl.nasa.gov/glossary/au.html 
//-----------------------------------------------------------------------

const Unit AU("AU", 1.49597870691e11 * units::m);

//-----------------------------------------------------------------------
/// Planck's constant, in J s
//-----------------------------------------------------------------------

const DoubleWithUnit planck(6.6260693e-34, 
			    units::J * units::s);

  }
}
#endif

