#ifndef WGS84_CONSTANT_H
#define WGS84_CONSTANT_H
#include "double_with_unit.h"

namespace FullPhysics {
//-----------------------------------------------------------------------
/// This contains various constants that describe the wgs84 ellipsoid.
//-----------------------------------------------------------------------

  namespace OldConstant {
//-----------------------------------------------------------------------
/// Equatorial radius.
//-----------------------------------------------------------------------
const DoubleWithUnit wgs84_a(6378137.0000, units::m);

//-----------------------------------------------------------------------
/// Polar radius in meters.
//-----------------------------------------------------------------------
const DoubleWithUnit wgs84_b(6356752.3142450, units::m);

//-----------------------------------------------------------------------
/// Eccentricity squared. From CRC.
//-----------------------------------------------------------------------
// This isn't actually used anywhere, and the compiler gives an
// annoying warning about this. So comment out, until we happen to
// need this somewhere. 
// const double wgs84_esq = ((wgs84_a * wgs84_a - wgs84_b * wgs84_b) / 
// 			  (wgs84_a * wgs84_a)).value;
  }
}
#endif
