#ifndef CONSTANT_H
#define CONSTANT_H
#include "double_with_unit.h"
namespace FullPhysics {
/****************************************************************//**
This class contains various constants. We have this in a separate
class rather than just having a header of constants so that we can
replace these constants. In particular, some of the constants depend
on the planet we are working on (e.g., Venus).
*******************************************************************/
class Constant : public Printable<Constant> {
  public:
    virtual ~Constant() {}
    virtual void print(std::ostream& Os) const = 0;

    //-----------------------------------------------------------------------
    /// Rayleigh depolarization factor. 
    //-----------------------------------------------------------------------

    virtual double rayleigh_depolarization_factor() const = 0;
  
    //-----------------------------------------------------------------------
    /// Rayleigh "a" value. This along with "b" are the wavelength
    /// dependence coefficients for the refractive index.
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_a() const = 0;
  
    //-----------------------------------------------------------------------
    /// Rayleigh "b" value. This along with "a" are the wavelength
    /// dependence coefficients for the refractive index.
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_b() const = 0;

    //-----------------------------------------------------------------------
    /// Molar weight of dry air
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit molar_weight_dry_air() const = 0;

    //-----------------------------------------------------------------------
    /// Molar weight of water
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit molar_weight_water() const = 0;

    //-----------------------------------------------------------------------
    /// Avogadro constant
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit avogadro_constant() const = 0;

    //-----------------------------------------------------------------------
    /// The polynomial used to find the distance to the sun
    //-----------------------------------------------------------------------
    virtual boost::array<double, 7> solar_distance_param() const = 0;

  };
}
#endif
