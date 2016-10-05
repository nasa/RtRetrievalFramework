#ifndef DEFAULT_CONSTANT_H
#define DEFAULT_CONSTANT_H
#include "constant.h"
namespace FullPhysics {
/****************************************************************//**
This class is an implementation of Constant that uses
hard coded values suitable for Earth. This is useful for things like
unit tests where we just need reasonable values.
*******************************************************************/
class DefaultConstant : public Constant {
  public:
    DefaultConstant() {}
    virtual ~DefaultConstant() {}
    virtual void print(std::ostream& Os) const 
    { Os << "Default constants for Earth"; }

    //-----------------------------------------------------------------------
    /// Rayleigh depolarization factor. 
    //-----------------------------------------------------------------------

    virtual double rayleigh_depolarization_factor() const 
    { return 0.02790; // depolarization factor; for air is 0.02790 (young, 1980)
    }
    //-----------------------------------------------------------------------
    /// Rayleigh "a" value. This along with "b" are the wavelength
    /// dependence coefficients for the refractive index.
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_a() const
    {
      // a and b are wavelength dependence coefficients for the
      // refractive index (allen, 1964) (note: wavelengths must be
      // in microns) For Earth, a = 2.871e-04, b = 5.67e-03
      return DoubleWithUnit(2.871e-04, "dimensionless");
    }
  
    //-----------------------------------------------------------------------
    /// Rayleigh "b" value. This along with "a" are the wavelength
    /// dependence coefficients for the refractive index.
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_b() const 
    {
      // a and b are wavelength dependence coefficients for the
      // refractive index (allen, 1964) (note: wavelengths must be
      // in microns) For Earth, a = 2.871e-04, b = 5.67e-03
      return DoubleWithUnit(5.67e-03, "dimensionless");
    }

    //-----------------------------------------------------------------------
    /// Molar weight of dry air
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit molar_weight_dry_air() const
    {
      return DoubleWithUnit(28.9644, "g / mol");
    }

    //-----------------------------------------------------------------------
    /// Molar weight of water
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit molar_weight_water() const
    {
      return DoubleWithUnit(18.01528, "g / mol");
    }

    //-----------------------------------------------------------------------
    /// Avogadro constant
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit avogadro_constant() const
    {
      return DoubleWithUnit(6.02214179e23, "mol^-1");
    }
    
    //-----------------------------------------------------------------------
    /// The polynomial used to find the distance to the sun
    //-----------------------------------------------------------------------
    virtual boost::array<double, 7> solar_distance_param() const
    {
      const boost::array<double, 7> sdp =
	{{0.98334, -1.82823e-5, 2.30179e-6, 6.62402e-9,
	  -1.33287e-10, 3.98445e-13, -3.54239e-16}};
      return sdp;
    }
  };
}
#endif
