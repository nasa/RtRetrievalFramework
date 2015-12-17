#ifndef HDF_CONSTANT_H
#define HDF_CONSTANT_H
#include "constant.h"
#include "hdf_file.h"
namespace FullPhysics {
  /******************************************************************
  This class is an implementation of Constant where all the values are
  read from an Hdf file.
  *******************************************************************/
  class HdfConstant : public Constant {
  public:
    HdfConstant(const boost::shared_ptr<HdfFile>& Hdf_file);
    virtual ~HdfConstant() {}
    virtual void print(std::ostream& Os) const 
    { Os << "HdfConstant"; }

    //-----------------------------------------------------------------------
    /// Rayleigh depolarization factor. 
    /// \todo is this really dimensionless?
    //-----------------------------------------------------------------------

    virtual double rayleigh_depolarization_factor() const 
    { return -1;
    }
    //-----------------------------------------------------------------------
    /// Rayleigh "a" value. This along with "b" are the wavelength
    /// dependence coefficients for the refractive index.
    /// \todo what are the units here?
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_a() const
    {
      return DoubleWithUnit(-1, "dimensionless");
    }
  
    //-----------------------------------------------------------------------
    /// Rayleigh "b" value. This along with "a" are the wavelength
    /// dependence coefficients for the refractive index.
    /// \todo what are the units here?
    //-----------------------------------------------------------------------

    virtual DoubleWithUnit rayleigh_b() const 
    {
      return DoubleWithUnit(-1, "dimensionless");
    }
    //-----------------------------------------------------------------------
    /// Molar weight of dry air
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit molar_weight_dry_air() const
    {
      return DoubleWithUnit(-1, "g / mol");
    }

    //-----------------------------------------------------------------------
    /// Molar weight of water
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit molar_weight_water() const
    {
      return DoubleWithUnit(-1, "g / mol");
    }

    //-----------------------------------------------------------------------
    /// Avogadro constant
    //-----------------------------------------------------------------------
    virtual DoubleWithUnit avogadro_constant() const
    { return DoubleWithUnit(-1, "mol^-1"); }

    //-----------------------------------------------------------------------
    /// The polynomial used to find the distance to the sun
    //-----------------------------------------------------------------------
    virtual boost::array<double, 7> solar_distance_param() const
    {
      const boost::array<double, 7> sdp =
	{{-1, -1, -1, -1, -1, -1, -1}};
      return sdp;
    }
  };
}
#endif
