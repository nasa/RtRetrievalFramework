#ifndef PRESSURE_LEVEL_INPUT_H
#define PRESSURE_LEVEL_INPUT_H
#include "printable.h"
#include "heritage_file.h"
#include "hdf_file.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  In a retrieval, there are typically two different pressure levels of
  interest. One is the pressure levels where various initial parameters
  are defined, e.g. Temperature read from an ECMWF file at specific
  pressure levels. The second set is the current pressure levels that
  define the layers used in the Radiative Transfer calculation.  The
  first set is fixed constant level, it is whatever was used when we
  initial read the input data. The second will potentially vary as we
  do a retrieval.

  This class captures the first, fixed set of pressure levels. This is
  little more than an Array of values.
*******************************************************************/

class PressureLevelInput : public Printable<PressureLevelInput> {
public:
//-----------------------------------------------------------------------
/// Constructor.
/// \param Press_level The pressure levels in ascending order. This
/// is in Pascals.
//-----------------------------------------------------------------------

  PressureLevelInput(const blitz::Array<double, 1>& Press_level)
    : press_level(Press_level.copy()) {}

//-----------------------------------------------------------------------
/// Constructor.
/// \param Heritage_file The HeritageFile to read. This takes the
/// column marked "Pressure".
//-----------------------------------------------------------------------

  PressureLevelInput(const HeritageFile& Heritage_file)
    : press_level(Heritage_file.data("Pressure")) {}
  
//-----------------------------------------------------------------------
/// Constructor.
/// \param Hdf_file The HdfFile to read. This reads the given HDF
/// group.
/// \param Hdf_group The HDF group to read.
//-----------------------------------------------------------------------

  PressureLevelInput(const HdfFile& Hdf_file, 
		     const std::string& Hdf_group = "Pressure")
    : press_level(Hdf_file.read_field<double, 1>(Hdf_group + "/Pressure")) {}

  virtual ~PressureLevelInput() {}

//-----------------------------------------------------------------------
/// Pressure levels that input data was defined on.
/// \return Pressure level in ascending order, in Pascals.
//-----------------------------------------------------------------------

  const blitz::Array<double, 1>& pressure_level() const {return press_level; }

//-----------------------------------------------------------------------
/// Print to a stream.
//-----------------------------------------------------------------------

  void print(std::ostream& Os) const;
private:
  blitz::Array<double, 1> press_level;
};
}
#endif
