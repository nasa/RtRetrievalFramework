#ifndef LEVEL_1B_OCO1_H
#define LEVEL_1B_OCO1_H
#include "level_1b_hdf.h"
#include "hdf_sounding_id.h"
#include "fp_exception.h"
#include "hdf_file.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This reads a Level 1B file that is in the HDF format.
*******************************************************************/
class Level1bOco: public Level1bHdf {
public:

  Level1bOco(const std::string& Fname, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);
  Level1bOco(const boost::shared_ptr<HdfFile>& Hfile, 
	     const boost::shared_ptr<HdfSoundingId>& Sounding_id);

//-----------------------------------------------------------------------
/// The acquisition mode. The data we process will be "Nadir", "Glint"
/// or "Target". 
//-----------------------------------------------------------------------

  const std::string& acquisition_mode() const
  {
    return acquisition_mode_;
  }

//-----------------------------------------------------------------------
/// True if we have the SpikeEOF data available.
//-----------------------------------------------------------------------

  bool has_spike_eof(int Spec_index) const;
  blitz::Array<double, 1> spike_eof(int Spec_index) const;

//-----------------------------------------------------------------------
/// True if the Level 1 data has the solar relative velocity. Older 
/// versions of OCO and also the simulated data do not have this
/// field, so this returns false.
//-----------------------------------------------------------------------

  bool has_solar_relative_velocity() const
  {
    return has_solar_relative_velocity_;
  }

//-----------------------------------------------------------------------
/// Return land fraction, as a percentage going from 0 to 100.
//-----------------------------------------------------------------------

  double land_fraction() const {return land_fraction_;}

//-----------------------------------------------------------------------
/// Distance from sounding location to the sun.
//-----------------------------------------------------------------------

  DoubleWithUnit solar_distance() const
  {
    if(!has_solar_relative_velocity())
      throw Exception("Level 1b file does not have the solar relative velocity");
    return solar_distance_;
  }

//-----------------------------------------------------------------------
/// Velocity of sun along the sounding location/sun vector.
/// Note this includes the rotation of the earth along with the
/// earth/sun velocity, so we can take this directly to calculate the
/// doppler shift.
//-----------------------------------------------------------------------

  DoubleWithUnit solar_velocity() const
  {
    if(!has_solar_relative_velocity())
      throw Exception("Level 1b file does not have the solar relative velocity");
    return solar_velocity_;
  }
  
protected:
  // For sub classes
  Level1bOco() {}

  virtual SpectralRange radiance_no_uncertainty(int Spec_index) const;

  DoubleWithUnit solar_distance_, solar_velocity_;
  double land_fraction_;
  bool has_solar_relative_velocity_;
  std::string acquisition_mode_;

private:
  void initialize();
};
}

#endif
