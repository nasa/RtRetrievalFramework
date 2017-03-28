#ifndef LEVEL_1B_H
#define LEVEL_1B_H
#include "fp_time.h"
#include "printable.h"
#include "double_with_unit.h"
#include "spectral_range.h"
#include <blitz/array.h>
#include <stdint.h>

namespace FullPhysics {
/****************************************************************//**
  This is used to read a Level 1B file.
*******************************************************************/
class Level1b: public Printable<Level1b> {
public:
  virtual ~Level1b() { };

//-----------------------------------------------------------------------
/// Number of spectrometers.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// Latitude
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Latitude.
//-----------------------------------------------------------------------

  virtual DoubleWithUnit latitude(int i) const = 0;

//-----------------------------------------------------------------------
/// Longitude
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Longitude
//-----------------------------------------------------------------------

  virtual DoubleWithUnit longitude(int i) const = 0;

//-----------------------------------------------------------------------
/// Sounding zenith
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Sounding zenith
//-----------------------------------------------------------------------

  virtual DoubleWithUnit sounding_zenith(int i) const = 0;

//-----------------------------------------------------------------------
/// Sounding azimuth
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Sounding azimuth
//-----------------------------------------------------------------------

  virtual DoubleWithUnit sounding_azimuth(int i) const = 0;

//-----------------------------------------------------------------------
/// Return stokes coefficients
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Stokes coefficients, with size 4.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> stokes_coefficient(int i) const = 0;

//-----------------------------------------------------------------------
/// Solar zenith
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Solar zenith angle
//-----------------------------------------------------------------------

  virtual DoubleWithUnit solar_zenith(int i) const = 0;

//-----------------------------------------------------------------------
/// Solar azimuth
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Solar azimuth angle
//-----------------------------------------------------------------------

  virtual DoubleWithUnit solar_azimuth(int i) const = 0;

//-----------------------------------------------------------------------
/// Realtive azimuth
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Relative azimuth angle between solar and sounding azimuth
//     angles.
//-----------------------------------------------------------------------

  virtual DoubleWithUnit relative_azimuth(int i) const;

//-----------------------------------------------------------------------
/// Altitude
/// \param i Spectrometer index (between 0 and number_spectrometer() -
///    1)
/// \return Altitude
//-----------------------------------------------------------------------

  virtual DoubleWithUnit altitude(int i) const = 0;

//-----------------------------------------------------------------------
/// Relative velocity.
/// \return Relative velocity
//-----------------------------------------------------------------------

  virtual DoubleWithUnit relative_velocity(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Returns coefficients for an equation describing the special domain
/// used to translate radiance value indexes to their corresponding 
/// spectral grid. (ie wavenumber, wavelength, etc)
/// The meaning of these coefficients will be specific to the instrument
/// that measured the data.
//-----------------------------------------------------------------------

  virtual ArrayWithUnit<double, 1> spectral_coefficient(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Time of sounding.
//-----------------------------------------------------------------------

  virtual Time time(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Sounding ID. By convention this is a 64 bit integer.
//-----------------------------------------------------------------------

  virtual int64_t sounding_id() const = 0;

//-----------------------------------------------------------------------
/// Exposure index. By convention this is 1 based.
//-----------------------------------------------------------------------

  virtual int exposure_index() const = 0;

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const {Os << "Level1b";}

//-----------------------------------------------------------------------
/// Radiance, for given spectral band. This returns the radiance with
/// associated units. It may or may not have a uncertainity with the
/// radiance. 
//-----------------------------------------------------------------------

  virtual SpectralRange radiance(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Calculate an approximation to the size of the continuum signal
/// where there is no significant atmosphere absorption. We
/// approximate this by finding the 10 highest radiance values and
/// averaging them.
///
/// Optionally takes a list of sample indexes. Will only uses these
/// sample indexes for the calculation when supplied.
//-----------------------------------------------------------------------

  virtual DoubleWithUnit signal(int Spec_index, const std::vector<int>& Sample_indexes = std::vector<int>()) const;
   
};
} // End of FullPhysics namespace
#endif
