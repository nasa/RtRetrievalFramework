#ifndef INSTRUMENT_H
#define INSTRUMENT_H
#include "printable.h"
#include "spectral_window.h"
#include "array_ad.h"
#include "state_vector.h"
#include "spectrum.h"
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace FullPhysics {
/****************************************************************//**
  This applies a instrument model to radiances.
*******************************************************************/

class Instrument : virtual public StateVectorObserver, 
		   public Observable<Instrument> {
public:
  virtual ~Instrument() {}

  virtual void add_observer(Observer<Instrument>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Instrument>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Clone an Instrument object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Instrument>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Instrument object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Instrument> clone() const = 0;

//-----------------------------------------------------------------------
/// Give number of spectrometers.
//-----------------------------------------------------------------------

  virtual int number_spectrometer() const = 0;

//-----------------------------------------------------------------------
/// Apply the instrument model to both the radiance and derivatives.
///
/// \param High_resolution_spectrum High resolution spectrum.
/// \param Pixel_list List of pixels to include in radiance
/// \param Spec_index Spectral index
/// \return Spectrum with instrument model applied.
//-----------------------------------------------------------------------

  virtual Spectrum apply_instrument_model(
    const Spectrum& High_resolution_spectrum,
    const std::vector<int>& Pixel_list,
    int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// This is the pixel wavenumber/wavelength for each pixel.
//-----------------------------------------------------------------------

  virtual SpectralDomain pixel_spectral_domain(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// Band name for given Spec_index.
//-----------------------------------------------------------------------

  virtual std::string band_name(int Spec_index) const = 0;

//-----------------------------------------------------------------------
/// In general, the name used in HDF files for a particular band is
/// similar but not identical to the more human readable band_name.
/// For example, with GOSAT we use the HDF field name "weak_co2", but
/// the band name is "WC-Band". This gives the HDF name to use.
///
/// The default implementation just returns the same string as the
/// band name.
//-----------------------------------------------------------------------

  virtual std::string hdf_band_name(int Spec_index) const
  { return band_name(Spec_index); }

//-----------------------------------------------------------------------
/// This is the half width of the ILS in wavenumber.
//-----------------------------------------------------------------------

  virtual DoubleWithUnit
  ils_half_width(int Spec_index) const = 0;
  virtual void ils_half_width(int Spec_index, DoubleWithUnit& half_width) = 0;
  virtual void print(std::ostream& Os) const {Os << "Instrument";}
};
}
#endif
