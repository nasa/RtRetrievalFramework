#ifndef ILS_H
#define ILS_H
#include "state_vector.h"
#include "double_with_unit.h"
#include "spectral_domain.h"

namespace FullPhysics {
/****************************************************************//**
  This class models an Instrument Line Shape (ILS). We convolve with
  high resolution data to produce a model of what we expect to observe
  in the Level 1b data.
*******************************************************************/
class Ils : virtual public StateVectorObserver,
		    public Observable<Ils> {
public:
  virtual ~Ils() {}

  virtual void add_observer(Observer<Ils>& Obs) 
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<Ils>& Obs) 
  { remove_observer_do(Obs, *this);}

//-----------------------------------------------------------------------
/// Apply the ILS. This version does not calculate a Jacobian.
///
/// \param High_resolution_wave_number The wave numbers going with the
///   high resolution radiance data. This is in cm^-1, and should be
///   ordered from smallest to largest wavenumber.
/// \param High_resolution_radiance The high resolution radiance
///   data. This is in w/m^2 / sr / cm^-1
/// \param Pixel_list List of instrument pixels to include in the
///   results. The order of the pixels is the same order that we
///   return our results in.
/// \return Radiance with ILS applied. This is in w/m^2 / sr / cm^-1.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const blitz::Array<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;

//-----------------------------------------------------------------------
/// Apply the ILS. This includes propagating the Jacobian from the
/// high resolution data, and adding in any dependence of the ILS on
/// the state vector elements (e.g., dispersion state vector elements).
///
/// \param High_resolution_wave_number The wave numbers going with the
///   high resolution radiance data. This is in cm^-1, and should be
///   ordered from smallest to largest wavenumber.
/// \param High_resolution_radiance The high resolution radiance
///   data and jacobian . This is in w/m^2 / sr / cm^-1
/// \param Pixel_list List of instrument pixels to include in the
///   results. The order of the pixels is the same order that we
///   return our results in.
/// \return Radiance with ILS applied, and Jacobian This is in w/m^2 /
///   sr / cm^-1. 
//-----------------------------------------------------------------------

  virtual ArrayAd<double, 1> apply_ils
  (const blitz::Array<double, 1>& High_resolution_wave_number,
   const ArrayAd<double, 1>& High_resolution_radiance,
   const std::vector<int>& Pixel_list) const = 0;

//-----------------------------------------------------------------------
/// Descriptive name of the band.
//-----------------------------------------------------------------------

  virtual std::string band_name() const = 0;

//-----------------------------------------------------------------------
/// Clone an Ils object. Note that the cloned version will *not*
/// be attached to and StateVector or Observer<Ils>, although you
/// can of course attach them after receiving the cloned object.
///
/// Because this isn't attached to the StateVector, one use of the
/// clone operator is to create a "frozen" Ils object.
//-----------------------------------------------------------------------

  virtual boost::shared_ptr<Ils> clone() const = 0;

//-----------------------------------------------------------------------
/// In general, the name used in HDF files for a particular band is
/// similar but not identical to the more human readable band_name.
/// For example, with GOSAT we use the HDF field name "weak_co2", but
/// the band name is "WC-Band". This gives the HDF name to use.
///
/// The default implementation just returns the same string as the
/// band name.
//-----------------------------------------------------------------------

  virtual std::string hdf_band_name() const { return band_name();}

//-----------------------------------------------------------------------
/// This is the pixel grid for each pixel.
//-----------------------------------------------------------------------

  virtual SpectralDomain pixel_grid() const = 0;

//-----------------------------------------------------------------------
/// This is the half width of the ILS.
//-----------------------------------------------------------------------

  virtual DoubleWithUnit ils_half_width() const = 0;

//-----------------------------------------------------------------------
/// Set the half width of the ILS.
//-----------------------------------------------------------------------

  virtual void ils_half_width(const DoubleWithUnit& half_width) = 0;

};
}
#endif
